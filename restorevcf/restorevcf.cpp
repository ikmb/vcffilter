/*
 * restorevcf.c
 *
 *  Created on: May 15, 2024
 *      Author: lwienbrandt
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "RestoreArgs.h"

#define BUFSIZE 1073741824

using namespace std;

char* findInfoField(char* info, const char* query) {
    char* ret = NULL;
    do {
        ret = strstr(info, query);
        if (ret != NULL && (ret == info || *(ret-1) == ';'))
            return ret;
    } while(ret != NULL);
    return NULL;
}

int ptrcmp (const void* a, const void* b) {
   return ( *(char**)a - *(char**)b );
}

int main (int argc, char **argv) {

    // parse args
    RestoreArgs args = RestoreArgs::parseArgs(argc, argv);

    bool fpass = args.fpass;
    bool rminfo = args.rminfo;
    bool keepaa = args.keepaa;
    size_t macfilter = args.macfilter;
    float aafilter = args.aafilter;
    float missfilter = args.missfilter;
    bool filterunk = args.filterunk;
    bool splitma = args.splitma;

    cerr << "Args:" << endl;
    cerr << "  fpass:         " << fpass << endl;
    cerr << "  rminfo:        " << rminfo << endl;
    cerr << "  keepaa:        " << keepaa << endl;
    cerr << "  macfilter:     " << macfilter << endl;
    cerr << "  aafilter:      " << aafilter << endl;
    cerr << "  missfilter:    " << missfilter << endl;
    cerr << "  filterunknown: " << filterunk << endl;
    cerr << "  splitma:       " << splitma << endl;

    size_t len = BUFSIZE;
    const size_t lenstart = len;
    char* line = (char*) malloc(len*sizeof(char));
    size_t nread = 0;
    size_t nprint = 0;
    size_t nskip = 0;
    size_t nsplit = 0;

    // parse header
    size_t nh = getline(&line, &len, stdin); // read header line
    if (nh > 0 && nh != (size_t)-1) { // contains data

        // overwrite newline character at the end of the line (prevents correct parsing below)
        *(line+nh-1) = '\0';

        // search for more args (separated by tab)
        char* chrend = strchr(line, '\t');
        // null terminate chromosome name
        if (chrend != NULL) {
            *chrend = '\0';
        }
        string chrom(line); // copy chromosome name

        // parse more args
        bool parsegq = false;
        char* arg = (chrend != NULL) ? chrend+1 : NULL;
        while (arg != NULL) {
            char* argend = strchr(arg, ';'); // separator for the args in the header is the semicolon char
            if (argend != NULL) {
                *argend = '\0';
            }
            if (strcmp(arg, "--gq") == 0) // --gq option was set -> restore GQ field
                parsegq = true;
            arg = (argend != NULL) ? argend+1 : NULL;
        }

        // reserve space for allele counters
        size_t nac = 10;
        size_t* ac = (size_t*) malloc(nac * sizeof(size_t));

        // parse rest of file
        ssize_t nline = 0;
        while((nline = getline(&line, &len, stdin)) != -1) {

            // genomic position
            char* pos = line;
            char* posend = strchr(pos, '\t'); // end of genomic position (exclusive, points to tab char)
            if (posend == NULL) // no tab char -> invalid line (can happen for the last line when it contains only a newline character, or when concatenating files the header of the new file) -> skip
                continue;
            nread++;

            // variant ID
            char* varid = posend+1;
            char* varidend = strchr(varid, '\t');

            // alleles
            char* refall = varidend+1;
            char* refallend = strchr(refall, '\t'); // end of first allele (exclusive)
            char* altall = refallend+1;
            char* altallend = strchr(refallend+1, '\t'); // end of alternative alleles (exclusive)

            // count number of alternative alles + check if unknown + prepare split
            size_t nalt = 0;
            int unkidx = -1; // points to the alt allele which is unknown (-1 if none)
            vector<char*> maaltalleles; // used to store the pointers to the alt alleles
            *altallend = '\0'; // null terminate the allele string (to stop following loop)
            for (char* t = refallend; t != NULL; t = strchr(t+1, ',')) { // will stop at altallend as we have null terminated the buffer there
                // t always points to the delimiter before the current allele (also at the beginning)!
                if (filterunk && *(t+1) == '*') { // found unknown allele
                    unkidx = nalt;
                }
                if (splitma && nalt >= 1) { // we have a multi-allelic variant here which we want to split
                    if (nalt == 1) { // second ma alt allele
                        // need to insert the first alt allele in our vector
                        maaltalleles.push_back(altall);
                    }
                    maaltalleles.push_back(t+1); // store beginning of this alternative allele
                    *t = '\0'; // need to null terminate the delimiter for printing later
                }
                nalt++;
            }

            // filter single unknown "*" alleles
            if (unkidx >= 0 && nalt == 1) { // filter is activated and an unknown allele was found which is the only alt allele
                nskip++;
                continue; // skip this line
            } // else if there was an unknown allele together in a multi-allelic context, we filter it later

            // prepare for ma splits
            bool masplitnow = splitma && nalt > 1;
            vector<bool> maaltfilter;
            if (masplitnow) { // multi-allelic variant that needs to be split
                nsplit++;
                *refallend = '\0'; // null terminate the ref allele field as we need to separate printing of alt alleles later
                // initialize the flags for which allele should be filtered -> default: not filtered
                maaltfilter.assign(nalt, false);
                if (unkidx >= 0) // if there's an unknown allele to be filtered, mark it already
                    maaltfilter[unkidx] = true;
            } else { // no ma splitting required
                *altallend = '\t'; // restore tab -> faster printing later
            }

            // initialize allele counters
            // (we use directly allocated memory here as it is significantly faster than a vector!)
            if (nalt > nac) {
                ac = (size_t*) realloc(ac, nalt * sizeof(size_t));
                nac = nalt;
            }
            for (size_t i = 0; i < nalt; i++)
                ac[i] = 0;
            size_t an = 0;

            // qual field
            char* qual = altallend+1; // start
            char* qualend = strchr(qual, '\t'); // end (exclusive)

            // filter field
            char* filter = qualend+1; // start of filter field
            char* filterend = strchr(filter, '\t'); // end of filter (exclusive)
            *filterend = '\0'; // null terminate filter field, as we are going to modify the following INFO fields for printing

            // FILTER == PASS filter
            if (fpass) {
                if (strcmp(filter, "PASS") != 0) { // not passed!
                    nskip += masplitnow ? nalt : 1;
                    continue; // skip this line
                }
            }

            // info
            char* info = filterend+1; // start
            char* infoend = strchr(info, '\t'); // end
            *infoend = '\0'; // null terminate INFO field

            // AAScore filter
            char* aa = NULL; // will be beginning of AAScore field, if we filter by AAScore or keep AAScores while removing the rest of INFO (as for split MA's)
            vector<char*> aaval; // pointers to the AAScore values
            if (aafilter > 0 || keepaa) {
                bool pass = !(aafilter > 0); // only required as false if we want to filter by AAScore
                aa = findInfoField(info, "AAScore=");
                char* aatmp = aa+8; // beginning of first value

                // iterate over all alt alleles
                for (size_t n = 0; n < nalt; n++) {
                    char* aaend;
                    if (n < nalt-1) { // more than one alt allele
                        aaend = strchr(aatmp, ','); // will be found in a proper VCF
                    } else { // last alt allele
                        aaend = strchr(aatmp, ';'); // end of AAScore field or null if at the end of INFO
                        if (aaend == NULL)
                            aaend = infoend;
                    }
                    if (keepaa){
                        aaval.push_back(aatmp);
                        *aaend = '\0'; // null terminate for printing separately
                    }
                    float aav = stof(string(aatmp, aaend - aatmp));
                    if (aav >= aafilter) { // value above threshold
                        pass = true;
                        if (!keepaa && !masplitnow) // if we do not keepaa explicitly and if we do not split MA's, we can stop here
                            break;
                    } else if (masplitnow && aafilter > 0) { // value below threshold and we are splitting MA here and a filter is active
                        // mark this allele to be filtered later
                        maaltfilter[n] = true;
                    }
                    aatmp = aaend + 1; // beginning of next value
                }
                if (!pass) { // all AAScores are below the threshold
                    nskip += masplitnow ? nalt : 1;
                    continue; // skip this line
                }
            }

            // genotypes
            char* gtstart = infoend+1; // start of genotypes (pointing at first gt char!)
            vector<char*> magts; // for MA splits, the beginning of the genotypes
            vector<vector<char*>> magtparts;  // for MA splits and large allele indices we need to replace more than one digit. we replace them with '\0' with this vector pointing to all replaced positions plus one (so the next part to print)
            if (masplitnow) {
                // copy gts for the multi-allelic splits (including null terminator!)
                size_t gtsize = (line + nline) - gtstart + 1;
                magts.resize(nalt, NULL); // reserve for all alt alleles
                magts[0] = gtstart; // for the first alternative allele, we keep the line buffer
                for (size_t a = 1; a < nalt; a++) { // for all other alt alleles, we reserve space and copy the gts, if we do not filter them anyway
                    if (!maaltfilter[a]) {
                        magts[a] = (char*) malloc(gtsize * sizeof(char));
                        memcpy(magts[a], gtstart, gtsize * sizeof(char));
                    }
                }
                if (nalt >= 10) {
                    magtparts.resize(nalt); // we need one for each alt allele!
                }
            }
            int gtflag = 1; // signalizes if the current field contains genotypes or not
            size_t ngtmiss = 0; // missing alleles counter
            for (char* gt = gtstart; *gt != '\0'; gt++) { // until the end of the line buffer

                if (gtflag && *gt >= '0' && *gt <= '9') { // points to valid haplotype
                    an++; // increase allele number
                    if (*gt >= '1') {
                        size_t gtpos = gt - gtstart;
                        size_t idx = *gt - '0';
                        while (*(gt+1) >= '0' && *(gt+1) <= '9') { // continue digit by digit
                            gt++;
                            idx = idx * 10 + (*gt - '0');
                            if (masplitnow) {
                                // allele index is >= 10, so we need to delete all digits after the first (replace with \0, but we store the following position for printing later)
                                // we take care of the first digit later
                                *gt = '\0';
                                magtparts[0].push_back(gt+1);
                                size_t pos = gt - gtstart;
                                for (size_t a = 1; a < nalt; a++) {
                                    if (!maaltfilter[a]) {
                                        *(magts[a]+pos) = '\0';
                                        magtparts[a].push_back(magts[a]+pos+1);
                                    }
                                }
                            }
                        }
                        ac[idx-1]++; // increase corresponding alt allele counter
                        if (masplitnow) {
                            // set current allele to '1' and all others to '0'
                            for (size_t a = 0; a < nalt; a++) {
                                if (!maaltfilter[a]) {
                                    if (a == idx-1)
                                        *(magts[a]+gtpos) = '1';
                                    else
                                        *(magts[a]+gtpos) = '0';
                                }
                            }
                        }
                    }
                } else if (gtflag && *gt == '.') { // missing gt
                    ngtmiss++;
                } else if (*gt == ':') { // double colon marks the end of a genotype
                    gtflag = 0;
                } else if (*gt == '\t') { // tab marks the end of a gt field
                    gtflag = 1;
                }

            }

            // GT missingness filter
            if (missfilter > 0) {
                if (ngtmiss / (float) an >= missfilter) {
                    // skip this line
                    if (!masplitnow) {
                        nskip++;
                        continue;
                    } else // need to continue here for a proper cleanup and we set all filters to '1'
                        maaltfilter.assign(nalt, true);
                }
            }

            // MAC filter
            if (macfilter) {
                bool pass = false;
                size_t mac = 0;
                // iterate over all alt alleles
                for (size_t n = 0; n < nalt; n++) {
                    // check minor(!) allele count
                    mac = (ac[n] <= an/2) ? ac[n] : an - ac[n];
                    if (mac >= (size_t) macfilter) {
                        pass = true;
                        if (!masplitnow) // if we are not splitting MA's here, we can stop
                            break;
                    } else if (masplitnow) { // value below threshold and we are splitting MA here
                        // mark this allele to be filtered later
                        maaltfilter[n] = true;
                    }
                }
                if (!pass) { // all MAC's are below the threshold
                    // skip this line
                    if (!masplitnow) {
                        nskip++;
                        continue;
                    } // else need to continue here for a proper cleanup, all filters were already set to '1'
                }
            }

            // *****************
            // print VCF line(s)
            // *****************

            size_t a = 0;
            do { // for each alt allele, if we split an MA, or only once if not

                if (!masplitnow || !maaltfilter[a]) { // only, if we do not filter this variant

                    // CHROM (chromosome name)
                    cout << chrom << "\t";

                    // POS, ID, ref allele + alt alleles, QUAL, FILTER if not splitting
                    cout << pos;

                    if (masplitnow) { // need to print correct alt allele now and then QUAL and FILTER
                        cout << "\t" << maaltalleles[a];
                        cout << "\t" << qual; // prints QUAL + FILTER
                    }

                    // INFO
                    // print self generated values
                    float anf = (float) an;
                    printf("\tAF=%.8f", ac[a]/anf); // AF of first alt allele
                    if (!masplitnow) {
                        for (size_t n = 1; n < nalt; n++)
                            printf(",%.8f", ac[n]/anf); // AF of further alleles if multi-allelic
                    }
                    printf(";AC=%lu", ac[a]); // AC of first alt allele
                    if (!masplitnow) {
                        for (size_t n = 1; n < nalt; n++)
                            printf(",%lu", ac[n]); // AC of further alleles if multi-allelic
                    }
                    printf(";AN=%lu", an); // AN

                    // original values
                    if (!rminfo) { // take over all original values -> no MA split possible here
                        if (info != infoend) { // info is not empty
                            cout << ";";

                            // find original values of the replaced ones above -> prefix them with "Org"
                            char* org[3];
                            org[0] = findInfoField(info, "AF=");
                            org[1] = findInfoField(info, "AC=");
                            org[2] = findInfoField(info, "AN=");
                            // sort...
                            qsort(org, 3, sizeof(char*), ptrcmp);

                            // print INFO and prefix above fields with 'Org'
                            char* infoit = info;
                            for (int i = 0; i < 3; i++) {
                                if (org[i] != NULL) {
                                    *(org[i]) = '\0'; // temporarily add null terminator (luckily, we know that we replace an 'A' here...)
                                    cout << infoit << "OrgA"; // we add the 'A' here, as we replaced it above with the null terminator...
                                    infoit = org[i]+1; // next char after the null terminator
                                }
                            }
                            // print remainder
                            cout << infoit;
                        }
                    } else if (keepaa) { // remove all original, but keep AAScore
                        if (aa != NULL) { // AAScore is present (already detected before)
                            cout << ";AAScore=" << aaval[a]; // print AAScore for current allele
                            if (!masplitnow) { // print all values for the other alleles as well
                                for (size_t n = 1; n < nalt; n++)
                                    cout << "," << aaval[n];
                            }
                        }
                    }

                    // FORMAT
                    if (parsegq)
                        cout << "\tGT:GQ\t";
                    else
                        cout << "\tGT\t";

                    // genotypes (all buffers end with newline!)
                    if (masplitnow) {
                        cout << magts[a];
                        if (nalt >= 10) { // there were indices >= 10 which have been replaced by \0
                            // print remaining parts after the replacement with \0
                            for (char* p : magtparts[a])
                                cout << p;
                        }
                    } else
                        cout << gtstart;

                    nprint++;

                } else { // this variant is filtered from a multi-allelic split
                    nskip++;
                }
                a++;

            } while(masplitnow && a < nalt); // for each alt allele, if we split an MA, or only once if not

            // cleanup for MA splits
            if (masplitnow) {
                for (size_t a = 1; a < nalt; a++) {
                    if (magts[a])
                        free(magts[a]);
                }
            }

        }

        cout << flush;
        free(ac);

    } // END contains data

    cerr << "Number of read variants: " << nread << endl;
    cerr << "Number of printed variants: " << nprint << endl;
    cerr << "Number of splitted variants: " << nsplit << endl;
    cerr << "Number of skipped variants (after split): " << nskip << endl;
    cerr << "Line buffer size: " << len;
    if (len != lenstart)
        cerr << " -> changed!!";
    cerr << endl;

    free(line);

}

