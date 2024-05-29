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

    cerr << "Args:" << endl;
    cerr << "  fpass:      " << fpass << endl;
    cerr << "  rminfo:     " << rminfo << endl;
    cerr << "  keepaa:     " << keepaa << endl;
    cerr << "  macfilter:  " << macfilter << endl;
    cerr << "  aafilter:   " << aafilter << endl;

    size_t len = BUFSIZE;
    const size_t lenstart = len;
    char* line = (char*) malloc(len*sizeof(char));
    size_t nline = 0;

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
        while(getline(&line, &len, stdin) != -1) {

            // genomic position
            char* pos = line;
            char* posend = strchr(pos, '\t'); // end of genomic position (exclusive, points to tab char)
            if (posend == NULL) // no tab char -> invalid line (can happen for the last line when it contains only a newline character, or when concatenating files the header of the new file) -> skip
                continue;
    //        *posend = '\0'; // null terminate the pos string

            // variant ID
            char* varid = posend+1;
            char* varidend = strchr(varid, '\t');

            // alleles
            char* refall = varidend+1;
            char* refallend = strchr(refall, '\t'); // end of first allele (exclusive)
            char* altall = refallend+1;
            char* altallend = strchr(refallend+1, '\t'); // end of alternative alleles (exclusive)

            // count number of alternative alles:
            size_t nalt = 0;
            *altallend = '\0'; // null terminate the allele string (to stop following loop)
            for (char* t = altall; t != NULL; t = strchr(t+1, ',')) // will stop at altallend as we have null terminated the buffer there
                nalt++;
            *altallend = '\t'; // restore tab

            // initialize allele counters
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
                if (strcmp(filter, "PASS") != 0) // not passed!
                    continue; // skip this line
            }

            // info
            char* info = filterend+1; // start
            char* infoend = strchr(info, '\t'); // end
            *infoend = '\0'; // null terminate INFO field

            // AAScore filter
            if (aafilter > 0) {
                bool pass = false;
                char* aa = findInfoField(info, "AAScore=");
                aa += 8; // forward to beginning of first value

                // iterate over all alt alleles
                for (size_t n = 0; n < nalt; n++) {
                    char* aaend;
                    if (n < nalt-1) { // more than one alt allele
                        aaend = strchr(aa, ','); // will be found in a proper VCF
                    } else { // last alt allele
                        aaend = strchr(aa, ';'); // end of AAScore field or null if at the end of INFO
                        if (aaend == NULL)
                            aaend = infoend;
                    }
                    float aav = stof(string(aa, aaend - aa));
                    if (aav >= aafilter) {
                        pass = true;
                        break;
                    }
                    aa = aaend + 1; // beginning of next value
                }
                if (!pass)
                    continue; // skip this line
            }

            // genotypes
            char* gtstart = infoend+1; // start of genotypes (pointing at first gt char!)
            int gtflag = 1; // signalizes if the current field contains genotypes or not
            for (char* gt = gtstart; *gt != '\0'; gt++) { // until the end of the line buffer

                if (gtflag && *gt >= '0' && *gt <= '9') { // points to valid haplotype
                    an++; // increase allele number
                    if (*gt >= '1') {
                        int idx = *gt - '0';
                        while (*(gt+1) >= '0' && *(gt+1) <= '9') { // continue digit by digit
                            gt++;
                            idx = idx * 10 + (*gt - '0');
                        }
                        ac[idx-1]++; // increase corresponding alt allele counter
                    }
                } else if (*gt == ':') { // double colon marks the end of a genotype
                    gtflag = 0;
                } else if (*gt == '\t') { // tab marks the end of a gt field
                    gtflag = 1;
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
                        break;
                    }
                }
                if (!pass)
                    continue; // skip this line
            }

            // print VCF line

            // CHROM (chromosome name)
            cout << chrom << "\t";

            // POS, ID, alleles, QUAL, FILTER
            cout << pos;

            // INFO
            // print self generated values
            float anf = (float) an;
            printf("\tAF=%.8f", ac[0]/anf); // AF of first alt allele
            for (size_t n = 1; n < nalt; n++)
                printf(",%.8f", ac[n]/anf); // AF of further alleles if multi-allelic
            printf(";AC=%lu", ac[0]); // AC of first alt allele
            for (size_t n = 1; n < nalt; n++)
                printf(",%lu", ac[n]); // AC of further alleles if multi-allelic
            printf(";AN=%lu", an); // AN

            // original values
            if (!rminfo) { // take over all original values
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
                char* aa = findInfoField(info, "AAScore=");
                if (aa != NULL) { // found AAScore
                    char* aaend = strchr(aa, ';'); // will be found or is already at the end
                    if (aaend != NULL)
                        *aaend = '\0'; // null terminate AAScore
                    cout << ";" << aa;
                }
            }

            // FORMAT
            if (parsegq)
                cout << "\tGT:GQ\t";
            else
                cout << "\tGT\t";

            // genotypes
            cout << gtstart; // rest of the line buffer containing the genotypes, ends with newline!

            nline++;
        }

        cout << flush;
        free(ac);

    } // END contains data

    cerr << "Number of variants: " << nline << endl;
    cerr << "Line buffer size: " << len;
    if (len != lenstart)
        cerr << " -> changed!!";
    cerr << endl;

    free(line);

}

