/*
 *    Copyright (C) 2025 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *
 *    This file is part of Vcffilter.
 *
 *    Vcffilter is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    Vcffilter is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Vcffilter. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "RemoveArgs.h"

// large buffer
#define BUFSIZE 1073741824

using namespace std;

int main (int argc, char **argv) {

    // parse args
    RemoveArgs args = RemoveArgs::parseArgs(argc, argv);

    string skipidfilename = args.skipidfilename;
    size_t macfilter = args.macfilter;
    float maffilter = args.maffilter;
    float missfilter = args.missfilter;

    cerr << "Args:" << endl;
    cerr << "  skip ID file:  " << skipidfilename << endl;
    cerr << "  macfilter:     " << macfilter << endl;
    cerr << "  maffilter:     " << maffilter << endl;
    cerr << "  missfilter:    " << missfilter << endl;
    cerr << endl;


    FILE* skipidfile = fopen(args.skipidfilename.c_str(), "r");

    // print usage
    if (!skipidfile) {
        cerr << " Unable to open file " << skipidfilename << endl;
        cerr << " Please provide a valid file with sample IDs to be removed.\n";
        exit(EXIT_FAILURE);
    }

    size_t len = BUFSIZE;
    char* line = (char*) malloc(len*sizeof(char));
    vector<string> skipids;
    vector<size_t> skipidxs;

    size_t nvars = 0;
    size_t nsamples = 0;
    size_t nskip = 0;

    size_t nh = getline(&line, &len, stdin); // read first line
    if (nh > 0 && nh != (size_t)-1) { // contains data

        // copy header until #CHROM line
        do {
            if (string(line,6).compare("#CHROM") == 0) { // found the #CHROM line (last line of header)
                break;
            }
            cout << line;

        } while((nh = getline(&line, &len, stdin)) != (size_t)-1);

        // add header line indicating the use of this tool
        cout << "##removesamples_command=";
        for (int i = 0; i < argc; i++)
            cout << argv[i] << " ";
        cout << endl;

        // load sample IDs that should be excluded from input file
        {
            cerr << "Reading sample IDs from " << argv[1] << endl;
            size_t sidlen = 128;
            char* sid = (char*) malloc(sidlen*sizeof(char));
            size_t nc;
            while((nc = getline(&sid, &sidlen, skipidfile)) != (size_t)-1) {
                skipids.push_back(string(sid, nc-1)); // discard newline char
            }
            fclose(skipidfile);
            free(sid);
            cerr << "Read " << skipids.size() << " sample IDs." << endl;
        }

        // parse the #CHROM line for all sample IDs
        {
            cerr << "Will remove the following samples found in the input stream:" << endl;

            // skip the first 9 tab characters (standard columns until first sample column)
            char* s = line;
            for (int t = 0; t < 9; t++) {
                s = strchr(s+1, '\t');
            }
            // print standard columns
            if (s)
                cout << string(line, s-line); // excluding tab character at end
            else // input contains no samples
                cout << string(line, string(line).size()-1); // discard newline at end

            // parse all sample IDs and check whether they should be excluded
            size_t curridx = 0;
            char* se = s; // s either points to tab or is null
            while (se) {
                s++;
                se = strchr(s, '\t');
                string sampleid;
                if (se) {
                    sampleid = string(s, se-s); // sample ID string excluding tab
                } else { // tab not found -> last sample
                    sampleid = string(s); // sample ID including newline char
                    sampleid = sampleid.substr(0, sampleid.size()-1); // removed newline char
                }
                // naive search of the sample in the complete list of sample IDs to be excluded
                bool include = true;
                for (const auto& skipid : skipids) {
                    if (skipid.compare(sampleid) == 0) {
                        include = false;
                        break;
                    }
                }
                // print included sample ID
                if (include) {
                    cout << "\t" << sampleid;
                } else { // store index of excluded sample
                    skipidxs.push_back(curridx);
                    cerr << "\t" << sampleid << endl;
                }
                curridx++;
                s = se;
            }
            // newline at end
            cout << endl;

            nsamples = curridx;
            cerr << "Read " << nsamples << " samples from header, of these " << skipidxs.size() << " will be skipped." << endl;
        }

        // parse rest of file
        cerr << "Processing..." << endl;

        vector<pair<char*,char*>> outptrs;
        outptrs.reserve(skipidxs.size()+3); // all sample blocks around the skipped samples + fields before INFO + FORMAT field after INFO
        ssize_t nline = 0;
        while((nline = getline(&line, &len, stdin)) != -1) {

            size_t ac = 0;
            size_t an = 0;
            size_t missc = 0;

            nvars++;
            if (nvars % 10000 == 0)
                cerr << "  " << nvars << " variants..." << endl;

            outptrs.clear();
            char* lineend = line + nline; // points to the end of the line (after the newline character)

            // print everything until beginning of INFO column (8th column)
            char* s = line;
            for (int t = 0; t < 7; t++) {
                s = strchr(s+1, '\t');
            }
            s++; // points to beginning of INFO column now
            outptrs.emplace_back(line, s); // will print the beginning until INFO including tab character at end

//            // DEBUG
//            cerr << string(line, s-line);

            s = strchr(s, '\t'); // move forward to end of INFO (skips INFO, points to tab before FORMAT now)

            char* se = strchr(s+1, '\t'); // is either NULL (no samples) or points to tab before first sample
            if (!se)
                se = lineend;

            // include the FORMAT field as first block for the output
            outptrs.emplace_back(s, se);
            s = se;

            // parse the rest of the line, count alleles and skip samples from list
            size_t curridx = 0;
            auto skipidxit = skipidxs.cbegin();
            size_t nextskipidx;
            // set next skip idx to number of samples if nothing to skip
            if (skipidxit != skipidxs.cend()) {
                nextskipidx = *skipidxit;
                skipidxit++; // always points to the next one after the current
            } else
                nextskipidx = nsamples;
            while (s < lineend) {

                char* sc = s; // tab before sample column
                while (curridx < nextskipidx) { // until we reach the next sample to skip (or until we handled the last sample)

                    sc++; // first character of sample column

                    // find end of current sample column first
                    char* sce = strchr(sc, '\t');
                    if (!sce)
                        sce = lineend;

                    // for further parsing, we need to temporarily replace the end by the null character (if we are at the end, it points to null anyway)
                    *sce = '\0';

                    // we assume the first character represents the first allele, all alleles different from 0 or missing (.) are counted
                    if (*sc != '.') {
                        if (*sc != '0')
                            ac++;
                        an++;
                    } else
                        missc++;

                    // diploid and phased?
                    char* sc2 = strchr(sc, '|');
                    if (sc2) { // found -> count
                        sc2++; // points to allele
                        if (*sc2 != '.') {
                            if (*sc2 != '0')
                                ac++;
                        } else
                            missc++;
                        an++;
                    } else { // not found -> search for unphased
                        sc2 = strchr(sc, '/');
                        if (sc2) { // found -> count
                            sc2++; // points to allele
                            if (*sc2 != '.') {
                                if (*sc2 != '0')
                                    ac++;
                            } else
                                missc++;
                            an++;
                        }
                    }

                    // restore the tab character
                    if (sce != lineend)
                        *sce = '\t';

                    sc = sce;
                    curridx++;
                } // END while (curridx < nextskipidx)

                // insert block for output
                outptrs.emplace_back(s, sc);
                s = sc;

                // skip the next sample (if we are not at the end yet)
                if (s < lineend) {
                    // find end of sample column
                    s = strchr(s+1, '\t');
                    if (!s)
                        s = lineend;
                    curridx++;
                }

                // set next index to skip
                if (skipidxit != skipidxs.cend()) {
                    nextskipidx = *skipidxit;
                    skipidxit++;
                } else
                    nextskipidx = nsamples;

            } // END while (s < lineend)

            // apply filters
            bool pass = true;
            if (macfilter && ac < macfilter)
                pass = false;
            if (maffilter && ((float)ac)/an < maffilter)
                pass = false;
            if (missfilter && ((float)missc)/an >= missfilter)
                pass = false;

            // print if filters passed
            if (pass) {

                auto outit = outptrs.cbegin();

                // print first fields (before INFO, block already contains the tab character at the end)
                cout << string(outit->first, outit->second - outit->first);
                outit++;

                // print INFO field now
                cout << "AC=" << ac << ";AN=" << an;

    //            // DEBUG
    //            cerr << "AC=" << ac << ";AN=" << an << endl;

                // print all remaining blocks (all blocks should begin with a tab, the last block should contain the newline character)
                for (; outit != outptrs.cend(); ++outit)
                    cout << string(outit->first, outit->second - outit->first);

            } else // if (!pass)
                nskip++;

        } // END while(getline)

        cout << flush;

    } // END contains data

    cerr << "Number of processed variants:             " << nvars << endl;
    cerr << " Of these skipped due to applied filters: " << nskip << endl;
    cerr << " Total variants in output:                " << nvars - nskip << endl;

    free(line);

}

