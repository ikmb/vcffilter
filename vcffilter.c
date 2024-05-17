/*
 * vcffilter.c
 *
 *  Created on: May 13, 2024
 *      Author: lwienbrandt
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define BUFSIZE 1073741824

int main (int argc, char **argv) {

    int parsegq = 0;
    if (argc > 1 && strcmp(argv[1], "-gq") == 0)
        parsegq = 1;

//    FILE *fp = NULL;
//    if (strcmp(argv[1], "-") != 0) {
//        fp = fopen(argv[1], "r");
//        if(fp == NULL) {
//            perror("Unable to open file!");
//            exit(1);
//        }
//    }

    size_t len = BUFSIZE;
    const size_t lenstart = len;
    char *line = malloc(len*sizeof(char));
    size_t nline = 0;

    //while(getline(&line, &len, fp) != -1) {
    while(getline(&line, &len, stdin) != -1) {
        if (*line != '#') { // skip header

            // genomic position
            char* posstart = strchr(line, '\t')+1; // start of genomic position (skipped chromosome name)
            char* posend = strchr(posstart, '\t'); // end of genomic position (exclusive, points to tab char)
            *posend = '\0'; // null terminate the pos string
            fputs(posstart, stdout); // print pos string

            // alleles
            char* allstart = strchr(posend+1, '\t'); // start of wild type allele (skipped variant ID), pointing at tab character!
            char* allend = strchr(allstart+1, '\t'); // end of first allele (exclusive)
            allend = strchr(allend+1, '\t'); // end of second allele (exclusive)
            *allend = '\0'; // null terminate the allele string
            fputs(allstart, stdout); // print the alleles

            // print first character of filter field (must be either P for PASS or L for LowXXX)
            char* filter = strchr(allend+1, '\t'); // start of filter field, pointing to tab! (skipped QUAL field)
            *(filter+2) = '\0'; // null terminate the filter field (overwrite some character following the first)
            fputs(filter, stdout); // print filter

            // parse FORMAT field for GQ if desired (skip INFO)
            char* fmt = strchr(filter+3, '\t'); // start of info
            fmt = strchr(fmt+1, '\t') + 1; // start of format (pointing at first char in format field)
            char* gtstart = strchr(fmt, '\t'); // start of genotypes (pointing at tab at beginning!)
            int gqidx = -1; // disabled GQ parsing until we find the GQ field
            if (parsegq) {
                // search for position of GQ
                int gqidxtmp = 0;
                char* fmt2 = fmt; // init
                *gtstart = ':'; // need to terminate the format field with ':' to be able to always end the following searches here
                while (fmt2 != gtstart) { // until we reached the end of the format field
                    fmt2 = strchr(fmt, ':'); // find end of format description field -> will not be NULL here as we terminated format above with ':'
                    // null terminate the actual format description field
                    *fmt2 = '\0';
                    if (strcmp(fmt, "GQ") == 0) { // found!
                        // leave while, gqidx is the index of the GQ field
                        gqidx = gqidxtmp;
                        break;
                    }
                    // else -> next format description field
                    gqidxtmp++;
                    fmt = fmt2+1;
                }
                //*gtstart = '\t'; // restore tab character at the beginning of genotypes -> not necessary!
            }

            // genotypes (start was already searched for above) -> assuming GT is the first field!
            int finished = 0;
            while (!finished) {
                char* gtallend = strchr(gtstart+1, '\t');
                if (gtallend == NULL) { // '\t' not found -> must be at the end of the line
                    gtallend = strchr(gtstart+1, '\0'); // points to termination char now
                    finished = 1; // finished hereafter
                }
                *gtallend = ':'; // terminate gt field with ':' (just in case there are only GTs and not anything else separated by ':')
                char* gtend = strchr(gtstart+1, ':'); // find end of GTs
                *gtstart = '\t'; // write tab at the beginning
                *gtend = '\0'; // null terminate GTs
                fputs(gtstart, stdout); // print genotype (including beginning '\t')

                if (gqidx >= 0) { // if we want to add the GQ field
                    int gqidxtmp = 0; // now still pointing at GT
                    while (gqidxtmp < gqidx) { // move to GQ field
                        gtstart = gtend+1;
                        gtend = strchr(gtstart, ':'); // note, we terminated the GT field with ':'
                        gqidxtmp++;
                    }
                    // now, gtstart points to beginning of GQ, gtend to char after GQ
                    *(gtstart-1) = ':'; // write ':' at the beginning
                    *gtend = '\0'; // null terminate GQ
                    fputs(gtstart-1, stdout);
                }

                gtstart = gtallend; // note, that gtstart points to a ':' or '\0' character!
            }

//            // old
//
//            char* gtend = gtstart; // just for the beginning
//            for (; gtstart != NULL; gtstart = strchr(gtend+1, '\t')) { // find next genotypes
//                gtend = strchr(gtstart+1,':'); // end of genotype (exclusive)
//                *gtend = '\0'; // null terminate genotype
//                fputs(gtstart, stdout); // print genotype (including beginning '\t')
//                if (gqidx >= 0) { // if we want to add the GQ field
//                    int gqidxtmp = 0; // now still pointing at GT
//                    while (gqidxtmp < gqidx) { // move to GQ field
//                        gtstart = gtend+1;
//                        gtend = strchr(gtstart, ':'); // TODO moves into next field if GQ is last field! no : at the end!
//                        gqidxtmp++;
//                    }
//                }
//            }


            printf("\n"); // newline at the end

        }
        nline++;
    }

    fprintf(stderr, "Number of lines: %lu\n", nline);
    fprintf(stderr, "Line buffer size: %lu", len);
    if (len != lenstart)
        fprintf(stderr, " -> changed!!\n");
    else
        fprintf(stderr, "\n");

//    fclose(fp);
    free(line);

}

