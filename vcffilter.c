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

    // parse args
    int parsegq = 0;
//    int parseaa = 0;

    char** cargv = argv+1; // to first arg
    int cargc = argc-1;
    while (cargc) {
        if (strcmp(*cargv, "--gq") == 0)
            parsegq = 1;
//        else if (strcmp(*cargv, "--aa") == 0)
//            parseaa = 1;
        cargc--;
        cargv++;
    }

    size_t len = BUFSIZE;
    const size_t lenstart = len;
    char *line = malloc(len*sizeof(char));
    size_t nline = 0;

    // skip header and print chromosome to output
    while(getline(&line, &len, stdin) != -1) {
        if (*line != '#') { // just found the first line after the header
            // copy chromosome name
            char* chromend = strchr(line, '\t');
            *chromend = '\0'; // null terminate chromosome name
            fputs(line, stdout); // print chromosome name
            *chromend = '\t'; // restore tab character for further processing below
            break;
        }
    }

    // print <args> after chromosome, using ';' as separator
    for (int i=1; i < argc; i++) {
        printf(";%s", argv[i]);
    }
    printf("\n");

    // parse rest of file
    do {
        // genomic position
        char* posstart = strchr(line, '\t')+1; // start of genomic position (skipped chromosome name)
        char* posend = strchr(posstart, '\t'); // end of genomic position (exclusive, points to tab char)
//        *posend = '\0'; // null terminate the pos string
//        fputs(posstart, stdout); // print pos string

        // variant ID
        char* varidstart = posend+1; // start
        char* varidend = strchr(varidstart, '\t'); // end (tab)

        // alleles
        char* allstart = varidend+1; // start of wild type allele
        char* allend = strchr(allstart, '\t'); // end of first allele (exclusive)
        allend = strchr(allend+1, '\t'); // end of second allele (exclusive)
//        *allend = '\0'; // null terminate the allele string
//        fputs(allstart, stdout); // print the alleles

        // QUAL column
        char* qual = allend+1; // start of qual
        char* qualend = strchr(qual, '\t'); // end of qual (tab)
//        *qualend = '\0'; // null terminate qual field
//        printf("\t");
//        fputs(qual, stdout); // print qual

        // FILTER column
        char* filter = qualend+1; // start of filter field, pointing to first char
        char* filterend = strchr(filter, '\t'); // end of filter field (tab)
//        *filterend = '\0'; // null terminate the filter field
//        printf("\t");
//        fputs(filter, stdout); // print filter

        // INFO column // TODO maybe add switch to be able to exclude INFO?
        char* info = filterend+1; // beginning of INFO column
        char* infoend = strchr(info, '\t'); // end of INFO (exclusive)
        *infoend = '\0'; // null terminate info field
        fputs(posstart, stdout); // print all fields from position to INFO (inclusive)

//        // parse for AAScore, if desired
//        if (parseaa) {
//            // search for AAScore:
//            // need to print a separate tab character anyway, even if we do not find the AAScore (otherwise, restore won't work)
//            printf("\t");
//            char* aa = info;
//            char* aaend = info; // init
//            *infoend = ';'; // terminate to stop the following while loop
//            while (aaend != infoend) { // until the end of the info field
//                aaend = strchr(aa, ';'); // find end of first field
//                if (aa[0] == 'A' && aa[1] == 'A') { // assuming no other info field starts with AA
//                    *aaend = '\0'; // null terminate AAScore
//                    // print AAScore
//                    fputs(aa, stdout);
//                    break; // found, no need to search further
//                }
//                aa = aaend+1;
//            }
//            //*infoend = '\t'; // restore tab -> not necessary
//        }

        // parse FORMAT field for GQ if desired
        char* fmt = infoend+1; // start of format, pointing at first char in format field!
        char* fmtend = strchr(fmt, '\t'); // end of format field (pointing at tab)
        int gqidx = -1; // disabled GQ parsing until we find the GQ field
        if (parsegq) {
            // search for position of GQ
            int gqidxtmp = 0;
            char* fmt2 = fmt; // init
            *fmtend = ':'; // need to terminate the format field with ':' to be able to always end the following searches here
            while (fmt2 != fmtend) { // until we reached the end of the format field
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
            *fmtend = '\t'; // restore tab character at the beginning of genotypes: in the case no GQ was found, the process below awaits a tab at the beginning!
        }

        // genotypes -> assuming GT is the first field!
        char* gtstart = fmtend; // tab at beginning

        if (gqidx < 0) { // if we do not want to add the GQ field
            // simple and fast version as long as we are not requesting the GQ field
            // (works only if there is additional info attached to each genotype separated by ':')

            char* gtend = gtstart; // just for the beginning
            for (; gtstart != NULL; gtstart = strchr(gtend+1, '\t')) { // find next genotypes
                gtend = strchr(gtstart+1,':'); // end of genotype (exclusive)
                *gtend = '\0'; // null terminate genotype
                fputs(gtstart, stdout); // print genotype (including beginning '\t')
            }

        } else { // add GQ field
            // slower more deperate version

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

                gtstart = gtallend; // note, that gtstart points to a ':' or '\0' character!
            }

        }

        printf("\n"); // newline at the end
        nline++;

    } while(getline(&line, &len, stdin) != -1);

    fprintf(stderr, "Number of variants: %lu\n", nline);
    fprintf(stderr, "Line buffer size: %lu", len);
    if (len != lenstart)
        fprintf(stderr, " -> changed!!\n");
    else
        fprintf(stderr, "\n");

    free(line);

}

