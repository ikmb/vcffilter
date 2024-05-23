/*
 * restorevcf.c
 *
 *  Created on: May 15, 2024
 *      Author: lwienbrandt
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define BUFSIZE 1073741824

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

    size_t len = BUFSIZE;
    const size_t lenstart = len;
    char *line = malloc(len*sizeof(char));
    size_t nline = 0;

    // parse header
    size_t nh = getline(&line, &len, stdin); // read header line
    if (nh == -1) // file does not contain data -> finished
        return 0;
    // overwrite newline character at the end of the line (prevents correct parsing below)
    *(line+nh-1) = '\0';

    // search for more args (separated by tab)
    char* chrend = strchr(line, '\t');
    // null terminate chromosome name
    if (chrend != NULL) {
        *chrend = '\0';
    }
    size_t chrlen = strlen(line); // length of chromosome name
    char* chrom = malloc((chrlen+1)*sizeof(char));
    strcpy(chrom, line); // copy chrom

    // parse more args
    int parsegq = 0;
//    int parseaa = 0;
    char* arg = (chrend != NULL) ? chrend+1 : NULL;
    while (arg != NULL) {
        char* argend = strchr(arg, ';'); // separator for the args in the header is the semicolon char
        if (argend != NULL) {
            *argend = '\0';
        }
        if (strcmp(arg, "--gq") == 0) // --gq option was set -> restore GQ field
            parsegq = 1;
//        if (strcmp(arg, "--aa") == 0) // --aa option was set -> restore AAScore field
//            parseaa = 1;
        arg = (argend != NULL) ? argend+1 : NULL;
    }

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
        int nalt = 0;
        *altallend = '\0'; // null terminate the allele string (to stop following loop)
        for (char* t = altall; t != NULL; t = strchr(t+1, ',')) // will stop at altallend as we have null terminated the buffer there
            nalt++;
        *altallend = '\t'; // restore tab

        // initialize allele counters
        size_t* ac = (size_t*) calloc(nalt, sizeof(size_t));
        size_t an = 0;

        // qual field
        char* qual = altallend+1; // start
        char* qualend = strchr(qual, '\t'); // end (exclusive)
//        *qualend = '\0'; // null terminate

        // filter field
        char* filter = qualend+1; // start of filter field
        char* filterend = strchr(filter, '\t'); // end of filter (exclusive)
//        *filterend = '\0'; // null terminate

        // info
        char* info = filterend+1; // start
        char* infoend = strchr(info, '\t'); // end
//        *infoend = '\0'; // null terminate

//        // AAScore (if desired)
//        char* aa = filterend+1; // points to AAScore (if enabled) or to first genotype
//        char* aaend = filterend;
//        if (parseaa) {
//            aaend = strchr(aa, '\t'); // end of AAScore
//            *aaend = '\0'; // null terminate AAScore
//        }

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

        // print VCF line

        // CHROM (chromosome name)
        fputs(chrom, stdout);
        printf("\t");

        // POS, ID, alleles, QUAL, FILTER
        *filterend = '\0'; // null terminate filter field, as we are going to modify the following INFO fields
        fputs(pos, stdout);

//        printf("\t.\t"); // ID (set to unknown)
//
//        fputs(posend+1, stdout); // REF ALT (alleles)
//
////        printf("\t.\t"); // QUAL (set to unknown)
//        printf("\t");
//        fputs(qual, stdout);
//        printf("\t");
//
//        fputs(filter, stdout); // FILTER

        // INFO
        float anf = (float) an;
        printf("\tAF=%.8f", ac[0]/anf); // AF of first alt allele
        for (int n = 1; n < nalt; n++)
            printf(",%.8f", ac[n]/anf); // AF of further alleles if multi-allelic
        printf(";AC=%lu", ac[0]); // AC of first alt allele
        for (int n = 1; n < nalt; n++)
            printf(",%lu", ac[n]); // AC of further alleles if multi-allelic
        printf(";AN=%lu", an); // AN
//        if (parseaa) { // restore AA Score
//            printf(";");
//            fputs(aa, stdout);
//        }
        printf(";");
        *infoend = '\0'; // null terminate INFO field

        // find original values of the replaced ones above -> prefix them with "Org"
        char* org[3];
        org[0] = findInfoField(info, "AF=");
        org[1] = findInfoField(info, "AC=");
        org[2] = findInfoField(info, "AN=");
        // sort...
        qsort(org, 3, sizeof(char*), ptrcmp);

        // print INFO and prefix above fields with 'Org'
        const char repl[] = "OrgA"; // we add the 'A' here, as we replace it below with the null terminator...
        char* infoit = info;
        for (int i = 0; i < 3; i++) {
            if (org[i] != NULL) {
                *(org[i]) = '\0'; // temporarily add null terminator (luckily, we know that we replace an 'A' here...)
                fputs(infoit, stdout);
                fputs(repl, stdout);
                infoit = org[i]+1; // next char after the null terminator
            }
        }
        // print remainder
        fputs(infoit, stdout);

        // FORMAT
        if (parsegq)
            printf("\tGT:GQ\t");
        else
            printf("\tGT\t");

        // genotypes
        fputs(gtstart, stdout); // rest of the line buffer containing the genotypes, ends with newline!
        fflush(stdout);

        free(ac);
        nline++;
    }

    fprintf(stderr, "Number of lines: %lu\n", nline);
    fprintf(stderr, "Line buffer size: %lu", len);
    if (len != lenstart)
        fprintf(stderr, " -> changed!!\n");
    else
        fprintf(stderr, "\n");

    free(chrom);
    free(line);

}

