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
    int parseaa = 0;
    char* arg = (chrend != NULL) ? chrend+1 : NULL;
    while (arg != NULL) {
        char* argend = strchr(arg, '\t');
        if (argend != NULL) {
            *argend = '\0';
        }
        if (strcmp(arg, "--gq") == 0) // --gq option was set -> restore GQ field
            parsegq = 1;
        if (strcmp(arg, "--aa") == 0) // --aa option was set -> restore AAScore field
            parseaa = 1;
        arg = (argend != NULL) ? argend+1 : NULL;
    }

    // parse rest of file
    while(getline(&line, &len, stdin) != -1) {

        // genomic position
        char* posend = strchr(line, '\t'); // end of genomic position (exclusive, points to tab char)
        if (posend == NULL) // no tab char -> invalid line (can happen for the last line when it contains only a newline character) -> skip
            continue;
        *posend = '\0'; // null terminate the pos string

        // alleles
        char* refallend = strchr(posend+1, '\t'); // end of first allele (exclusive)
        char* altallend = strchr(refallend+1, '\t'); // end of alternative alleles (exclusive)
        *altallend = '\0'; // null terminate the allele string
        // count number of alternative alles:
        int nalt = 0;
        for (char* t = refallend+1; t != NULL; t = strchr(t+1, ',')) // will stop at altallend as we have null terminated the buffer there
            nalt++;
        // initialize allele counters
        size_t* ac = (size_t*) calloc(nalt, sizeof(size_t));
        size_t an = 0;

        // filter field
        char* filter = altallend+1; // start of filter field
        char* filterend = strchr(filter, '\t'); // end of filter (exclusive)
        *filterend = '\0'; // null terminate filter string

        // AAScore (if desired)
        char* aa = filterend+1; // points to AAScore (if enabled) or to first genotype
        char* aaend = filterend;
        if (parseaa) {
            aaend = strchr(aa, '\t'); // end of AAScore
            *aaend = '\0'; // null terminate AAScore
        }

        // genotypes
        char* gtstart = aaend+1; // start of genotypes (pointing at first gt char!)
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

        fputs(chrom, stdout); // CHROM (chromosome name)
        printf("\t");

        fputs(line, stdout); // POS (genomic position)

        printf("\t.\t"); // ID (set to unknown)

        fputs(posend+1, stdout); // REF ALT (alleles)

        printf("\t.\t"); // QUAL (set to unknown)

        fputs(filter, stdout); // FILTER

        // INFO
        float anf = (float) an;
        printf("\tAF=%.8f", ac[0]/anf); // AF of first alt allele
        for (int n = 1; n < nalt; n++)
            printf(",%.8f", ac[n]/anf); // AF of further alleles if multi-allelic
        printf(";AC=%lu", ac[0]); // AC of first alt allele
        for (int n = 1; n < nalt; n++)
            printf(",%lu", ac[n]); // AC of further alleles if multi-allelic
        printf(";AN=%lu", an); // AN
        if (parseaa) { // restore AA Score
            printf(";");
            fputs(aa, stdout);
        }

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

