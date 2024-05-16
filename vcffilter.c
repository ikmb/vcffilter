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

            // genotypes (skip INFO and FORMAT)
            char* gtstart = strchr(filter+3, '\t'); // start of info
            gtstart = strchr(gtstart+1, '\t'); // start of format
            gtstart = strchr(gtstart+1, '\t'); // start of genotypes (pointing at tab at beginning!)
            char* gtend = gtstart; // just for the beginning
            for (; gtstart != NULL; gtstart = strchr(gtend+1, '\t')) { // always take next genotypes
                gtend = strchr(gtstart+1,':'); // end of genotype (exclusive)
                *gtend = '\0'; // null terminate genotype
                fputs(gtstart, stdout); // print genotype
            }

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

