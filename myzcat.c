/*
 *    Copyright (C) 2024 by Lars Wienbrandt,
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

// using a large buffer
#define BUFSIZE 1073741824

/* compress or decompress from fin (command line argument) to stdout */
int main(int argc, char **argv)
{
    if(argc <= 1) { // <= (number of expected CLI arguments)
        fprintf(stderr, "Usage: %s <input files>\n", argv[0]);
        return -1;
    }

    char* buf = malloc(BUFSIZE * sizeof(char));

    for (int i = 1; i < argc; i++) {
        gzFile fin = gzopen(argv[i], "rb");
        size_t n;
        while((n = gzread(fin, buf, BUFSIZE)) > 0) {
            fwrite(buf, sizeof(char), n, stdout);
        }
        gzclose(fin);
    }

    free(buf);

    return 0;
}

