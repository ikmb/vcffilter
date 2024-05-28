RM := rm -rf

all: vcffilter myzcat

vcffilter:
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"vcffilter.d" -MT"vcffilter.o" -o "vcffilter.o" "../vcffilter.c"
	gcc  -o "vcffilter" "./vcffilter.o"

myzcat:
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"myzcat.d" -MT"myzcat.o" -o "myzcat.o" "../myzcat.c"
	gcc  -o "myzcat" "./myzcat.o" -lz

clean:
	-$(RM) vcffilter* myzcat*

