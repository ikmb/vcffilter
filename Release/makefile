RM := rm -rf

all: vcffilter myzcat restorevcf

vcffilter:
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"vcffilter.d" -MT"vcffilter.o" -o "vcffilter.o" "../vcffilter.c"
	gcc  -o "vcffilter" "./vcffilter.o"

myzcat:
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"myzcat.d" -MT"myzcat.o" -o "myzcat.o" "../myzcat.c"
	gcc  -o "myzcat" "./myzcat.o" -lz
	
restorevcf:
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"restorevcf.d" -MT"restorevcf.o" -o "restorevcf.o" "../restorevcf.c"
	gcc  -o "restorevcf" "./restorevcf.o"
	
clean:
	-$(RM) vcffilter* myzcat* restorevcf*

