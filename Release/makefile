RM := rm -rf

all: vcffilter myzcat restorevcf removesamples

vcffilter:
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"vcffilter.d" -MT"vcffilter.o" -o "vcffilter.o" "../vcffilter.c"
	gcc  -o "vcffilter" "./vcffilter.o"

myzcat:
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"myzcat.d" -MT"myzcat.o" -o "myzcat.o" "../myzcat.c"
	gcc  -o "myzcat" "./myzcat.o" -lz

restorevcf:
	$(MAKE) -C ../restorevcf/Release all

removesamples:
	$(MAKE) -C ../removesamples/Release all

clean:
	$(RM) vcffilter* myzcat*
	$(MAKE) -C ../restorevcf/Release clean
	$(MAKE) -C ../removesamples/Release clean
