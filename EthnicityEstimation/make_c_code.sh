gcc -c parseVCF.c

gcc -o vcf2eigen vcf2eigen.c parseVCF.o -lz -lm
