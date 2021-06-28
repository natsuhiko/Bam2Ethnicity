
CC=gcc
CFLAGS := ${CFLAGS} -std=gnu99 -fpic -g -O2
LDFLAGS := ${LDFLAGS} -lhts -lgsl -lgslcblas -llapack -lm -lz -Wl,-rpath,/software/team170/miniconda3/lib

PROGRAM = vcf2eigen
OBJS = parseVCF.o vcf2eigen.o

.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $^ $(LDFLAGS) $(PROGFLAG)

.c.o:
	$(CC) $(CFLAGS) -c $< $(PROGFLAG)



