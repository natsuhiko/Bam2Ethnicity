
CC=gcc
CFLAGS := $(CFLAGS) -std=gnu99 -I/usr/include -fpic -g -O2
LDFLAGS := $(LDFLAGS)  -lm -lz




.SUFFIXES: .c .o

.PHONY: all
all: qcFilterBam.o pasteFiles.o zpaste.o countAS.o
	$(CC) -o qcFilterBam qcFilterBam.o $(LDFLAGS)
	$(CC) -o pasteFiles pasteFiles.o $(LDFLAGS)
	$(CC) -o zpaste zpaste.o $(LDFLAGS)
	$(CC) -o countAS countAS.o $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm qcFilterBam pasteFiles zpaste countAS *.o
