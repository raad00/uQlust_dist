CC = gcc

# gcc options
	CFLAGS = -O3 -lm 

	FLAGS = -O3 
LDFLAGS= -shared
SOURCES = $(shell echo *.c)
HEADERS = $(shell echo *.h)
OBJECTS = $(SOURCES:.c=.o)

TARGET = libRNA.so

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)


$(TARGET) : analyze.o fpair.o fpair_sub.o nrutil.o pair_type.o rnaview.o
	$(CC) -shared $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)
	cp $(TARGET) ..

analyze.o: rna_header.h analyze.c
	$(CC)  -c -fPIC $(CFLAGS) analyze.c
fpair.o:  rna_header.h fpair.c
	$(CC) -c -fPIC $(CFLAGS) fpair.c
fpair_sub.o: rna_header.h fpair_sub.c
	$(CC) -c -fPIC $(CFLAGS)  fpair_sub.c
nrutil.o: nrutil.h nrutil.c
	$(CC) -c -fPIC $(CFLAGS) nrutil.c
pair_type.o: rna_header.h pair_type.c
	$(CC) -c -fPIC $(CFLAGS) pair_type.c
rnaview.o: rna_header.h rnaview.c
	$(CC) -c -fPIC $(CFLAGS) rnaview.c

	
