ODEPS=bwa/bwtindex.o bwa/is.o bwa/bwt_gen.o bwa/QSufSort.o
CFLAGS=
LFLAGS=-L bwa/ -lbwa -lz -pthread
OBJS=IndexBWA.cpp

compile:
	g++ main.cpp $(CFLAGS) $(OBJS) $(LFLAGS) $(ODEPS)

clean:
	make -C bwa clean
	rm -f *.fa.*
	rm -f *.out

deps:
	make -C bwa
