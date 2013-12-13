ODEPS=bwa/bwtindex.o bwa/is.o bwa/bwt_gen.o bwa/QSufSort.o
CFLAGS=
LFLAGS=-L bwa/ -lbwa -lz -pthread
OBJS=IndexBWA.cpp IndexProtein.cpp
FILE=main_protein

compile:
	g++ $(FILE).cpp $(CFLAGS) $(OBJS) $(LFLAGS) $(ODEPS)

clean:
	make -C bwa clean
	rm -f data/*.fa.*
	rm -f data/*.fasta.*
	rm -f *.out

deps:
	make -C bwa

tags:
	find . -type f -iname "*.[ch]*" | xargs etags -a
