test:	test.o
	mpic++ -O0 -g -o test test.o -lpng

test.o:	test.cpp diamond.h PngWriter.hpp
	mpic++ -O0 -g -std=c++11 -I/usr/lib/openmpi/include -c $< 

clean:
	rm -f *.o test
