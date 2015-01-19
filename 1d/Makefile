test:	test.o
	g++ -O0 -g -o test test.o -lpng

test.o:	test.cpp diamond.h PngWriter.hpp
	g++ -O0 -g -std=c++11 -c $< 

clean:
	rm -f *.o test
