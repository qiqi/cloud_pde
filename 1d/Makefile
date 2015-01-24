test:	test.o
	mpic++ -O0 -g -o test test.o -L/master/home/qiqi/lib -lpng

test.o:	test.cpp classic_scheme.h diamond.h PngWriter.hpp
	mpic++ -O0 -g -std=c++11 -I/master/home/qiqi/include/libpng16 -c $< 

run:	test
	mpirun ./test

clean:
	rm -f *.o test
