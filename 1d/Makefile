default:	kuramoto_classic kuramoto_diamond

kuramoto_classic:	kuramoto_classic.o
	mpic++ -O3 -o kuramoto_classic kuramoto_classic.o -L/master/home/qiqi/lib -lpng

kuramoto_diamond:    kuramoto_diamond.o
	mpic++ -O3 -o kuramoto_diamond kuramoto_diamond.o -L/master/home/qiqi/lib -lpng

kuramoto_classic.o:  kuramoto_classic.cpp kuramoto.h classic_scheme.h diamond.h PngWriter.hpp
	mpic++ -O3 -std=c++11 -I/master/home/qiqi/include/libpng16 -c $< 

kuramoto_diamond.o:  kuramoto_diamond.cpp kuramoto.h diamond_scheme.h diamond.h PngWriter.hpp
	mpic++ -O3 -std=c++11 -I/master/home/qiqi/include/libpng16 -c $< 

run:	kuramoto_classic kuramoto_diamond
	mpirun ./kuramoto_classic 8
	mpirun ./kuramoto_classic 16
	mpirun ./kuramoto_classic 32
	mpirun ./kuramoto_classic 64
	mpirun ./kuramoto_classic 128
	mpirun ./kuramoto_classic 256
	mpirun ./kuramoto_classic 512
	mpirun ./kuramoto_classic 1024
	mpirun ./kuramoto_classic 2048
	mpirun ./kuramoto_classic 4096
	mpirun ./kuramoto_classic 8192
	mpirun ./kuramoto_classic 16384
	mpirun ./kuramoto_classic 32768
	mpirun ./kuramoto_diamond 8
	mpirun ./kuramoto_diamond 16
	mpirun ./kuramoto_diamond 32
	mpirun ./kuramoto_diamond 64
	mpirun ./kuramoto_diamond 128
	mpirun ./kuramoto_diamond 256
	mpirun ./kuramoto_diamond 512
	mpirun ./kuramoto_diamond 1024
	mpirun ./kuramoto_diamond 2048
	mpirun ./kuramoto_diamond 4096
	mpirun ./kuramoto_diamond 8192
	mpirun ./kuramoto_diamond 16384
	mpirun ./kuramoto_diamond 32768

clean:
	rm -f *.o kuramoto_diamond kuramoto_classic
