all: bin/frac.mpi bin/uni-frac

bin/frac.mpi:	bin/frac.o bin/EasyBMP.o
	mpic++ -std=c++11 -Wall -O3 bin/frac.o bin/EasyBMP.o -o bin/frac.mpi
	
bin/frac.o:	frac.cpp
	mpic++ -std=c++11 -Wall -O3 -c frac.cpp -o bin/frac.o
	
bin/EasyBMP.o:	EasyBMP.cpp
	mpic++ -std=c++11 -Wall -O3 -c EasyBMP.cpp -o bin/EasyBMP.o
	
bin/uni-frac:	bin/uni-frac.o bin/EasyBMP.o
	g++ -std=c++11 -Wall -O3 bin/uni-frac.o bin/EasyBMP.o -o bin/uni-frac
	
bin/uni-frac.o:	uni-frac.cpp
	g++ -std=c++11 -Wall -O3 -c uni-frac.cpp -o bin/uni-frac.o
	
clean:
	rm bin/*.o bin/frac.mpi bin/uni-frac
