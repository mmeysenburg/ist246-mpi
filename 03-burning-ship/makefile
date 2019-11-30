bin/frac.mpi:	bin/frac.o bin/EasyBMP.o
	mpic++ -std=c++11 -Wall -O3 bin/frac.o bin/EasyBMP.o -o bin/frac.mpi
	
bin/frac.o:	frac.cpp
	mpic++ -std=c++11 -Wall -O3 -c frac.cpp -o bin/frac.o
	
bin/EasyBMP.o:	EasyBMP.cpp
	mpic++ -std=c++11 -Wall -O3 -c EasyBMP.cpp -o bin/EasyBMP.o
	
clean:
	rm bin/*.o bin/frac.mpi
