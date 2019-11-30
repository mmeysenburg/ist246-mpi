all: bin/mpi-deadlock1.mpi bin/mpi-deadlock2.mpi

bin/mpi-deadlock1.mpi:	deadlock1.cpp
	mpic++ -std=c++11 -O3 deadlock1.cpp -o bin/mpi-deadlock1.mpi

bin/mpi-deadlock2.mpi:	deadlock2.cpp
	mpic++ -std=c++11 -O3 deadlock2.cpp -o bin/mpi-deadlock2.mpi

clean:
	rm bin/mpi-deadlock1.mpi bin/mpi-deadlock2.mpi
