#include <cstdlib>
#include <ctime>
#include <iostream>
#include <mpi.h>

/**
 * mpi-hello.cpp: simple MPI program to illustrate setting up and shutting
 *     down the MPI communicator. Each processor in the communicator gets the
 *     local time and prints it out to the standard output.
 */
int main(int argc, char** argv) {
	using namespace std;

	int rank;		// processor's rank ID
	int nProcs;		// number of processors in the communicator

	// initialize MPI constructs
	MPI_Init(&argc, &argv);					// set up the communicator
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	// which processor am I?
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);	// how many processors are there?

	// get local time and print it
	time_t epochTime;
	time(&epochTime);
	struct tm *timeInfo;
	timeInfo = localtime(&epochTime);
	cout << "Processor " << rank << " / " << nProcs << " time: " <<
		asctime(timeInfo);

	// shut down the communicator
	MPI_Finalize();

	return EXIT_SUCCESS;
}