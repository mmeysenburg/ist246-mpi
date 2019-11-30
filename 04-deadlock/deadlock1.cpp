#include<cstdlib>
#include <ctime> 
#include<iostream>
#include <mpi.h>

/**
 * Deadlock example -- all processors try to send to another process at the 
 * same time, resulting in a deadlock situation. Have your CTRL+C or scancel
 * command ready!
 * 
 * \param argC number of command-line arguments
 *
 * \param ppArgv command-line arguments, as an array of C strings
 *
 * \return EXIT_SUCCESS if everthing went well
 */
int main(int argC, char** ppArgv) {
    using namespace std;

	// Tag for MPI send and receive operations
	int MSG_NUMBERS = 0;

	int rank;    // this processor's rank in the communicator
	int nProcs;  // number of processors in the communicator

	// initialize MPI constructs
	MPI_Init(&argC, &ppArgv);                // set up the communicator
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);    // which processor am I?
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);  // how many processors are there?

	// create and fill an array to send
	srand(time(0));
	int n = 1000000;
	int *pSend = new int[n];
	for (int i = 0; i < n; i++) {
		pSend[i] = rand();
	}

	// create an array to receive
	int *pRecv = new int[n];

	// send this processor's array to all other processes
	for (int i = 0; i < nProcs; i++) {
		if (i != rank) {
			cout << "Sending " << n << " numbers to processor " << i << endl;
			MPI_Send(pSend, n, MPI::INT, i, MSG_NUMBERS, MPI_COMM_WORLD);
		}
	}

	// receive data from other processors
	MPI_Status status;
	for (int i = 0; i < nProcs; i++) {
		if (i != rank) {
			cout << "Receiving " << n << " numbers from processor " << i << endl;
			MPI_Recv(pRecv, n, MPI::INT, i, MSG_NUMBERS, MPI_COMM_WORLD, &status);
		}
	}

	// shut down MPI 
	MPI_Finalize();

	// free memory
	delete[] pSend;
	delete[] pRecv;

    return EXIT_SUCCESS;
}