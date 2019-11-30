#include<cstdlib>
#include <ctime> 
#include<iostream>
#include <mpi.h>

/**
 * Simple deadlock avoidance example. Even numbered processors send to odds,
 * then receive from odds, while odd numbered processors receive from evens,
 * then send to evens.
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

	if (rank % 2 == 0) {
		// even processes send to odds, then receive from odds

		// send this processor's array to odd processes
		for (int i = 1; i < nProcs; i += 2) {
			if (i != rank) {
				cout << "Sending " << n << " numbers to processor " << i << endl;
				MPI_Send(pSend, n, MPI::INT, i, MSG_NUMBERS, MPI_COMM_WORLD);
			}
		}

		// receive data from odd numbered processors
		MPI_Status status;
		for (int i = 1; i < nProcs; i += 2) {
			if (i != rank) {
				cout << "Receiving " << n << " numbers from processor " << i << endl;
				MPI_Recv(pRecv, n, MPI::INT, i, MSG_NUMBERS, MPI_COMM_WORLD, &status);
			}
		}
	} else {
		// odd processors reveive from evens, then send to evens

		// receive data from even numbered processors
		MPI_Status status;
		for (int i = 0; i < nProcs; i += 2) {
			if (i != rank) {
				cout << "Receiving " << n << " numbers from processor " << i << endl;
				MPI_Recv(pRecv, n, MPI::INT, i, MSG_NUMBERS, MPI_COMM_WORLD, &status);
			}
		}

		// send this processor's array to even processes
		for (int i = 0; i < nProcs; i += 2) {
			if (i != rank) {
				cout << "Sending " << n << " numbers to processor " << i << endl;
				MPI_Send(pSend, n, MPI::INT, i, MSG_NUMBERS, MPI_COMM_WORLD);
			}
		}
	}

	// shut down MPI 
	MPI_Finalize();

	// free memory
	delete[] pSend;
	delete[] pRecv;

    return EXIT_SUCCESS;
}