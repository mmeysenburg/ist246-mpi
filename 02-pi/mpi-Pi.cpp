#include <chrono>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <random>

int main(int argc, char** argv) {
	using namespace std;

	int rank;		// processor's rank ID
	int nProcs;		// number of processors in the communicator

	int N = 100000000;	// number of darts each processor throws
	int TAG_COUNT = 0;		// tag for sending results back to master processor

	// initialize MPI constructs
	MPI_Init(&argc, &argv);					// set up the communicator
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	// which processor am I?
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);	// how many processors are there?

	if (rank != 0) {
		// non-master processors all throw N darts and count how many land in 
		// the upper-right quadrant of a unit circle centered at (0, 0)

		// each processor needs its own random number generator, using system
		// time plus rank so they're all seeded differently. Here we use the 
		// Mersenne Twister PRNG, producing a uniform distribution of numbers 
		// in (-1, 1)
		mt19937 prng(chrono::system_clock::now().time_since_epoch().count() + rank);
		uniform_real_distribution<long double> dis(-1.0, 1.0);

		// throw the darts
		int count = 0;
		long double x, y;
		for (int i = 0; i < N; i++) {
			x = dis(prng); y = dis(prng);

			// in the circle?
			if ((x * x + y * y) <= 1.0) {
				count++;
			}
		}

		// send this processor's count to the master processor
		MPI_Send(&count, 1, MPI::INT, 0, TAG_COUNT, MPI_COMM_WORLD);
	} else {
		// master process accumulates each slave's count, then completes the estimate
		// of pi
		MPI_Status status;
		int count, accumulator = 0;
		for (int i = 1; i < nProcs; i++) {
			MPI_Recv(&count, 1, MPI::INT, i, TAG_COUNT, MPI_COMM_WORLD, &status);
			accumulator += count;
		}

		// create the estimate
		long double pi = 4.0L * accumulator / ((nProcs - 1) * N);

		// report results
		cout << "Estimate of pi: " << pi << endl; 
	}

	// shut down the communicator
	MPI_Finalize();

	return EXIT_SUCCESS;
}