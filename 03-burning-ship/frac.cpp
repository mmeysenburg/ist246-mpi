#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include <mpi.h>
#include "EasyBMP.h"

/*
* Program to create a Burning Ship fractal, using MPI.
* This program uses the EasyBMP C++ Bitmap Library
* (see http://easybmp.sourceforge.net/).
*/


/**
 * Function to map from pixel coordinates to fractacl coordinates. 
 *
 * \param pixelCoord Pixel coord value (row or column) to map to fractal coord
 *
 * \param maxPixelCoord Maximum value of pixelCoord; min is assumed to be 0
 *
 * \param fracRangeLow Low value in fractal coord range
 *
 * \param fracRangeHigh High value in fractal coord range
 *
 * \return pixelCoord mapped to the range [fracRangeLow, fracRangeHigh]
 */
inline double mapToRange(double pixelCoord, double maxPixelCoord, 
	double fracRangeLow, double fracRangeHigh) {

    return fracRangeLow + (fracRangeHigh - fracRangeLow) /
		maxPixelCoord * pixelCoord;
}

/**
 * Test a point to see if it is in the Burning Ship fractal or not.
 *
 * \param x X value of point to be tested
 *
 * \param y Y value of point to be tested
 * 
 * \return Number of iterations before fractal value gets too large,
 * in [0, 255]
 */
int testPoint(double x, double y) {
    double a = 0.0, b = 0.0, ta = 0.0;
    int i = 0;

    while(i < 255) {
        a = a < 0.0 ? -a : a;
        b = b < 0.0 ? -b : b;

        ta = (a * a) - (b * b) + x;
        b = 2 * a * b + y;
        a = ta;
        
        if (a < -2.0 || a > 2.0 || b < -2.0 || b > 2.0) {
            return i;
        }

        i++;
    }

    return i;
}

/**
 * Application entry point. Usage, to produce a fractal image of size 
 * <width> x <height>, looking at the box bounded by (<x0>, <yo>) and 
 * (<x1>, <y1>):
 *
 * mpirun bin/frac.mpi <width> <height> <x0> <x1> <y0> <y1>
 *
 * If no parameters are supplied, or the wrong number of parameters, the app
 * defaults to:
 *
 * <width> = 1920, <height> = 1080
 *
 * <x0> = -2.0, <x1> = 1.5, <y0> = -2.0, <y1> = 0.5,
 *
 * which produces the classic Burning Ship image. 
 *
 * \param argC number of command-line arguments
 *
 * \param ppArgv command-line arguments, as an array of C strings
 *
 * \return EXIT_SUCCESS if everthing went well
 */
int main(int argC, char** ppArgv) {
    using namespace std;

	// Tags for MPI send and receive operations
	int MSG_IMG_SIZE = 0;        // image width and height
	int MSG_FRACTAL_COORDS = 1;  // fractal box bounding coordinates
	int MSG_SLICE_COORDS = 2;    // bonding rows for sub-processor slices
	int MSG_SLICE = 3;           // points in sub-processor's slice

	int rank;    // this processor's rank in the communicator
	int nProcs;  // number of processors in the communicator

	// initialize MPI constructs
	MPI_Init(&argC, &ppArgv);                // set up the communicator
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);    // which processor am I?
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);  // how many processors are there?

	// arrays to hold size of image and bounding coordinates of the fractal;
	// initialized to hold the default values
	int pImgSize[2] = { 1920, 1080 };
	double pFracCoords[4] = { -2.0, 1.5, -2.0, 0.5 };

	// buffer for sub-processor slice row start and finish
	int pSliceBounds[2];

	//-------------------------------------------------------------------------
	// root process algorithm: send image size, fractal bounds, and slice info
	// to sub-processors, receive their slices, assemble them into a bitmap 
	// image, and then output the image
	//-------------------------------------------------------------------------
	if (rank == 0) {

		// if the correct number of command-line arguments are provided, 
		// use them for image size and fractal bounds; otherwise, the defaults
		// are used
		if (argC == 7) {
			pImgSize[0] = atoi(ppArgv[1]);
			pImgSize[1] = atoi(ppArgv[2]);
			pFracCoords[0] = atof(ppArgv[3]);
			pFracCoords[1] = atof(ppArgv[4]);
			pFracCoords[2] = atof(ppArgv[5]);
			pFracCoords[3] = atof(ppArgv[6]);
		}

		// determine the number of rows in the slices and set first slice
		// staring and ending rows
		int sliceHeight = pImgSize[1] / (nProcs - 1);
		pSliceBounds[0] = 0;
		pSliceBounds[1] = sliceHeight - 1;

		// dispatch information to the first nProcs - 2 subprocessors, who will
		// all have identically-sized slices
		for (int proc = 1; proc < nProcs - 1; proc++) {
			// send overall image size and fractal coordinates (in case the 
			// default values were changed via command-line parameters)
			MPI_Send(pImgSize, 2, MPI::INT, proc, MSG_IMG_SIZE, MPI_COMM_WORLD);
			MPI_Send(pFracCoords, 4, MPI::DOUBLE, proc, MSG_FRACTAL_COORDS, 
				MPI_COMM_WORLD);

			// send slice starting and ending row values
			MPI_Send(pSliceBounds, 2, MPI::INT, proc, MSG_SLICE_COORDS, 
				MPI_COMM_WORLD);

			// update slice starting and ending row values for next processor
			pSliceBounds[0] = pSliceBounds[1] + 1;
			pSliceBounds[1] += sliceHeight;
		}

		// the final sub-processor's slice may be of a different size, so send its 
		// values; ending row is just last row in the image
		pSliceBounds[1] = pImgSize[1] - 1;

		MPI_Send(pImgSize, 2, MPI::INT, nProcs - 1, MSG_IMG_SIZE, MPI_COMM_WORLD);
		MPI_Send(pFracCoords, 4, MPI::DOUBLE, nProcs - 1, MSG_FRACTAL_COORDS, 
			MPI_COMM_WORLD);
		MPI_Send(pSliceBounds, 2, MPI::INT, nProcs - 1, MSG_SLICE_COORDS, 
			MPI_COMM_WORLD);

		// while the sub-processors work, create a blank bitmap image
		BMP image;
		image.SetSize(pImgSize[0], pImgSize[1]);
		image.SetBitDepth(32);

		// now, collect data from the sub-processors and create the image

		// first, create a buffer to receive data, big enough to hold the 
		// largest slice we might receive
		int lastSliceHeight = pSliceBounds[1] - pSliceBounds[0] + 1;
		int maxNumPoints = (lastSliceHeight > sliceHeight) ?
			lastSliceHeight : sliceHeight;
		maxNumPoints *= pImgSize[0];

		int* pPoints = new int[maxNumPoints];
		int numPoints;

		// p represents a single pixel in the image; its color values will be
		// adjusted based on the values computed with the testPoint() 
		// function
		RGBApixel p;
		p.Green = (ebmpBYTE)0;  // in this coloring scheme, we will only change
		p.Blue = (ebmpBYTE)0;   // the red componend of the color
		p.Alpha = (ebmpBYTE)0;

		MPI_Status status;     // used to determine number of points we get
		int col = 0, row = 0;  // image pixel coordinage

		// receive calculated points from subprocessors, from top to bottom
		// of the image
		for (int proc = 1; proc < nProcs; proc++) {

			MPI_Recv(pPoints, maxNumPoints, MPI::INT, proc, MSG_SLICE, MPI_COMM_WORLD, 
				&status);

			// see how many points we actually received (since the last slice may be 
			// smaller than the rest
			MPI_Get_count(&status, MPI::INT, &numPoints);

			// now use the calculated points to set pixel colors in the output image
			for (int m = 0; m < numPoints; m++) {
				p.Red = (ebmpBYTE)pPoints[m];

				image.SetPixel(col, row, p);

				// update image coordinates
				col++;
				if (col >= pImgSize[0]) {
					col = 0;
					row++;
				}
			}
		}

		// free the receive buffer memory
		delete[] pPoints;

		// finally, write the image to disk
		image.WriteToFile("out.bmp");
	} else {

		//----------------------------------------------------------------------
		// sub-processor algorithm: receive the coordinate data, make a 
		// buffer to hold fractal points, values, compute, and send data back
		//---------------------------------------------------------------------
		MPI_Status status;
		// get image size
		MPI_Recv(pImgSize, 2, MPI::INT, 0, MSG_IMG_SIZE, MPI_COMM_WORLD, 
			&status);
		// get fractal coordinate bounds
		MPI_Recv(pFracCoords, 4, MPI::DOUBLE, 0, MSG_FRACTAL_COORDS, 
			MPI_COMM_WORLD, &status);
		// get first and last row numbers for this processor
		MPI_Recv(pSliceBounds, 2, MPI::INT, 0, MSG_SLICE_COORDS, 
			MPI_COMM_WORLD, &status);

		//cout << rank << ":\t" << pImgSize[0] << "\t" << pImgSize[1] << "\t" <<
		//	pFracCoords[0] << "\t" << pFracCoords[1] << "\t" << pFracCoords[2] << "\t" << pFracCoords[3] << "\t" <<
		//	pSliceBounds[0] << "\t" << pSliceBounds[1] << endl;;

		// create a send buffer large enough to hold all the points this 
		// processor will calculate
		int numPoints = (pSliceBounds[1] - pSliceBounds[0] + 1) * pImgSize[0];
		int *pPoints = new int[numPoints];
		//cout << rank << ":\t" << numPoints << endl;

		double x, y;                // fractal x, y coords
		int col = 0;                // image col coord
		int row = pSliceBounds[0];  // image row coord

		// calculate the points!
		for (int m = 0; m < numPoints; m++) {
			// map image to fractal coordinates
			x = mapToRange(col, pImgSize[0], pFracCoords[0], pFracCoords[1]);
			y = mapToRange(row, pImgSize[1], pFracCoords[2], pFracCoords[3]);

			// get fractal point value
			pPoints[m] = testPoint(x, y);

			// update image coordinates
			col++;
			if (col >= pImgSize[0]) {
				col = 0;
				row++;
			}
		}

		// send stripe of fractal back
		MPI_Send(pPoints, numPoints, MPI::INT, 0, MSG_SLICE, MPI_COMM_WORLD);

		// free this processor's send buffer
		delete[] pPoints;
	}
    
	// shut down MPI and exit the program
	MPI_Finalize();

    return EXIT_SUCCESS;
}