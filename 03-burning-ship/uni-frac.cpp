#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include "EasyBMP.h"

/*
 * Program to create a Burning Ship fractal, using a single process.
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

	while (i < 255) {
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
 * bin/uni-frac <width> <height> <x0> <x1> <y0> <y1>
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

	// arrays to hold size of image and bounding coordinates of the fractal;
	// initialized to hold the default values
	int pImgSize[2] = { 1920, 1080 };
	double pFracCoords[4] = { -2.0, 1.5, -2.0, 0.5 };

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

	// create a blank bitmap image
	BMP image;
	image.SetSize(pImgSize[0], pImgSize[1]);
	image.SetBitDepth(32);

	RGBApixel p;
	p.Green = (ebmpBYTE)0;  // in this coloring scheme, we will only change
	p.Blue = (ebmpBYTE)0;   // the red componend of the color
	p.Alpha = (ebmpBYTE)0;
	double x, y;
	for (int row = 0; row < pImgSize[1]; row++) {
		for (int col = 0; col < pImgSize[0]; col++) {
			// map image to fractal coordinates
			x = mapToRange(col, pImgSize[0], pFracCoords[0], pFracCoords[1]);
			y = mapToRange(row, pImgSize[1], pFracCoords[2], pFracCoords[3]);

			p.Red = (ebmpBYTE)testPoint(x, y);
			image.SetPixel(col, row, p);
		}
	}

	// finally, write the image to disk
	image.WriteToFile("out.bmp");

	return EXIT_SUCCESS;
}