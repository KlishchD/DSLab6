#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <chrono>

#define MATRIX_SIZE 3000

#define FOX_ALGORITHM 0
#define SERIAL_ALGORITHM 1

#define ALGORITHM 1

void printarr(float* arr, int n) {
	fprintf(stdout, "\n");
	for (int row = 0; row < n * n; row++) {

		if (row % n == 0 && row != 0) {
			fprintf(stdout, "\n");
		}

		fprintf(stdout, "%6.2f ", arr[row]);
	}

	fprintf(stdout, "\n\n");

	return;
}
