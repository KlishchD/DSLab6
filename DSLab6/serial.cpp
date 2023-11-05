#include "header.h"

#if ALGORITHM == SERIAL_ALGORITHM

int main() {
	int n = MATRIX_SIZE;

	float* a = new float[n * n];
	float* b = new float[n * n];
	float* c = new float[n * n];

	srand(0);
	float num;
	for (int i = 0; i < n * n; i++) {
		num = (float)rand() / RAND_MAX * 2.0 - 1.0;
		a[i] = num;
		num = (float)rand() / RAND_MAX * 2.0 - 1.0;
		b[i] = num;
	}

	printarr(a, n);
	printarr(b, n);

	auto start = std::chrono::system_clock::now();

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double value = 0;
			for (int k = 0; k < n; ++k) {
				value += a[i * n + k] * b[k * n + j];
			}
			c[i * n + j] = value;
		}
	}

	auto end = std::chrono::system_clock::now();

	printarr(c, n);

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s"	<< std::endl;
}

#endif