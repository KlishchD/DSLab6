#include "header.h"

#if ALGORITHM == FOX_ALGORITHM

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

int main(int argc, char** argv) {
    int comm_sz;        
    int my_rank;        
    int n = MATRIX_SIZE;
    int master = 0;     
    int tag = 0;        
    float* a = new float[n * n];     
    float* b = new float[n * n];     
    float* matrix = new float[n * n];
	double start_time;
	double finish_time;
	double final_time;

    srand(0);
    float num;
    for (int i = 0; i < n * n; i++) {
        num = (float)rand() / RAND_MAX * 2.0 - 1.0;
        a[i] = num;
        num = (float)rand() / RAND_MAX * 2.0 - 1.0;
        b[i] = num;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	printf("My rank %d \n", my_rank);

	if (my_rank == 0) {
		printarr(a, n);          
		printarr(b, n);
	}

    int np = pow(comm_sz, 0.5);
    int nr = n / np;           
    MPI_Request request;       
    MPI_Status status;         

    if (n % np != 0) {         
        MPI_Finalize();
        fprintf(stderr, "Square Root of Processes does not divide Number of elements.\n");

		delete[] a;
		delete[] b;
		delete[] matrix;

        return 0;
    }

    start_time = MPI_Wtime();

    if (my_rank == master) {

        fprintf(stdout, "\nSize: %d\nProcesses: %d\n", n, comm_sz);

        int i = 0;
        int k = 0;
        for (int row_num = 0; row_num < np; row_num++) {
            for (int col_num = 0; col_num < np; col_num++) {
                k = col_num * nr + row_num * n * nr;
                for (int j = 0; j < nr * nr; j++) {
                    if (j % nr == 0 && j != 0) {
                        k = k + (n - nr);
                    }

                    matrix[i] = a[k];

                    k++;
                    i++;
                }
            }
        }

        i = 0;
        k = 0;
        for (int row_num = 0; row_num < np; row_num++) {
            for (int col_num = 0; col_num < np; col_num++) {
                k = col_num * nr + row_num * n * nr;
                for (int j = 0; j < nr * nr; j++) {
                    if (j % nr == 0 && j != 0) {
                        k = k + (n - nr);
                    }

                    a[i] = b[k];

                    k++;
                    i++;
                }
            }
        }
    }

    if (my_rank == 0) {
        printf("Matrix");
        printarr(matrix, n);    
        printf("B");
        printarr(b, n);
    }

    float* rank_a = new float[nr * nr];
    float* rank_b = new float[nr * nr];
    float* local_a = new float[nr * nr];
    float* local_b = new float[nr * nr];

    MPI_Scatter(a, (nr * nr), MPI_FLOAT, rank_a, (nr * nr), MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, (nr * nr), MPI_FLOAT, rank_b, (nr * nr), MPI_FLOAT, 0, MPI_COMM_WORLD);

    float* result = new float[nr * nr];
    for (int i = 0; i < nr * nr; i++) {
        result[i] = 0.0;               
    }

    int* source = new int[np * np];
    for (int i = 0; i < np; ++i) {
        source[i] = i * (np + 1);
    }
    for (int i = 1; i < np; ++i) {
        for (int j = 0; j < np; ++j) {
            if ((source[(i - 1) * np + j] + 1) >= np * (j + 1)) {
                source[i * np + j] = np * (j);
            }
            else {
                source[i * np + j] = source[(i - 1) * np + j] + 1;
            }
        }
    }

    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {

            for (int k = 0; k < nr * nr; k++) {  
                local_b[k] = rank_b[k];          
            }

            int low = (source[i * np + j] / np);
            int high = (source[i * np + j] / np);
            low = low * np;
            high = (high + 1) * np;

            if (my_rank == source[i * np + j]) { 
                for (int k = low; k < high; k++) {
                    if (my_rank != k) {
                        MPI_Send(rank_a, nr * nr, MPI_FLOAT, k, tag, MPI_COMM_WORLD);
                    }
                    else { 
                        for (int k = 0; k < nr * nr; k++) {
                            local_a[k] = rank_a[k];
                        }
                    }
                }
            }
            else {
                if (my_rank >= low && my_rank < high) {
                    MPI_Recv(local_a, nr * nr, MPI_FLOAT, source[i * np + j], tag, MPI_COMM_WORLD, &status);
                }
            }

            for (int x = 0; x < nr; x++) {
                for (int y = 0; y < nr; y++) {
                    for (int z = 0; z < nr; z++) {
                        result[x * nr + y] = result[x * nr + y] + rank_a[x * nr + z] * local_b[z * nr + y];
                    }
                }
            }

            int destination = my_rank - np;
            int source = my_rank + np;
            if (destination < 0) {
                destination = np * np + destination;
            }
            if (source >= np * np) {
                source = source - np * np;
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Isend(rank_b, nr * nr, MPI_FLOAT, destination, 0, MPI_COMM_WORLD, &request);
            MPI_Irecv(rank_b, nr * nr, MPI_FLOAT, source, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    finish_time = MPI_Wtime();
    final_time = finish_time - start_time;

    if (my_rank == master) {
        fprintf(stdout, "\n");
        fprintf(stdout, "Total Time Elapsed is %.10f seconds\n", final_time);
    }

    int greeting = 1;
    if (my_rank == master) {
        printarr(result, nr);
        MPI_Send(&greeting, 1, MPI_INT, my_rank + 1, tag, MPI_COMM_WORLD);
    }
    else if (my_rank == comm_sz - 1) {
        MPI_Recv(&greeting, 1, MPI_INT, my_rank - 1, tag, MPI_COMM_WORLD, &status);
        printarr(result, nr);
    }
    else {
        MPI_Recv(&greeting, 1, MPI_INT, my_rank - 1, tag, MPI_COMM_WORLD, &status);
        printarr(result, nr);
        MPI_Send(&greeting, 1, MPI_INT, my_rank + 1, tag, MPI_COMM_WORLD);
    }

	delete[] a;
	delete[] b;
	delete[] matrix;

    delete[] rank_a;
	delete[] rank_b;
	delete[] local_a;
	delete[] local_b;

    delete[] result;
    delete[] source;

    MPI_Finalize();

    return 0;
}

#endif 