#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

void serial_matrix_computation(float** matrix_A, float** matrix_B, float** matrix_C, int R_indexA, int C_indexA, int R_indexB, int C_indexB, int n){
	
	//printf("In serial computation\n");	
	for (int i = R_indexA; i < (R_indexA + n); i++){
		for (int k = C_indexA; k < (C_indexA + n); k++){
			for (int j = C_indexB; j < (C_indexB + n); j++){
				matrix_C[i][j] = matrix_C[i][j] + (matrix_A[i][k]*matrix_B[k][j]);
			}
		}
	}
}


void parallel_recursive_multiplication(float** matrix_A, float** matrix_B, float** matrix_C, int R_indexA, int C_indexA, int R_indexB, int C_indexB, int n){
	if (n <= 32){
		serial_matrix_computation(matrix_A, matrix_B, matrix_C, R_indexA, C_indexA, R_indexB, C_indexB, n);
	}
	else{
		//printf("In else block\n");
		int size = n/2;
		
		cilk_spawn parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA, C_indexA, R_indexB, C_indexB, size);
		
		cilk_spawn parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA, C_indexA, R_indexB, C_indexB + size, size);
		
		cilk_spawn parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA + size, C_indexA, R_indexB, C_indexB, size);
		
		parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA + size, C_indexA, R_indexB, C_indexB + size, size);
		
		cilk_sync;
		
		cilk_spawn parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA, C_indexA + size, R_indexB + size, C_indexB, size);
		
		cilk_spawn parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA, C_indexA + size, R_indexB + size, C_indexB + size, size);
		
		cilk_spawn parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA + size, C_indexA + size, R_indexB + size, C_indexB, size);
		
		parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, R_indexA + size, C_indexA + size, R_indexB + size, C_indexB + size, size);
		
		cilk_sync;
		
	}
}



int main()
{
	int n = pow(2, 10);
	//__cilkrts_set_param("nworkers", "32");
	int numth = __cilkrts_get_nworkers();
	printf("number of active threads is %d\n", numth);

	//matrix A creation
	float** matrix_A = malloc(n * sizeof(float*));
	for (int i = 0; i < n; i++){
		matrix_A[i] = malloc(n * sizeof(float));
	}

	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			float x = rand()%10;
			matrix_A[i][j] = x/3;
		}
	}

	//matrix B creation
	float** matrix_B = malloc(n * sizeof(float*));
	for (int i = 0; i < n; i++){
		matrix_B[i] = malloc(n * sizeof(float));
	}

	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			float x = rand()%10;
			matrix_B[i][j] = x/3;
		}
	}


	//matrix C creation
	float** matrix_C = malloc(n * sizeof(float*));
	for (int i = 0; i < n; i++){
		matrix_C[i] = malloc(n * sizeof(float));
	}

	//matrix C initialization
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			matrix_C[i][j] = 0;
		}
	}


	//papi events
	int retval;
	#define num_events 1
	int eventset = PAPI_NULL;
	long long values[num_events];
	int PAPI_events[num_events] = {PAPI_L2_TCM};

	if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT){
			printf("PAPI_version_mismatch");
			exit(1);
	}

	values[0] = 0;
	//values[1] = 0;

	retval = PAPI_query_event(PAPI_L1_TCM);
	if (retval != PAPI_OK)
			printf("L1_TCM not present\n");

	retval = PAPI_query_event(PAPI_L2_TCM);
	if (retval != PAPI_OK)
			printf("L2_TCM not present\n");


	if ((retval = PAPI_create_eventset(&eventset)) != PAPI_OK){
			printf("PAPI_event_creation_error");
			exit(1);
	}

	if ((retval = PAPI_add_events(eventset, PAPI_events, num_events)) != PAPI_OK){
			printf("PAPI_event_addition_error\n");
			exit(1);
	}

	if ((retval = PAPI_start(eventset)) != PAPI_OK){
			printf("PAPI_start_error");
			exit(1);
	}

	//matrix_A*matrix_B multiplication
	struct timeval begin, end;
        gettimeofday(&begin, 0);

	parallel_recursive_multiplication(matrix_A, matrix_B, matrix_C, 0, 0, 0, 0, n);
	
	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double totaltime = seconds + microseconds*1E-6;
	printf("num of threads is %d and Time taken for matrix multiplication is %.20f\n", numth, totaltime);

	if ((retval = PAPI_stop(eventset, values)) != PAPI_OK){
			printf("PAPI_stop_error");
			exit(1);
	}

	printf("L1orL2 Cache Misses is %lld\n", values[0]);
	//printf("L2 Cache Misses is %lld\n", values[1]);
}
