#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <time.h>

int main()
{
	int n = pow(2, 10);

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
	
	retval = PAPI_query_event(PAPI_L3_TCM);
	if (retval != PAPI_OK)
		printf("L3_TCM not present\n");


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
	clock_t start = clock();

	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				matrix_C[i][j] = matrix_C[i][j] + (matrix_A[i][k]*matrix_B[k][j]);
			}
		}
	}

	clock_t end = clock();
	float seconds = (float)(end - start)/CLOCKS_PER_SEC;
	printf("Time taken for matrix multiplication is %f\n", seconds);

	if ((retval = PAPI_stop(eventset, values)) != PAPI_OK){
		printf("PAPI_stop_error");
		exit(1);
	}

	printf("L2 Cache Misses is %lld\n", values[0]);
	//printf("L2 Cache Misses is %lld\n", values[1]);

}
