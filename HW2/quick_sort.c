#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <pthread.h>
#include <omp.h>

void Insertion_sort(int* A, int q, int r){
	
	if (r <= q) return;
	//printf("In serial computation\n");	
	int i, key, j;
	int n = r - q + 1;
    for (i = 1; i < n; i++) {
        key = A[q + i];
        j = i + q - 1;
 
        /* Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
        while (j >= q && A[j] > key) {
            A[j + 1] = A[j];
            j = j - 1;
        }
        A[j + 1] = key;
    }
}


int* Prefix_sum(int* input_array, int len){
	
	int* output_array = malloc(len * sizeof(int));
	
	if (len == 1){
		output_array[0] = input_array[0];
	}
	else{
		int new_len = len/2;
		int* new_array = malloc(new_len * sizeof(int));

		#pragma omp parallel for
		for (int i = 0; i < new_len; i++){
			new_array[i] = input_array[2*(i+1) - 1] + input_array[2*i];
		}

		int* new_array2 = Prefix_sum(new_array, new_len);

		#pragma omp parallel for
		for (int i = 0; i < len; i++){
			if (i == 0){
				output_array[0] = input_array[0];
			}
			else if (i % 2 == 1){
				output_array[i] = new_array2[(i - 1)/2];
			}
			else{
				output_array[i] = new_array2[(i - 2)/2] + input_array[i];
 			}
		}	
	}
	
	return output_array;	
}

int Par_Partition(int* A, int q, int r, int x){

    int n = r - q + 1;
    if (n == 1) return q;

    int* B = malloc(n * sizeof(int));
    int* lt = malloc(n * sizeof(int));
    int* gt = malloc(n * sizeof(int));

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        B[i] = A[q + i];
        if (B[i] < x)
			lt[i] = 1;
		else	
			lt[i] = 0;

		if (B[i] > x)
			gt[i] = 1;
		else
			gt[i] = 0;
    }

	lt = Prefix_sum(lt, n);
	gt = Prefix_sum(gt, n);

	// for (int i = 0; i < n; i++){
	// 	printf(" %d", lt[i]);
	// }
	// printf("\n");

	// for (int i = 0; i < n; i++){
	// 	printf(" %d", gt[i]);
	// }
	// printf("\n");


	int k = q + lt[n-1];
	A[k] = x;

	#pragma omp parallel for
	for (int i = 0; i < n; i++){
		if (B[i] < x)
			A[q + lt[i] - 1] = B[i];
		else if (B[i] > x)
			A[k + gt[i]] = B[i];
	}

	return k;
}



void Par_Randomized_Quicksort(int* A, int q, int r){

	if (r < q) return;

	int n = r - q + 1;

	//printf("q and r index is %d %d\n", q, r);

	if (n <= 4096){
		//call insertion sort
		Insertion_sort(A, q, r);
	}
	else{
		int indexofx = (rand() % n) + q;
		int x = A[indexofx];

		//printf("index of x is %d\n", indexofx);
		//printf("x is %d\n", x);

		int k = Par_Partition(A, q, r, x);

		//printf("pivot index k is %d\n", k);

		// for (int i = q; i < n; i++){
        // 	printf(" %d", A[q + i]);
    	// }

    	//printf("\n");

		cilk_spawn Par_Randomized_Quicksort(A, q, k - 1);

		Par_Randomized_Quicksort(A, k + 1, r);

		cilk_sync;

	}
}


int main()
{
	__cilkrts_set_param("nworkers", "68");
	int numth = __cilkrts_get_nworkers();
	printf("number of active threads for cilk is %d\n", numth);

	omp_set_num_threads(68);

	printf("number of max threads allowed for openmp is %d\n", omp_get_max_threads());

	for (int startn = 12; startn <= 24; startn++){

	int n = pow(2,startn);

	//matrix A creation
	int* Input_Array = malloc(n * sizeof(int));
	

	for (int i = 0; i < n; i++){
		//Input_Array[i] = rand()%10;
		Input_Array[i] = n - i;
	}

    // for (int i = 0; i < n; i++){
    //     printf(" %d", Input_Array[i]);
    // }

    // printf("\n");

	
	struct timeval begin, end;
    gettimeofday(&begin, 0);
   
   /* 
    pthread_t threadid[numth];
	int indexes[numth];
    for (int i=0; i<numth; i++) {
		indexes[i] = i;
		
        pthread_create(&threadid[i], NULL, thread_func, (void*)(indexes + i));
    }
	
    for (int i = 0; i < numth; i++)
       pthread_join(threadid[i], NULL); */

    //int* ouput_array = Prefix_sum(Input_Array, n);

	//int k = Par_Partition(Input_Array, 0, n - 1, 3);

	//Insertion_sort(Input_Array, 4, 7);

	Par_Randomized_Quicksort(Input_Array, 0, n - 1);
	
	gettimeofday(&end, 0);

	// for (int i = 0; i < n; i++){
    //     printf(" %d", Input_Array[i]);
    // }

    // printf("\n");

	free(Input_Array);

	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double totaltime = seconds + microseconds*1E-6;
	printf("n val is %d and Time taken for parallel quick sort is %.20f\n", startn, totaltime);
	}

}