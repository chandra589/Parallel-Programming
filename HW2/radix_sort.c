#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <omp.h>

#define NumofProcessors 68

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

    //update the input array using the ouput array
    #pragma omp parallel for
    for (int i = 0; i < len; i++){
        input_array[i] = output_array[i];
    }
	
	return output_array;	
}

void Counting_Rank(int* S, int n, int d, int* r){

    int row = pow(2, d);
    int col = NumofProcessors;  //number of processors (p)
    int* f = malloc((row * col)*sizeof(int));
    int* r1 = malloc((row * col)*sizeof(int));
    int* js = malloc(col * sizeof(int));
    int* je = malloc(col * sizeof(int));
    int* ofs = malloc(col * sizeof(int));

    int noverp = n/col; //gives floor 

    #pragma omp parallel for
    for (int i = 0; i < col; i++){
        for (int j = 0; j < row; j++){
            f[j * col + i] = 0;
        }

        js[i] = i * noverp;
        je[i] = (i + 1) * noverp;

        for (int j = js[i]; j < je[i]; j++){
            f[S[j]*col + i] = f[S[j] * col + i] + 1;
        }
    }


    // int* dum = f;
    // printf("before prefix_sum");
    // for (int i = 0; i < row; i++){
    //     for (int j = 0; j < col; j++){
    //         printf("%d", *dum);
    //         dum++;
    //     }
    //     printf("\n");
    // }

    for (int j = 0; j < row; j++){
        Prefix_sum(f + (j*col), col);
    }

    // dum = f;
    // printf("after prefix_sum");
    // for (int i = 0; i < row; i++){
    //     for (int j = 0; j < col; j++){
    //         printf("%d", *dum);
    //         dum++;
    //     }
    //     printf("\n");
    // }

    #pragma omp parallel for
    for (int i = 0; i < col; i++){
        ofs[i] = 1;
        for (int j = 0; j < row; j++){
            r1[j*col + i] = (i == 0) ? ofs[i] : (ofs[i] + f[j*col + i - 1]);
            ofs[i] = ofs[i] + f[j*col + col - 1];
        }

        for (int j = js[i]; j < je[i]; j++){
            r[j] = r1[S[j] * col + i];
            r1[S[j] * col + i] = r1[S[j] * col + i] + 1;
        }
    }
}


void Radix_Sort(int* A, int n, int b){

    int* s = malloc(n * sizeof(int));
    int* r = malloc(n * sizeof(int));
    int* B = malloc(n * sizeof(int));

    int p = NumofProcessors; //number of processors

    double dum = log(n/(p*log(n)));

    //printf("dum val is %f\n", dum);

    double dum1 = ceil(dum);

    //printf("dum1 val is %f\n", dum1);

    int d = (int)(dum1);

    if (d <= 0) d = 1;   //because dum can be negative if n is small and p is high

    //printf("d value is %d\n", d);

    for (int k = 0; k < b; k+=d){

        int q = (k + d <= b) ? d : (b - k);

        //printf("q val and d val is %d and %d\n", q, d);
        
        #pragma omp parallel for
        for (int i = 1; i <= n; i++){
            int bitval = (((1 << q) - 1) & (A[i - 1] >> k));
            //int bitval = (A[i - 1] >> (k + i - 1)) & 1;
            //printf("%d - %d, ", A[i - 1], bitval);
            s[i - 1] = bitval;
        }

        // printf("\n");

        // for (int i = 0; i < n; i++){
        //     printf("%d", s[i]);
        // }

        // printf("\n");

        // printf("in here\n");

        Counting_Rank(s, n, q, r);

        // for (int i = 0; i < n; i++){
        //     printf("%d", r[i]);
        // }

        // printf("\n");
    
        #pragma omp parallel for
        for (int i = 1; i <= n; i++){
            B[r[i - 1] - 1] = A[i - 1];
        }

        #pragma omp parallel for
        for (int i = 1; i <= n; i++){
            A[i - 1] = B[i - 1];
        }

        // for (int i = 0; i < n; i++){
        //     printf("%d", A[i]);
        // }

        // printf("\n");
    }
}


int main()
{
    //int bits = 24;
    
	// __cilkrts_set_param("nworkers", "4");
	// int numth = __cilkrts_get_nworkers();
	// printf("number of active threads is %d\n", numth);

    omp_set_num_threads(NumofProcessors);
    printf("number of max threads allowed is %d\n", omp_get_max_threads());

    for (int bits = 12; bits <= 24; bits++){

    int n = pow(2,bits);

	//matrix A creation
	int* Input_Array = malloc(n * sizeof(int));
	

	for (int i = 0; i < n; i++){
		//Input_Array[i] = rand()%15;
		Input_Array[i] = n - i - 1;
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

	//Par_Randomized_Quicksort(Input_Array, 0, n - 1);

    //int* rank = malloc(n * sizeof(int));
    //Counting_Rank(Input_Array, n, 3, rank);

    Radix_Sort(Input_Array, n, bits);

    // for (int i = 0; i < n; i++){
    //     printf(" %d", Input_Array[i]);
    // }

    // printf("\n");

    free(Input_Array);
	
	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double totaltime = seconds + microseconds*1E-6;
	printf("n value is %d and Time taken for radix sort is %.20f\n", bits, totaltime);
    }

}