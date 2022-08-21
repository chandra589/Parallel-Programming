#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <omp.h>
#include <stdbool.h>

#define NumofProcessors 1

typedef struct edge{
    int u;
    int v;
    double wt;
}Edge;

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

void Insertion_sort_Edges(Edge* A, int q, int r){
	
	if (r <= q) return;
	//printf("In serial computation\n");	
	int i, j;
    double keywt;
    int keyu, keyv;
	int n = r - q + 1;
    for (i = 1; i < n; i++) {
        keywt = A[q + i].wt;
        keyu = A[q + i].u;
        keyv = A[q + i].v;

        j = i + q - 1;
 
        /* Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
        while (j >= q && A[j].wt > keywt) {
            A[j + 1].wt = A[j].wt;
            A[j + 1].u = A[j].u;
            A[j + 1].v = A[j].v;
            j = j - 1;
        }
        A[j + 1].wt = keywt;
        A[j + 1].u = keyu;
        A[j + 1].v = keyv;
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


void Priority_Con_write(int n, Edge* E, int numofEdges, int* R){

    int* A = malloc(numofEdges * sizeof(int));
    double dum = ceil(log(numofEdges)) + 1;
    int k = (int)(dum);

    //#pragma omp parallel for
    for (int i = 1; i <= numofEdges; i++){
        // if (E[i - 1].u == E[i - 1].v)
        //     A[i - 1] = 0;
        // else
        //     A[i - 1] = (E[i - 1].u << k) + i;
    
        A[i - 1] = (E[i - 1].u << k) + i;

        printf("%d ",A[i - 1]);
    }

    printf("\n");
    printf("Step - 3\n");

    //
    dum = ceil(log(n));
    int dum1 = (int)(dum);

    printf("num of bits to radix sort is %d\n", k + dum1);

    Radix_Sort(A, numofEdges, k + dum1 + 1);

    for (int i = 1; i <= numofEdges; i++){
        printf("%d ",A[i - 1]);
    }

    printf("\n");
    printf("Step - 4\n");

    #pragma omp parallel for
    for (int i = 1; i <= numofEdges; i++){
        int u = A[i - 1] >> k;
        int j = A[i - 1] - (u << k);
        if ((i == 1) || (u != (A[i - 2] >> k)))
            R[u - 1] = j; 
    }

}

void MST_Priority_CW(int n, Edge* E, int numofEdges, int* MST){

    int* L = malloc(n * sizeof(int));
    int* C = malloc(n * sizeof(int));
    int* R = malloc(n * sizeof(int));

    //sort the edges based on edge weights 
    Insertion_sort_Edges(E, 0, numofEdges - 1);

    printf("Sorted Edges\n");
    for (int i = 0; i < numofEdges; i++){
        printf("%d %d %f\n", E[i].u, E[i].v, E[i].wt);
    }

    #pragma omp parallel for
    for (int v = 0; v < n; v++){
        L[v] = v + 1;
    }

    bool F = (numofEdges > 0) ? true : false;

    //int numoflevels = numofEdges;
    //int prenumlevels = numoflevels - 1;

    while (F){
        #pragma omp parallel for
        for (int v = 0; v < n; v++){
            C[v] = rand()%2;  //Head - 0 and Tail - 1
        }

        printf("Step -1\n");

        for (int v = 0; v < n; v++){
            printf("%d ", C[v]);
        }
        printf("\n");

        //Priority-Write
        Priority_Con_write(n, E, numofEdges, R);


        for (int i = 0; i < n; i++){
            printf("%d\n", R[i]);
        }

        printf("Step-2\n");

        //#pragma omp parallel for
        for (int i = 1; i <= numofEdges; i++){
            int u = E[i - 1].u;
            int v = E[i - 1].v;
            if (C[u - 1] == 1 && C[v - 1] == 0 && R[u - 1] == i){
                L[u-1] = v;
                MST[i-1] = 1;
            }
        }


        printf("L vector\n");
        for (int i = 0; i < n; i++){
            printf("%d ", L[i]);
        }
        printf("\n");

        printf("MST vector\n");
        for (int i = 0; i < numofEdges; i++){
            printf("%d ", MST[i]);
        }
        printf("\n");

        //#pragma omp parallel for
        for (int i = 0; i < numofEdges; i++){
            E[i].u = L[E[i].u - 1];
            E[i].v = L[E[i].v - 1];
        }


        printf("new levels\n");
        for (int i = 0; i < numofEdges; i++){
            printf("%d %d\n", E[i].u, E[i].v);
        }
        printf("\n");

        exit(0);
        F = false;

        // numoflevels = 0;
        // #pragma omp parallel for shared(numoflevels, E) reduction(+: numoflevels)
        // for (int i = 0; i < numofEdges; i++){
        //     if (E[i].u != E[i].v)
        //         numoflevels += 1;
        // }

        #pragma omp parallel for
        for (int i = 0; i < numofEdges; i++){
            if (E[i].u != E[i].v){
                F = true;
            }
        }

    }

}

int main()
{   
	// __cilkrts_set_param("nworkers", "4");
	// int numth = __cilkrts_get_nworkers();
	// printf("number of active threads is %d\n", numth);

    omp_set_num_threads(NumofProcessors);
    printf("number of max threads allowed is %d\n", omp_get_max_threads());

    FILE* f;
    if((f = fopen("test.txt", "r")) == NULL)
       exit(1);

    int n, numofEdges;

    if(fscanf(f, "%d %d", &n, &numofEdges) != 2)
        exit(1);

    Edge* edges = malloc(numofEdges * sizeof(Edge));

    int u, v;
    double w =0;
    for (int i = 0; i < numofEdges; i++){
        if(fscanf(f, "%d %d %lf", &u, &v, &w) != 3)
            exit(1);
        edges[i].u = u; edges[i].v = v; edges[i].wt = w;
    }

    fclose(f);
	
    printf("vertices and edges are %d and %d\n", n, numofEdges);

    for (int i = 0; i < numofEdges; i++){
        printf("%d %d %f\n", edges[i].u, edges[i].v, edges[i].wt);
    }


    int* MST = malloc(numofEdges * sizeof(int));
    for (int i = 0; i < n; i++){
        MST[i] = 0;
    }
	
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

    MST_Priority_CW(n, edges, numofEdges, MST);
	
	gettimeofday(&end, 0);

    for (int i = 0; i < numofEdges; i++){
        printf(" %d", MST[i]);
    }

    printf("\n");

	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double totaltime = seconds + microseconds*1E-6;
	printf("Time taken for finding MST is %.20f\n", totaltime);

}