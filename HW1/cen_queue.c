#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <pthread.h>

int arraysize;

float** matrix_A = NULL;
float** matrix_B = NULL;
float** matrix_C = NULL;

void serial_matrix_computation(int R_indexA, int C_indexA, int R_indexB, int C_indexB, int n){
	
	//printf("In serial computation\n");	
	for (int i = R_indexA; i < (R_indexA + n); i++){
		for (int k = C_indexA; k < (C_indexA + n); k++){
			for (int j = C_indexB; j < (C_indexB + n); j++){
				matrix_C[i][j] = matrix_C[i][j] + (matrix_A[i][k]*matrix_B[k][j]);
			}
		}
	}
}

typedef struct job{
    int R_indexA;
    int C_indexA;
    int R_indexB;
    int C_indexB;
    int n;
    struct job *next;
    //struct job *prev;
}Job;

typedef struct workqueue{
    int numjobs;
    pthread_mutex_t mutex; 
    Job *qbottom; //head
    //Job *qtop; //tail
}CQUEUE;

#define num_threads 68

CQUEUE *dequeue = NULL;
pthread_barrier_t barr;

void AddjobtoQueue(int ra, int ca, int rb, int cb, int n){
    if (pthread_mutex_lock(&dequeue->mutex) < 0) return;

    dequeue->numjobs += 1;
    //printf("current job is %d and size is %d\n", dequeue->numjobs, n);

    if (dequeue->qbottom == NULL)
	{
		dequeue->qbottom = malloc(sizeof(Job));
		dequeue->qbottom->R_indexA = ra;
        dequeue->qbottom->C_indexA = ca;
        dequeue->qbottom->R_indexB = rb;
        dequeue->qbottom->C_indexB = cb;
        dequeue->qbottom->n = n;
		dequeue->qbottom->next = NULL;
	}
	else
	{
        Job* newnode = malloc(sizeof(Job));
        newnode->R_indexA = ra;
        newnode->C_indexA = ca;
        newnode->R_indexB = rb;
        newnode->C_indexB = cb;
        newnode->n = n;
        newnode->next = dequeue->qbottom;
        dequeue->qbottom = newnode;

	}

    if (pthread_mutex_unlock(&dequeue->mutex) < 0) return;
}


Job* GetJobfromQueue(){

    if (pthread_mutex_lock(&dequeue->mutex) < 0) return NULL;
    if (dequeue->qbottom == NULL){
        if (pthread_mutex_unlock(&dequeue->mutex) < 0) return NULL;
        return NULL;
    }

    Job* retnode = dequeue->qbottom;
    dequeue->qbottom = dequeue->qbottom->next;

    if (pthread_mutex_unlock(&dequeue->mutex) < 0) return NULL;
    return retnode;

}

void *thread_func(void* arg){
    //pthread_detach(pthread_self());
	int numcheck = 2 * num_threads * log(num_threads + 1);
    int count = 0;
    int res;
    while(1){
        
        Job* item = GetJobfromQueue();
        if (item == NULL){
            if (count == numcheck){
                res = pthread_barrier_wait(&barr);
                break;
            }
            count++;
            continue;
        }
		
		count = 0;


        int ra = item->R_indexA;
        int ca = item->C_indexA;
        int rb = item->R_indexB;
        int cb = item->C_indexB;
        int n = item->n;

        //printf("ra is %d, ca is %d, rb is %d, cb is %d\n", ra, ca, rb, cb);

        if (n <= 32)
            serial_matrix_computation(ra, ca, rb, cb, n);
        else{
            int size = n/2;
            AddjobtoQueue(ra, ca, rb, cb, size);
            AddjobtoQueue(ra, ca, rb, cb + size, size);
            AddjobtoQueue(ra + size, ca, rb, cb, size);
            AddjobtoQueue(ra + size, ca, rb, cb + size, size);
        }
    }

    count = 0;
    if (res == PTHREAD_BARRIER_SERIAL_THREAD){
        //serial thread add remaining jobs
        int size = arraysize/2;
        AddjobtoQueue(0, 0 + size, 0 + size, 0, size);
        AddjobtoQueue(0, 0 + size, 0 + size, 0 + size, size);
        AddjobtoQueue(0 + size, 0 + size, 0 + size, 0, size);
        AddjobtoQueue(0 + size, 0 + size, 0 + size, 0 + size, size);
    }
	
	while(1){
	
		Job* item = GetJobfromQueue();
		if (item == NULL){
			if (count == numcheck){
				break;
			}
			count++;
			continue;
		}
		
		count = 0;


		int ra = item->R_indexA;
		int ca = item->C_indexA;
		int rb = item->R_indexB;
		int cb = item->C_indexB;
		int n = item->n;

		//printf("ra is %d, ca is %d, rb is %d, cb is %d\n", ra, ca, rb, cb);

		if (n <= 32)
			serial_matrix_computation(ra, ca, rb, cb, n);
		else{
			int size = n/2;
			AddjobtoQueue(ra, ca + size, rb + size, cb, size);
			AddjobtoQueue(ra, ca + size, rb + size, cb + size, size);
			AddjobtoQueue(ra + size, ca + size, rb + size, cb, size);
			AddjobtoQueue(ra + size, ca + size, rb + size, cb + size, size);
		}
	}

    return NULL;
}


int main()
{
	int n = pow(2, 12);
    arraysize = n;
    int numth = num_threads;
	//__cilkrts_set_param("nworkers", "32");
	//int numth = __cilkrts_get_nworkers();
	//printf("number of active threads is %d\n", numth);

	//matrix A creation
	matrix_A = malloc(n * sizeof(float*));
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
	matrix_B = malloc(n * sizeof(float*));
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
	matrix_C = malloc(n * sizeof(float*));
	for (int i = 0; i < n; i++){
		matrix_C[i] = malloc(n * sizeof(float));
	}

	//matrix C initialization
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			matrix_C[i][j] = 0;
		}
	}

    dequeue = malloc(sizeof(CQUEUE));
    dequeue->numjobs = 0;
    dequeue->qbottom = NULL;

    int ret = pthread_mutex_init(&dequeue->mutex, NULL);
    if (ret < 0) return 1;

    pthread_barrier_init(&barr, NULL, numth); 

    int size = n/2;
    AddjobtoQueue(0, 0, 0, 0, size);
    AddjobtoQueue(0, 0, 0, 0 + size, size);
    AddjobtoQueue(0 + size, 0, 0, 0, size);
    AddjobtoQueue(0 + size, 0, 0, 0 + size, size);

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
   
    pthread_t threadid[numth];
    for (int i=0; i<numth; i++) {
        pthread_create(&threadid[i], NULL, thread_func, NULL);
    }
	
    for (int i = 0; i < numth; i++)
       pthread_join(threadid[i], NULL);

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

    long totalfpo = 2 * pow(n, 3);
	double flops = totalfpo/totaltime;
	double gflops = flops/pow(10, 9);
	printf("GFLOPS is %f", gflops);
}