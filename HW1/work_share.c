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
    struct job *prev;
}Job;

typedef struct workqueue{
    int numjobs;
    pthread_mutex_t mutex; 
    Job *head; //head
    Job *tail; //tail
}WORK_QUEUE;

#define num_threads 68
WORK_QUEUE allqueues[num_threads];

pthread_barrier_t barr;


void insertFirst(int ra, int ca, int rb, int cb, int n, int index) {
	
	if (pthread_mutex_lock(&allqueues[index].mutex) < 0) return;

	//create a link
	Job* newnode = malloc(sizeof(Job));
	newnode->R_indexA = ra;
	newnode->C_indexA = ca;
	newnode->R_indexB = rb;
	newnode->C_indexB = cb;
	newnode->n = n;
	newnode->next = NULL;
	newnode->prev = NULL; 

	int numjbs = allqueues[index].numjobs;

	if (numjbs == 0) {
		//make it the last link
		allqueues[index].tail = newnode;
	}
	else {
		//update first prev link
		allqueues[index].head->prev = newnode;
	}

	//point it to old first link
	newnode->next = allqueues[index].head;

	//point first to new first link
	allqueues[index].head = newnode;

	allqueues[index].numjobs += 1;
	
	if (pthread_mutex_unlock(&allqueues[index].mutex) < 0) return;
}


void insertLast(int ra, int ca, int rb, int cb, int n) {

    int queue_num = rand()%num_threads;

    if (pthread_mutex_lock(&allqueues[queue_num].mutex) < 0) return;

    //printf("Inserting in queue %d\n", queue_num);

    //create a link
	Job* newnode = malloc(sizeof(Job));
	newnode->R_indexA = ra;
	newnode->C_indexA = ca;
	newnode->R_indexB = rb;
	newnode->C_indexB = cb;
	newnode->n = n;
	newnode->next = NULL;
	newnode->prev = NULL; 

    int numjbs = allqueues[queue_num].numjobs;
	
    if(numjbs == 0) {
        //make it the last link
        allqueues[queue_num].tail = newnode;
        allqueues[queue_num].head = newnode;
    } else {
        //make link a new last link
        allqueues[queue_num].tail->next = newnode;     
        
        //mark old last node as prev of new link
        newnode->prev = allqueues[queue_num].tail;
    }

   //point last to new last node
   allqueues[queue_num].tail = newnode;

   allqueues[queue_num].numjobs += 1;

   if (pthread_mutex_unlock(&allqueues[queue_num].mutex) < 0) return;
}

Job* deleteFirst(int index) {
	
	if (pthread_mutex_lock(&allqueues[index].mutex) < 0) return NULL;
	
	int numjbs = allqueues[index].numjobs;

    //printf("In deleteFirst - queue is %d and numjobs is %d\n", index, numjbs);
	
    if (numjbs == 0){
        if (pthread_mutex_unlock(&allqueues[index].mutex) < 0) return NULL;
        return NULL;
    }

	//save reference to first link
	Job* tempLink = allqueues[index].head;

	//if only one link
	if (allqueues[index].head->next == NULL){
		allqueues[index].tail = NULL;
	}
	else {
		allqueues[index].head->next->prev = NULL;
	}

	allqueues[index].head = allqueues[index].head->next;

	allqueues[index].numjobs -= 1;
	
	if (pthread_mutex_unlock(&allqueues[index].mutex) < 0) return NULL;
	
	return tempLink;
}

Job* deleteLast(int index) {

	if (pthread_mutex_lock(&allqueues[index].mutex) < 0) return NULL;
	
	int numjbs = allqueues[index].numjobs;
	
    if (numjbs == 0){
        if (pthread_mutex_unlock(&allqueues[index].mutex) < 0) return NULL;
        return NULL;
    }

	//save reference to last link
	Job* tempLink = allqueues[index].tail;

	//if only one link
	if (allqueues[index].head->next == NULL) {
		allqueues[index].head = NULL;
	}
	else {
		allqueues[index].tail->prev->next = NULL;
	}

	allqueues[index].tail = allqueues[index].tail->prev;

	allqueues[index].numjobs -= 1;
	
	if (pthread_mutex_unlock(&allqueues[index].mutex) < 0) return NULL;
	
	return tempLink;
}

void *thread_func(void* arg){
	
	int index = *((int*) arg);
	
    int numcheck = 2 * num_threads * log(num_threads + 1);
	int count = 0;
    int res;
    while(1){
        
        Job* item = deleteFirst(index);
        if (item == NULL){
            count++;
            if (count == numcheck){
                res = pthread_barrier_wait(&barr);
                break;
            }
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
            insertLast(ra, ca, rb, cb, size);
            insertLast(ra, ca, rb, cb + size, size);
            insertLast(ra + size, ca, rb, cb, size);
            insertLast(ra + size, ca, rb, cb + size, size);
        }
    }

    count = 0;
	
    if (res == PTHREAD_BARRIER_SERIAL_THREAD){
        //serial thread add remaining jobs to queue one
        //printf("adding 2nd half of jobs");
        int size = arraysize/2;
        insertFirst(0, 0 + size, 0 + size, 0, size, 0);
        insertFirst(0, 0 + size, 0 + size, 0 + size, size, 0);
        insertFirst(0 + size, 0 + size, 0 + size, 0, size, 0);
        insertFirst(0 + size, 0 + size, 0 + size, 0 + size, size, 0);
    }
	
    //printf("2nd num of jobs is %d\n", allqueues[index].numjobs);
	
	while(1){
	
		Job* item = deleteFirst(index);
		if (item == NULL){
            count++;
            if (count == numcheck){
                break;
            }
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
			insertLast(ra, ca + size, rb + size, cb, size);
			insertLast(ra, ca + size, rb + size, cb + size, size);
			insertLast(ra + size, ca + size, rb + size, cb, size);
			insertLast(ra + size, ca + size, rb + size, cb + size, size);
		}
	}

    return NULL;
}


int main()
{
	int n = pow(2, 15);
    arraysize = n;
    int numth = num_threads;
	
	//Initialize all the queues & the mutex for each queue.
	for (int i = 0; i < numth; i++){
		int ret = pthread_mutex_init(&allqueues[i].mutex, NULL);
		if (ret < 0) return 1;
		allqueues[i].head = NULL;
		allqueues[i].tail = NULL;
		allqueues[i].numjobs = 0;
	}

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

    pthread_barrier_init(&barr, NULL, numth); 

	//Inserting the first batch of jobs in the first queue.
    int size = n/2;
    insertFirst(0, 0, 0, 0, size, 0);
    insertFirst(0, 0, 0, 0 + size, size, 0);
    insertFirst(0 + size, 0, 0, 0, size, 0);
    insertFirst(0 + size, 0, 0, 0 + size, size, 0);

	//papi events
	int retval;
	#define num_events 1
	int eventset = PAPI_NULL;
	long long values[num_events];
	int PAPI_events[num_events] = {PAPI_L1_TCM};

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
	int indexes[numth];
    for (int i=0; i<numth; i++) {
		indexes[i] = i;
        pthread_create(&threadid[i], NULL, thread_func, (void*)(indexes + i));
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