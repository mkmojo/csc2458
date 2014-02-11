#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#ifndef _REENTRANT
#define _REENTRANT		/* basic 3-lines for threads */
#endif
#include <pthread.h>

// #define DEBUG
#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}

typedef struct 
{
    int nsize;
    int task_id;
}thread_arg;

int task_num;
struct timeval start, finish;
double **matrix, *X, *R;

/* Pre-set solution. */

double *X__;

/* Barrier for synchronization */
pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
void barrier (int expect)
{
    static int arrived = 0;

    pthread_mutex_lock (&mut);	//lock

    arrived++;
    if (arrived < expect)
        pthread_cond_wait (&cond, &mut);
    else {
        arrived = 0;		// reset the barrier before broadcast is important
        pthread_cond_broadcast (&cond);
    }

    pthread_mutex_unlock (&mut);	//unlock
}
/* Initialize the matirx. */

int initMatrix(const char *fname)
{
    FILE *file;
    int l1, l2, l3;
    double d;
    int nsize;
    int i, j;
    double *tmp;
    char buffer[1024];

    printf("%s\n", fname);
    if ((file = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "The matrix file open error\n");
        exit(-1);
    }
    
    /* Parse the first line to get the matrix size. */
    fgets(buffer, 1024, file);
    sscanf(buffer, "%d %d %d", &l1, &l2, &l3);
    nsize = l1;
#ifdef DEBUG
    fprintf(stdout, "matrix size is %d\n", nsize);
#endif

    /* Initialize the space and set all elements to zero. */
    matrix = (double**)malloc(nsize*sizeof(double*));
    assert(matrix != NULL);
    tmp = (double*)malloc(nsize*nsize*sizeof(double));
    assert(tmp != NULL);    
    for (i = 0; i < nsize; i++) {
        matrix[i] = tmp;
        tmp = tmp + nsize;
    }
    for (i = 0; i < nsize; i++) {
        for (j = 0; j < nsize; j++) {
            matrix[i][j] = 0.0;
        }
    }

    /* Parse the rest of the input file to fill the matrix. */
    for (;;) {
	fgets(buffer, 1024, file);
	sscanf(buffer, "%d %d %lf", &l1, &l2, &d);
	if (l1 == 0) break;

	matrix[l1-1][l2-1] = d;
#ifdef DEBUG
	fprintf(stdout, "row %d column %d of matrix is %e\n", l1-1, l2-1, 
            matrix[l1-1][l2-1]);
#endif
    }

    fclose(file);
    return nsize;
}

/* Initialize the right-hand-side following the pre-set solution. */

void initRHS(int nsize)
{
    int i, j;

    X__ = (double*)malloc(nsize * sizeof(double));
    assert(X__ != NULL);
    for (i = 0; i < nsize; i++) {
        X__[i] = i+1;
    }

    R = (double*)malloc(nsize * sizeof(double));
    assert(R != NULL);
    for (i = 0; i < nsize; i++) {
        R[i] = 0.0;
        for (j = 0; j < nsize; j++) {
            R[i] += matrix[i][j] * X__[j];
        }
    }
}

/* Initialize the results. */

void initResult(int nsize)
{
    int i;

    X = (double*)malloc(nsize * sizeof(double));
    assert(X != NULL);
    for (i = 0; i < nsize; i++) {
        X[i] = 0.0;
    }
}

/* Get the pivot - the element on column with largest absolute value. */

void getPivot(int nsize, int currow)
{
    int i, pivotrow;

    pivotrow = currow;
    for (i = currow+1; i < nsize; i++) {
        if (fabs(matrix[i][currow]) > fabs(matrix[pivotrow][currow])) {
            pivotrow = i;
        }
    }

    if (fabs(matrix[pivotrow][currow]) == 0.0) {
        fprintf(stderr, "The matrix is singular\n");
        exit(-1);
    }
    
    if (pivotrow != currow) {
#ifdef DEBUG
	fprintf(stdout, "pivot row at step %5d is %5d\n", currow, pivotrow);
#endif
        for (i = currow; i < nsize; i++) {
            SWAP(matrix[pivotrow][i],matrix[currow][i]);
        }
        SWAP(R[pivotrow],R[currow]);
    }
}

extern char *optarg;

void errexit (const char *err_str)
{
    fprintf (stderr, "%s", err_str);
    exit (1);
}


void set_message(thread_arg *message, int nsize, int ii)
{
    message->nsize = nsize;
    message->task_id = ii + 1;
}

int set_begin(int nsize, int cnt_row, int task_id)
{
    int i = cnt_row;
    /* Factorize the rest of the matrix. */
    //rows that need to be dealt with, starting from (i + 1)_th. 
    //because i_th has been normalized. 
    int total_num_row = nsize-i-1; 
    
    //make sure that all row will be covered
    int section_num_row = (total_num_row - 1) / task_num + 1; 

    //task_id is 1 based index. 
    int begin = (i + 1) + (task_id - 1) * section_num_row;
    //printf("task_id %d, nsize %d, i %d, total_num_row %d, 
    //section_num_row %d, begin %d end %d\n",task_id, nsize, i,
    //total_num_row, section_num_row, begin, end);
    return begin;
}

int set_end(int nsize, int cnt_row, int task_id)
{
    int i = cnt_row;
    /* Factorize the rest of the matrix. */
    //rows that need to be dealt with, starting from (i + 1)_th. 
    //because i_th has been normalized. 
    int total_num_row = nsize-i-1; 
    
    //make sure that all row will be covered
    int section_num_row = (total_num_row - 1) / task_num + 1; 

    //task_id is 1 based index. 
    int begin = (i + 1) + (task_id - 1) * section_num_row;
    int end   = begin + section_num_row - 1; 
    if(end >= nsize)
        end = nsize - 1;
    //printf("task_id %d, nsize %d, i %d, total_num_row %d, 
    //section_num_row %d, begin %d end %d\n",task_id, nsize, i,
    //total_num_row, section_num_row, begin, end);
    return end;
}

/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */

void *computeGauss_row_version(void *arg)
{
    thread_arg *message = (thread_arg *) arg;
    int i, j, k, nsize, task_id, begin, end;
    double pivotval;
    pthread_attr_t attr;
    pthread_t *tid; 

    task_id = message->task_id;
    nsize = message->nsize;
    for (i = 0; i < nsize; i++) {
        if(task_id == 1){
            getPivot(nsize,i);

            /* Scale the main row. */
            pivotval = matrix[i][i];
            if (pivotval != 1.0) {
                matrix[i][i] = 1.0;
                for (j = i + 1; j < nsize; j++) {
                    matrix[i][j] /= pivotval; 
                }
                R[i] /= pivotval;
            }
            barrier(task_num);
        }
        else
            barrier(task_num);


       begin = set_begin(nsize, i, task_id);
       end = set_end(nsize, i, task_id);
       for (j = begin; j <= end; j++) {
            pivotval = matrix[j][i];
            matrix[j][i] = 0.0;
            for (k = i + 1; k < nsize; k++) {
                matrix[j][k] -= pivotval * matrix[i][k];
            }
            R[j] -= pivotval * R[i];
        }
        barrier(task_num);
    }
}

/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */
void *computeGauss_col_version(void *arg)
{
    thread_arg *message = (thread_arg *) arg;
    int i, j, k, nsize, task_id, begin, end;
    double pivotval;
    pthread_attr_t attr;
    pthread_t *tid; 

    task_id = message->task_id;
    nsize = message->nsize;
    for (i = 0; i < nsize; i++) {
        if(task_id == 1){
            getPivot(nsize,i);

            /* Scale the main row. */
            pivotval = matrix[i][i];
            if (pivotval != 1.0) {
                matrix[i][i] = 1.0;
                for (j = i + 1; j < nsize; j++) {
                    matrix[i][j] /= pivotval; 
                }
                R[i] /= pivotval;
            }
            barrier(task_num);
        }
        else
            barrier(task_num);

       begin = set_begin(nsize, i, task_id);
       end = set_end(nsize, i, task_id);
       for (j = i+1; j < nsize; j++) {
            pivotval = matrix[j][i];
            matrix[j][i] = 0.0;
            for (k = begin; k <= end; k++) {
                matrix[j][k] -= pivotval * matrix[i][k];
            }
            if(task_id == 1)
                R[j] -= pivotval * R[i];
        }

        barrier(task_num);
    }
}
/* Solve the equation. */

void solveGauss(int nsize)
{
    int i, j;

    X[nsize-1] = R[nsize-1];
    for (i = nsize - 2; i >= 0; i --) { 
        X[i] = R[i];
        for (j = nsize - 1; j > i; j--) {
            X[i] -= matrix[i][j] * X[j];
        }
    }

#ifdef DEBUG
    fprintf(stdout, "X = [");
    for (i = 0; i < nsize; i++) {
        fprintf(stdout, "%.6f ", X[i]);
    }
    fprintf(stdout, "];\n");
#endif
}



int main(int argc, char *argv[])
{
    int c, i, ii;
    int nsize = 0;
    double error;
    pthread_attr_t attr;
    thread_arg *messages;
    pthread_t *tid;
    char fname[1024];
    strcpy(fname,argv[1]);

    while ((c = getopt (argc, argv, "p:")) != -1)
        switch (c) {
        case 'p':
            task_num = atoi (optarg);
            break;
        }

    printf ("%d tasks\n", task_num);
    //printf ("%d\n",argc);

    nsize = initMatrix(fname);
    initRHS(nsize);
    initResult(nsize);

    //start timer
    gettimeofday(&start, 0);

    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    messages = (thread_arg *) malloc(sizeof(thread_arg)*task_num);
    tid = (pthread_t *) malloc(sizeof(pthread_t) * task_num);
    if(!messages || !tid)
        errexit("out of shared memory");

    for( ii = 0; ii < task_num; ii++){
        set_message(&messages[ii],nsize,ii);
        pthread_create(&tid[ii], &attr, computeGauss_row_version, &messages[ii]);
    }
    
    //wait for all thread to finish
    for( ii = 0; ii < task_num; ii++){
        pthread_join(tid[ii], NULL);
    }
    //end timer
    gettimeofday(&finish, 0);

    solveGauss(nsize);
    
    fprintf(stdout, "Time:  %f seconds\n", 
            (finish.tv_sec - start.tv_sec) + 
            (finish.tv_usec - start.tv_usec)*0.000001);

    error = 0.0;
    for (i = 0; i < nsize; i++) {
        double error__ = (X__[i]==0.0) ? 1.0 : fabs((X[i]-X__[i])/X__[i]);
        if (error < error__) {
            error = error__;
        }
    }
    fprintf(stdout, "Error: %e\n", error);

    return 0;
}
