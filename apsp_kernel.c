/*Contains the variations of the main computations functions of APSP*/
/*Floyd-Warshall algorithm implementations*/


//#define ENABLE_PREFETCH //Tested, but no difference
#define _GNU_SOURCE
#include <sched.h>
#include <sys/time.h>
#include "apsp_kernel.h"
#include <omp.h>
#include <pthread.h>
#include <stdio.h>
#include <limits.h>
#include <smmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>

#define INF 9999999
#define UNROLL_8 8
#define UNROLL_4 4


pthread_barrier_t bar;

struct thread_par {
    int n, start, end, id;
    void *dst;
};



/*Sequential prototype*/
/*Adapted from https://www.programmingalgorithms.com/algorithm/floyd%E2%80%93warshall-algorithm?lang=C*/
void apsp_prototype(int n, int *dist){
    int i,j,k;
    for(k=0; k<n; k++){
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                int ij=i*n+j;
                int ik=i*n+k;
                int kj=k*n+j;

                if(dist[ik]+dist[kj]< dist[ij]){
                    dist[ij] = dist[ik]+dist[kj];
                }

            }
        }
    }
}


/*Use this as the reference implementation*/
void apsp_seq(int n, int (*restrict dist)[n][n]){
    int i,j,k;

    for(k=0; k<n; k++){
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                int sum = (*dist)[i][k] +(*dist)[k][j];
                if(sum<(*dist)[i][j]){
                    (*dist)[i][j]=sum;
                }
            }
        }
    }
}




/*https://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax*/
/*Works, but apparently is slower than seq due to sophisticated hardware branch predictor*/
void apsp_seq_bit(int n, int (*restrict dist)[n][n]){
    int i,j,k;
    int int_size = sizeof(int);
    for(k=0; k<n; k++){
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                int cur_weight = (*dist)[i][j]; //as y
                int sum = (*dist)[i][k] +(*dist)[k][j]; //as x
                int r = cur_weight ^ ((sum ^ cur_weight) & -(sum < cur_weight));
                (*dist)[i][j] = r;

            }
        }
    }
}



/*Very bad Performance due to thread creation/destruction on each iteration of k.
 Commented code used for testing execution times for each iteration of k */
void apsp_omp(int n, int (*restrict dist)[n][n]){
    int i,j,k;
    omp_set_num_threads(16);

    for(k=0; k<n; k++){
        #pragma omp parallel for proc_bind(spread)
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                int sum = (*dist)[i][k] +(*dist)[k][j];
                (*dist)[i][j] = (sum<(*dist)[i][j]) ? sum : (*dist)[i][j];
            }
        }
    }
}




/*---------------------------------------------------------------------------------------------*/
/*                                       PTHREAD CODE                                          */

/*Created to test some differences of sse vs avx*/
void *calculate_paths_sse(struct thread_par *par){
    int i,j,k,l,ik;
    int N=par->n;
    int start_i = par->start;
    int end_i = par->end;
    int unroll_end = (int)(N/UNROLL_4) * UNROLL_4;
    int l_start = unroll_end; 
    int leftover = N % UNROLL_4; //Exact leftover
    int (*restrict dist)[N][N] = par->dst;
    const int core_id = par->id;
    const pthread_t pid = pthread_self();

    cpu_set_t cpuset; //Create empty bitset
    CPU_ZERO(&cpuset); //Init empty
    CPU_SET(core_id, &cpuset); //add cpu
    /*Pin the threads*/
    pthread_setaffinity_np(pid,sizeof(cpu_set_t), &cpuset);

    for(k = 0; k < N; ++k){
        for(i=start_i; i<=end_i; ++i){
            for(j = 0; j < unroll_end; j+=UNROLL_4){
                __m128i ij_values = _mm_load_si128(&(*dist)[i][j]);
                __m128i ik_values = _mm_set1_epi32((*dist)[i][k]);
                __m128i kj_values = _mm_load_si128(&(*dist)[k][j]);
                __m128i sums = _mm_add_epi32 (ik_values, kj_values);
                __m128i results = _mm_min_epi32(sums, ij_values);
                _mm_store_si128(&(*dist)[i][j], results);
            }
            for(l=l_start; l<N; ++l){//Deal with leftover if any
                int sum = (*dist)[i][k] +(*dist)[k][l];
                (*dist)[i][l] = (sum<(*dist)[i][l]) ? sum : (*dist)[i][l];             
            }

        }
        pthread_barrier_wait(&bar); //Make sure threads move to the next iteration of k together
    }
    return NULL;
}



/*Our FINAL version of our implementation. Results and graphs are based on this version. Unrolling by 8 to use AVX. Assumes N>8*/
void *calculate_paths_simd(struct thread_par *par){
    int i,j,k,l,ik,in,kn;
    int N=par->n;
    int start_i = par->start;
    int end_i = par->end;
    int unroll_end = (int)(N/UNROLL_8) * UNROLL_8;
    int l_start = unroll_end; 
    int leftover = N % UNROLL_8; //Exact leftover

    const int core_id = par->id;
    pthread_t pid = pthread_self();

    cpu_set_t cpuset; //Create empty bitset
    CPU_ZERO(&cpuset); //Init empty
    CPU_SET(core_id, &cpuset); //add id
    /*Pin the threads*/
    pthread_setaffinity_np(pid,sizeof(cpu_set_t), &cpuset);

    int *dist = par->dst;

    for(k = 0; k < N; ++k){
        kn=k*N;
        for(i=start_i; i<=end_i; ++i){
            ik = i*N+k;
            in= i*N;
            if(dist[ik]<INF){
                for(j = 0; j < unroll_end; j+=UNROLL_8){
                    int ij = in+j;       
                    int kj = kn+j;
                
                    __m256i ij_values = _mm256_load_si256(&dist[ij]);
                    __m256i ik_values = _mm256_set1_epi32(dist[ik]);
                    __m256i kj_values = _mm256_load_si256(&dist[kj]);

                    __m256i sums = _mm256_add_epi32 (ik_values, kj_values);
                    __m256i results = _mm256_min_epi32(sums, ij_values);

                    _mm256_storeu_si256(&dist[ij], results);
                }
                for(l=l_start; l<N; ++l){//Deal with leftover if any 
                    int kl=kn+l;
                    int il=in+l;
                    int sum = dist[ik] + dist[kl];
                    dist[il] = (sum<dist[il]) ? sum : dist[il];        
                }
            }
        }
        pthread_barrier_wait(&bar); //Make sure threads move to the next iteration of k together
    }
    return NULL;
}





/*The thread routine v1. We used this as a basis to create our pthread+simd version
 gcc reported vectorisation of this loop*/
void *calculate_paths(struct thread_par *par){
    int i,j,k;
    int N=par->n;
    int start_i = par->start;
    int end_i = par->end;

    const int core_id = par->id;
    const pthread_t pid = pthread_self();

    cpu_set_t cpuset; //Create empty bitset
    CPU_ZERO(&cpuset); //Init empty
    CPU_SET(core_id, &cpuset); //add cpu
    /*Pin the threads*/
    pthread_setaffinity_np(pid,sizeof(cpu_set_t), &cpuset);

    int (*restrict dist)[N][N] = par->dst;

    for(k = 0; k < N; k++){
        for(i=start_i; i<=end_i; i++){
            for(j = 0; j < N; j++){ //GCC reports vectorisation of this loop
                int sum = (*dist)[i][k] + (*dist)[k][j];
                (*dist)[i][j] = (sum<(*dist)[i][j]) ? sum : (*dist)[i][j];
            }
        }
        pthread_barrier_wait(&bar); //Make sure threads move to the next iteration of k together
    }

    return NULL;
}





/*Thread initialisation function*/
void apsp_pthread(int n, int (*restrict dist)[n][n], int thread_num){
    if(n<16)thread_num=2; //Index calculation breaks if thread_num > n

    pthread_t threads[thread_num];
    pthread_barrier_init(&bar, NULL, thread_num);
    struct thread_par thread_p[thread_num];
    /*Assign the thread parameters*/
    int ci=0;
    int rows_per_thread = (int) n/thread_num;
    int remaining_rows = n % thread_num;
    /*Calculates index ranges for threads. If n is not perfectly divisible by thread num then
      the extra rows are assigned on a 1 extra row per thread basis, starting from 1st thread.
      Most fair distribution we can use, instead of dumping all the extra work onto a single
      thread(load imbalance)*/
    for (int i=0; i<thread_num; i++){
        if (remaining_rows==0){
            thread_p[i].start = i * rows_per_thread;
            thread_p[i].end = thread_p[i].start + rows_per_thread-1;
        } else {
            if(i<=remaining_rows){
                ci = i;
            }
            thread_p[i].start = (i * rows_per_thread)+ci;
            thread_p[i].end = (i < remaining_rows) ? thread_p[i].start + rows_per_thread : thread_p[i].start + rows_per_thread-1;
        }
        thread_p[i].n = n;
        thread_p[i].dst = (void*)dist;  
        thread_p[i].id = i;
    }



    for (int t=0; t<thread_num; t++){
        pthread_create(&threads[t], NULL, (void*)&calculate_paths , &thread_p[t]);
        //pthread_create(&threads[t], NULL, (void*)&calculate_paths_simd , &thread_p[t]);
        //pthread_create(&threads[t], NULL, (void*)&calculate_paths_sse , &thread_p[t]);
    }

    for (int t=0; t<thread_num; t++){
        pthread_join(threads[t],NULL);
    }

}
