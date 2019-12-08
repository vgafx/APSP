#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include "output.h"

/*Count ops for Macro. We use ints for apsp*/
#define OPS_PER_ITERATION 2

#define INF 9999999


//!ops might overflow if size is too large.
inline  int64_t compute_ops(int size, double time){
    if(time<1.0)time=1.0; //avoid division by 0 for fast executions
    int64_t time64 = (int64_t) time; //Expected loss of percision here
    int64_t size_cubed = (int64_t) size;
    int64_t orig_size = (int64_t) size;

    size_cubed *=size_cubed;
    size_cubed *=orig_size;
    size_cubed *=OPS_PER_ITERATION;
    size_cubed /=time64;

    return size_cubed;
}



void print_distances(int *distances, int N){
    int i,j;

    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            (distances[i*N+j]==INF)? printf("INF \t") : printf("%d \t", distances[i*N+j]);
        }
        printf("\n");
    }
}



void report_output(struct results *r, int size, double time){
    int64_t new_ops = compute_ops(size,time);
    printf("Output: \n");
    printf("Total time:     %.6f \t\n",r->total_app_time);
    printf("Algorithm time: %.6f \t\n",r->algorithm_time);
    printf("OPS/s:          %" PRIu64 "\n", new_ops);
}