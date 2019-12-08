#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include "apsp_kernel.h"
#include "input.h"
#include "output.h"


int main(int argc, char **argv){
    struct timeval total_before, total_after, before, after;
    gettimeofday(&total_before, NULL);
    struct parameters p;
    struct results r;
    int *dist, *val;


    //Read cmd line arguements and init everything
    read_parameters(&p, argc, argv);
    omp_set_num_threads(p.nthreads);

    int size=p.n;
    dist = p.input_matrix;
    r.size = size;

    /*If validate set to true, create a seperate copy of the input*/
    if (p.validate){
        printf("Validate set to true. Preparing 2nd copy of input for seq...\n");
        val = malloc(size * size * sizeof(int));
        for(int i=0; i<size; i++){
            for(int j=0; j<size; j++){
                int idx = i * size +j;
                val[idx] = dist[idx];
            }
        }
    }


    int (*restrict distances)[size][size] = (int (*)[size][size])dist;

    if(p.print_output){printf("Input matrix:\n"); print_distances(p.input_matrix, size); printf("\n");}//Is here for debugging!
    
    /*Call compute function*/
    gettimeofday(&before, NULL);
    //apsp_prototype(size, dist);
    //apsp_seq(size,distances);
    //apsp_simd(size,distances);
    //apsp_seq_bit(size,distances);
    //apsp_omp(size,distances);
    //apsp_omp(size,distances);
    apsp_pthread(size,distances,p.nthreads);
    gettimeofday(&after, NULL);

    double algorithm_time = ((after.tv_sec + (after.tv_usec / 1000000.0)) -
            (before.tv_sec + (before.tv_usec / 1000000.0)));

    if(p.print_output){printf("Output matrix:\n"); print_distances(dist, size);}

    gettimeofday(&total_after, NULL);
    double total_time = ((total_after.tv_sec + (total_after.tv_usec / 1000000.0)) -
            (total_before.tv_sec + (total_before.tv_usec / 1000000.0)));

    /*Only validate when seq time + par time < DAS5 time limit, or when preserve is used */
    if (p.validate){
        printf("Validate set to true. Running sequential version...\n");
        int (*restrict validate)[size][size] = (int(*)[size][size])val;
        apsp_seq(size, validate);
        int faults=0;
        for(int i=0; i<size; i++){
            for(int j=0; j<size; j++){
                if((*validate)[i][j] != (*distances)[i][j]){
                    //printf("Wrong value at i=%d , j=%d \n",i,j);
                    faults++;
                }
            }
        }
        if(faults==0){
            printf("Validation successful!\n");
            free(val);
        } else {
            printf("Validate error: Total number of inconsistencies: %d\n", faults);
            free(val);
        }
    }

    r.total_app_time = total_time;
    r.algorithm_time = algorithm_time;
    report_output(&r, size, algorithm_time);

    free(dist);
    return 0;
}
