#ifndef INPUT_H_H
#define INPUT_H_H

#include <stddef.h>

/*Input parameters*/
struct parameters {
    /*Total number of nodes in the graph*/
    size_t n;

    /*Total number of edges in graph. */
    size_t e;

    /* number of threads */
    size_t nthreads;

    /*For printing general statistics. 0=No print 1=print. Default=1*/
    size_t print_results;

    /*For printing the resulting adjecency matrix. Default=0*/
    size_t print_output;

    /*Default:0 . If set to 1, will run sequential version after the parallel version to compare results*/
    /*WARNING: it will slow down total execution time. Consider the 15 min time limit on DAS5*/
    int validate;

    /*Input*/
    int *input_matrix;

};


void read_parameters(struct parameters* p, int argc, char **argv);


#endif