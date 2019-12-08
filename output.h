#ifndef OUTPUT_H
#define OUTPUT_H

#include "input.h"
#include <stddef.h>

struct results{
    double total_app_time;
    double algorithm_time;
    double ops;
    double size;
};


void print_distances(int *distances, int N);

void report_output(struct results *r, int size, double time);

#endif