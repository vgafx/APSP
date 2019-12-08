#ifndef APSP_KERNEL_H
#define APSP_KERNEL_H


void apsp_prototype(int n, int *dist);

void apsp_seq(int n, int (*restrict dist)[n][n]);

void apsp_simd(int n, int (*restrict dist)[n][n]);

void apsp_seq_bit(int n, int (*restrict dist)[n][n]);

void apsp_omp(int n, int (*restrict dist)[n][n]);

void apsp_pthread(int n, int (*restrict dist)[n][n], int thread_num);

#endif