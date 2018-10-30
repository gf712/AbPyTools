#ifndef ABPYTOOLS_OPS_H
#define ABPYTOOLS_OPS_H

double norm_op(const double *A, int p, int size);
void exp_op(const double *A, double *B, int size);
double l2_norm_op(const double *A, int size);
double l1_norm_op(const double *A, int size);

void add_op(const double *A, const double *B, double *C, int size);
void subtract_op(const double *A, const double *B, double *C, int size);
void multiply_op(const double *A, const double *B, double *C, int size);

#endif // ABPYTOOLS_OPS_H