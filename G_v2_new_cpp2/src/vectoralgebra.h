/*
 * vectoralgebra.h
 *
 *  Created on: Dec 4, 2018
 *      Author: alessandrobaroni
 */

#ifndef VECTORALGEBRA_H_
#define VECTORALGEBRA_H_

double dot_3d_product(double*a, double*b);
double dot_product(double *a,double* b,int length);
void matrix_prod_vec( int length,double (*m)[4],double* p, double* final);
void scal_prod_vec(int size,double a,double*init,double* final);

#endif /* VECTORALGEBRA_H_ */
