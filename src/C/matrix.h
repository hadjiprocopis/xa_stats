/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: matrix.h,v 1.7 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef matrix_h
#define matrix_h

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
EXTERN int verbose;
EXTERN int debug;
#include "utils.h"
#include "files.h"
#include "solution.h"

MATRIX *
mem_allocate_matrix(int n, int m);

MATRIX *
transpose(MATRIX *m,
	  int w);

double
dot(MATRIX *m1,
    MATRIX *m2,
    int mm1,
    int mm2);

MATRIX *
vector_to_matrix(VECTOR *v,
		 int w);

double
determinant(MATRIX *m,
	    int w);

MATRIX *
cofactor(MATRIX *m,
	 int t,
	 int w);

MATRIX * 
multiply1(MATRIX *m1,
	  MATRIX *m2,
	  int w);

MATRIX * 
multiply(MATRIX *m1,
	 MATRIX *m2,
	 int w);

MATRIX * 
add(MATRIX *m1,
    MATRIX *m2,
    int w);

MATRIX * 
subtract(MATRIX *m1,
	 MATRIX *m2,
	 int w);

MATRIX *
inv(MATRIX *m,
    int w);

MATRIX *
diag(VECTOR *v);

MATRIX *
fill(int m,
     int n,
     double f);

MATRIX *
identity(int n);

MATRIX *
remove_cols(MATRIX *m,
	    int *filter,
	    int col);

MATRIX *
augment_cols(MATRIX *m1,
	     MATRIX *m2,
	     int w);

MATRIX *
augment_rows(MATRIX *m1,
	     MATRIX *m2);

MATRIX *
duplicate(MATRIX *m);

MATRIX *
copy(MATRIX *m,
     int n1,
     int m1,
     int n2,
     int m2);

MATRIX *
kronecker(MATRIX *m1,
	  MATRIX *m2);

MATRIX *
cross(MATRIX *m1,
      MATRIX *m2,
      int w);

MATRIX *
invert(MATRIX *m,
       int w);

MATRIX *
pinvert(MATRIX *mat, VECTOR *vect);

MATRIX *
g2invert(MATRIX *m,
	 int k,
	 int **rank,
	 int *r,
	 int w);

void
sweep(MATRIX *m,
      int k1,
      int k2,
      int **rank,
      int *r);

void
mat_rank(MATRIX *m,
	 int **rank,
	 int *r);

MATRIX *
independent(MATRIX *m);

void
dump_matrix(MATRIX *mat);

#endif
