/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: matrix.c,v 1.7 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#define EXTERN
#include "matrix.h"

MATRIX *
mem_allocate_matrix(int n,
		    int m)
{

  int i;
  MATRIX *matrix;

  matrix = TYPE_ALLOC(MATRIX);
  matrix->n = n;
  matrix->m = m;
  matrix->data = ARRAY_ALLOC(n, double *);
  for (i=0;i<n;i++) {
    matrix->data[i] = ARRAY_ALLOC(m, double);
  }

  return matrix;

}

MATRIX *
transpose(MATRIX *m,
	  int w)
{

  int i, j;
  MATRIX *mt;

  mt = mem_allocate_matrix(m->m, m->n);

  for (i=0;i<m->n;i++) {
    for (j=0;j<m->m;j++) {
      mt->data[j][i]=m->data[i][j];
    }
  }

  if (w==1) {
    free(m);
  }

  return mt;

}

double
dot(MATRIX *m1,
    MATRIX *m2,
    int mm1,
    int mm2)
{

  if (m1->n != m2->n) {
    exit_on_error("dot:\nnot the same number of columns.");
  }

  int i;
  double d=0.0;

  for (i=0;i<m1->n;i++) {
    d = d + (m1->data[i][mm1] * m2->data[i][mm2]);
  }
  
  return d;

}

MATRIX *
vector_to_matrix(VECTOR *v,
		 int w)
{

  int i;
  MATRIX *m;

  m = mem_allocate_matrix(v->n, 1);

  for (i=0;i<v->n;i++) {
    m->data[i][0]=v->data[i];
  }

  if (w==1) {
    free(v);
  }

  return m;

}

double
determinant(MATRIX *m,
	    int w)
{

  int i, j, k, l;
  double det = 0;
  MATRIX *m1;

  if (m->n < 1) {
    exit_on_error("error in determinant");
  } else if (m->n == 1) {
    det = m->data[0][0];
  } else if (m->n == 2) {
    det = m->data[0][0] * m->data[1][1] - m->data[1][0] * m->data[0][1];
  } else {
    det = 0;
    for (i=0;i<m->n;i++) {
      m1 = mem_allocate_matrix(m->n-1, m->n-1);
      for (j=1;j<m->n;j++) {
	k = 0;
	for (l=0;l<m->n;l++) {
	  if (l == i) {
	    continue;
	  }
	  m1->data[j-1][k] = m->data[j][l];
	  k++;
	}
      }
      det += pow(-1.0,i+2.0) * m->data[0][i] * determinant(m1,1);
    }
  }

  if (w==1) {
    free(m);
  }

  return det;

}

MATRIX *
cofactor(MATRIX *m,
	 int t,
	 int w)
{
   int i,j,ii,jj,i1,j1;
   double det;
   MATRIX *m1, *c; 

   m1 = mem_allocate_matrix(m->n-1, m->n-1);
   c = mem_allocate_matrix(m->n, m->m);

   for (j=0;j<m->n;j++) {
     for (i=0;i<m->n;i++) {
       // Form the adjoint a_ij
       i1 = 0;
       for (ii=0;ii<m->n;ii++) {
	 if (ii == i) {
	   continue;
	 }
	 j1 = 0;
	 for (jj=0;jj<m->n;jj++) {
	   if (jj == j) {
	     continue;
	   }
	   m1->data[i1][j1] = m->data[ii][jj];
	   j1++;
	 }
	 i1++;
       }
       // Calculate the determinate
       if (m->n == 1) {
	 det = 1;
       } else {
	 det = determinant(m1,0);
       }
       // Fill in the elements of the
       // cofactor
       if (t == 0) {
	 c->data[i][j] = pow(-1.0,i+j+2.0) * det;
       } else {
	 // if we transpose it
	 c->data[j][i] = pow(-1.0,i+j+2.0) * det;
       }
     }
   }
   free(m1);
   if (w==1) {
     free(m);
   }

   return c;

}

MATRIX * 
multiply1(MATRIX *m1,
	  MATRIX *m2,
	  int w)
{

  int i, j, k;
  double val;
  MATRIX *m;

  // Simple method

  if (m1->m != m2->n) {
    exit_on_error("multiply:\nmatrices not compatible");
  }  

  m = mem_allocate_matrix(m1->n, m2->m);

  for(i=0;i<m1->n;i++) {
    for(j=0;j<m2->m;j++) {
      val=0.0;
      for(k=0;k<m2->n;k++) {
	val += m1->data[i][k]*m2->data[k][j];
      }
      m->data[i][j]=val;
    }
  }

  if (w==1) {
    free(m1);
    free(m2);
  }

  return m;

}

MATRIX * 
multiply(MATRIX *m1,
	 MATRIX *m2,
	 int w)
{

  int i, j, k;
  double row[m1->n], col[m2->m];
  MATRIX *m;

  // Winograd Method

  if (m1->m != m2->n) {
    exit_on_error("multiply:\nmatrices not compatible");
  }  

  m = mem_allocate_matrix(m1->n, m2->m);

  if (m1->m%2) {
    for(i=0;i<m1->n;i++) {
      row[i]=0;
      for(j=1;j<m1->m;j+=2) {
	row[i] += m1->data[i][j] * m1->data[i][j-1]; 
      }
    }
    for(i=0;i<m2->m;i++) {
      col[i]=0;
      for(j=1;j<m2->n;j+=2) {
	col[i] += m2->data[j][i] * m2->data[j-1][i]; 
      }
    }
    for(i=0;i<m1->n;i++) {
      for(j=0;j<m2->m;j++) {
	m->data[i][j] = m1->data[i][m1->m-1] * m2->data[m2->n-1][j] - row[i] - col[j]; 
	for(k=1;k<m2->n;k+=2) {
	  m->data[i][j] += (m1->data[i][k-1]+m2->data[k][j]) * (m1->data[i][k]+m2->data[k-1][j]);
	}
      }
    }
  } else {
    for(i=0;i<m1->n;i++) {
      row[i]=0;
      for(j=1;j<m1->m;j+=2) {
	row[i] += m1->data[i][j] * m1->data[i][j-1]; 
      }
    }
    for(i=0;i<m2->m;i++) {
      col[i]=0;
      for(j=1;j<m2->n;j+=2) {
	col[i] += m2->data[j][i] * m2->data[j-1][i]; 
      }
    }
    for(i=0;i<m1->n;i++) {
      for(j=0;j<m2->m;j++) {
	m->data[i][j] = - row[i] - col[j]; 
	for(k=1;k<m2->n;k+=2) {
	  m->data[i][j] += (m1->data[i][k-1]+m2->data[k][j]) * (m1->data[i][k]+m2->data[k-1][j]);
	}
      }
    }
  }

  if (w==1) {
    free(m1);
    free(m2);
  }

  return m;

}

MATRIX * 
add(MATRIX *m1,
    MATRIX *m2,
    int w)
{

  int i, j;
  MATRIX *m;

  if (m1->n != m2->n || m1->m != m2->m) {
    exit_on_error("add:\nmatrices not of the same dimenssions");
  }  

  m = mem_allocate_matrix(m1->n, m2->m);

  for(i=0;i<m1->n;i++) {
    for(j=0;j<m1->m;j++) {
      m->data[i][j] = m1->data[i][j]+m2->data[i][j];
    }
  }

  if (w==1) {
    free(m1);
    free(m2);
  }

  return m;

}

MATRIX * 
subtract(MATRIX *m1,
	 MATRIX *m2,
	 int w)
{
  
  int i, j;
  MATRIX *m;
  
  if (m1->n != m2->n || m1->m != m2->m) {
    exit_on_error("subtract:\nmatrices not of the same dimenssions");
  }  

  m = mem_allocate_matrix(m1->n, m2->m);

  for(i=0;i<m1->n;i++) {
    for(j=0;j<m1->m;j++) {
      m->data[i][j] = m1->data[i][j]-m2->data[i][j];
    }
  }

  if (w==1) {
    free(m1);
    free(m2);
  }

  return m;

}

MATRIX *
inv(MATRIX *m,
    int w)
{

  int i, j;
  double det;
  MATRIX *im;

  det = determinant(m, 0);

  im = cofactor(m, 0, 1);

  for(i=0;i<im->n;i++) {
    for(j=0;j<im->m;j++) {
      im->data[i][j] /= det;
    }
  }

  if (w==1) {
    free(m);
  }

  return im;

}

MATRIX *
diag(VECTOR *v)
{

  // I haven't tested this routine

  int i;
  MATRIX *m;
  MATRIX * fill(int n, int m, double d);

  m = fill(v->n, v->n, (double)0.0);
  
  for(i=0;i<m->n;i++) {
    m->data[i][i]=v->data[i];
  }

  free(v);

  return m;

}

MATRIX *
fill(int n,
     int m,
     double f)
{

  int i, j;
  MATRIX *mat;

  mat = mem_allocate_matrix(n, m);

  for(i=0;i<mat->n;i++) {
    for(j=0;j<mat->m;j++) {
      mat->data[i][j]=f;
    }
  }

  return mat;

}

MATRIX *
identity(int n)
{

  // I haven't tested this routine

  int i;
  MATRIX *m;
  MATRIX * fill(int n, int m, double d);

  m = fill(n, n, 0.0);
  
  for(i=0;i<m->n;i++) {
    m->data[i][i]=1.0;
  }

  return m;

}

MATRIX *
remove_cols(MATRIX *m,
	    int *filter,
	    int col)
{

  int i, j ,jj;
  MATRIX *mm;

  mm = mem_allocate_matrix(m->n, col);

  for(i=0;i<m->n;i++) {
    jj=0;
    for(j=0;j<m->m;j++) {
      if (filter[j]>0) {
	mm->data[i][jj] = m->data[i][j];
	jj++;
      }
    }
  }

  free(m);

  return mm;

}

MATRIX *
augment_cols(MATRIX *m1,
	     MATRIX *m2,
	     int w)
{

  int i, j, jj;
  MATRIX *m;

  if (m1->n != m2->n) {
    exit_on_error("augment_cols:\nnot the same number of rows in matrices");
  }

  m = mem_allocate_matrix(m1->n, m1->m+m2->m);

  for(i=0;i<m1->n;i++) {
    for(j=0;j<m1->m;j++) {
      m->data[i][j]=m1->data[i][j];
    }
  }
  for(i=0;i<m2->n;i++) {
    jj=m1->m;
    for(j=0;j<m2->m;j++) {
      m->data[i][jj]=m2->data[i][j];
      jj++;
    }
  }

  if (w==1) {
    free(m1);
    free(m2);
  }

  return m;

}

MATRIX *
augment_rows(MATRIX *m1,
	     MATRIX *m2)
{

  // I haven't tested this routine

  int i, j, ii;
  MATRIX *m;

  if (m1->m != m2->m) {
    exit_on_error("augment_rows:\nnot the same number of cols in matrices");
  }

  m = mem_allocate_matrix(m1->n+m2->n, m1->m+m2->m);

  for(i=0;i<m1->n;i++) {
    for(j=0;j<m1->m;j++) {
      m->data[i][j]=m1->data[i][j];
    }
  }
  ii=m1->n;
  for(i=0;i<m2->n;i++) {
    for(j=0;j<m2->m;j++) {
      m->data[ii][j]=m2->data[i][j];
    }
    ii++;
  }

  free(m1);
  free(m2);

  return m;

}

MATRIX *
duplicate(MATRIX *m)
{

  int i, j;
  MATRIX *dm;

  dm = mem_allocate_matrix(m->n, m->m);

  for(i=0;i<dm->n;i++) {
    for(j=0;j<dm->m;j++) {
      dm->data[i][j]=m->data[i][j];
    }
  }

  return dm;

}

MATRIX *
copy(MATRIX *m,
     int n1,
     int m1,
     int n2,
     int m2)
{

  int i, j, ii, jj;
  MATRIX *cp;

  cp = mem_allocate_matrix(n2-n1,m2-m1);

  ii=0;
  for(i=n1;i<n2;i++) {
    jj=0;
    for(j=m1;j<m2;j++) {
      cp->data[ii][jj]=m->data[i][j];
      jj++;
    }
    ii++;
  }

  return cp;

}

MATRIX *
kronecker(MATRIX *m1,
	  MATRIX *m2)
{

  // I haven't tested this routine

  int i, j, k, l, ii, jj;
  MATRIX *m;

  m = mem_allocate_matrix(m1->n*m2->n, m1->m*m2->m);

  for(i=ii=0;i<m1->n;i++) {
    for(j=jj=0;j<m1->m;j++) {
      for(k=0;k<m2->n;k++) {
	for(l=0;l<m2->m;l++) {
	  m->data[ii+k][jj+l] = m1->data[i][j]*m2->data[k][l];
	}
      }
    }
  }

  free(m1);
  free(m2);

  return m;

}

MATRIX *
cross(MATRIX *m1,
      MATRIX *m2,
      int w)
{

  int i, j, k, jj;
  int ncols=0, cols=m1->m*m2->m;
  double zero[cols];
  MATRIX *mt, *m;

  if (m1->n != m2->n) {
    exit_on_error("cross:\nnot the same number of rows in matrices");
  }

  mt = mem_allocate_matrix(m1->n, cols);

  for(i=0;i<cols;i++) {
    zero[i]=0;
  }

  for(i=0;i<m1->n;i++) {
    jj=0;
    for(j=0;j<m1->m;j++) {
      for(k=0;k<m2->m;k++) {
	mt->data[i][jj] = m1->data[i][j]*m2->data[i][k];
	zero[jj] += abs(mt->data[i][jj]);
	jj++;
      }
    }
  }

  if (w==1) {
    free(m1);
    free(m2);
  }

  // Skip the null columns
  for(i=0;i<cols;i++) {
    if (zero[i]>0) {
      ncols++;
    }
  }

  m = mem_allocate_matrix(mt->n, ncols);

  for(i=0;i<mt->n;i++) {
    jj=0;
    for(j=0;j<cols;j++) {
      if (zero[j]>0) {
	m->data[i][jj] = mt->data[i][j];
	jj++;
      }
    }
  }

  free(mt);

  return m;

}

MATRIX *
invert(MATRIX *m,
       int w)
{

  int i, j, *indx;
  double **mat, d, *col;
  MATRIX *im;

  // Matrix inversion by lu decomposition

  void ludcmp(double **mat, int n, int *indx, double *d);
  void lubksb(double **mat, int n, int *indx, double *col);

  im = mem_allocate_matrix(m->n, m->m);

  mat = ARRAY_ALLOC(m->n, double *);
  for(i=0;i<m->n;i++) {
    mat[i] = ARRAY_ALLOC(m->m, double);
    for(j=0;j<m->m;j++) {
      mat[i][j] = m->data[i][j];
    }
  }

  col = ARRAY_ALLOC(m->n, double);
  indx = ARRAY_ALLOC(m->n, int);
  ludcmp(mat,m->n,indx,&d);
  for(i=0;i<m->n;i++) {
    for(j=0;j<m->m;j++) {
      col[j]=0.0;
    }
    col[i]=1.0;
    lubksb(mat,m->n,indx,col);
    for(j=0;j<m->m;j++) {
      im->data[i][j] = col[j];
    }
  }
  free(mat);
  free(col);
  free(indx);
  if (w==1) {
    free(m);
  }
  
  return im;

}

MATRIX *
pinvert(MATRIX *matrix,
	VECTOR *vector)
{

  int i, j;
  double **mat, *y, *w, **v, *x, min, max=0.0;
  MATRIX *im;

  // I haven't tested this routine

  void svdcmp(double **a, int m, int n, double *w, double **v);
  void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x);

  im = mem_allocate_matrix(matrix->n, matrix->m);

  mat = ARRAY_ALLOC(matrix->n, double *);
  y = ARRAY_ALLOC(matrix->n, double);
  w = ARRAY_ALLOC(matrix->n, double);
  v = ARRAY_ALLOC(matrix->n, double *);
  x = ARRAY_ALLOC(matrix->m, double);
  for(i=0;i<matrix->n;i++) {
    mat[i] = ARRAY_ALLOC(matrix->m, double);
    v[i] = ARRAY_ALLOC(matrix->n, double);
    for(j=0;j<matrix->m;j++) {
      mat[i][j] = matrix->data[i][j];
    }
    y[i] = vector->data[i];
  }

  svdcmp(mat,matrix->m,matrix->n,w,v);

  for(i=0;i<matrix->n;i++) {
    if (w[i] > max) {
      max = w[i];
    }
  }
  min = max*(1.0e-7);
  for(i=0;i<matrix->n;i++) {
    if (w[i] < min) {
      w[i] = 0.0;
    }
  }

  svbksb(mat,w,v,matrix->m,matrix->n,y,x);

  for(i=0;i<matrix->n;i++) {
    for(j=0;j<matrix->m;j++) {
      im->data[i][j] = mat[i][j]*w[i];
    }
  }

  free(mat);
  free(y);
  free(w);
  free(v);
  free(x);

  return im;

}

MATRIX *
g2invert(MATRIX *m,
	 int k,
	 int **rank,
	 int *r,
	 int w)
{

  MATRIX *im;

  void sweep(MATRIX *m, int k1, int k2, int **rank, int *r);

  im = duplicate(m);
  sweep(im,0,k,rank,r);

  if (w==1) {
    free(m);
  }

  return im;

}

void
sweep(MATRIX *m,
      int k1,
      int k2,
      int **rank,
      int *r)
{

  int i, j, k;
  double d;

  if (m->n != m->m) {
    exit_on_error("sweep:\nnot a square matrix");
  }
  
  if (k2 < k1) {
    k=k1;
    k1=k2;
    k2=k;
  }

  *r=0;
  for(i=0;i<m->m;i++) {
    (*rank)[i]=1;
  }

  for(k=k1;k<k2;k++) {
    if (fabs(m->data[k][k])<EPS) {
      (*rank)[k]=0;
      *r += 1;
      for(i=0;i<k2;i++) {
	m->data[i][k]=0;
	m->data[k][i]=0;
      }
    } else {
      d=1/m->data[k][k];
      m->data[k][k]=d;
      for(i=0;i<m->n;i++) {
	if (i != k) {
	  m->data[i][k] *= -d; 
	  m->data[k][i] *= d; 
	}
      }
      for(i=0;i<m->n;i++) {
	if (i != k) {
	  for(j=0;j<m->n;j++) {
	    if (j != k) {
	      m->data[i][j] += m->data[i][k]*m->data[k][j]/d;
	    }
	  }
	}
      }
    }
  }
  *r = m->m - *r;

}

void
mat_rank(MATRIX *m,
	 int **rank,
	 int *r)
{
  
  MATRIX *s;

  void sweep(MATRIX *m, int k1, int k2, int **rank, int *r);

  s = transpose(m,0);
  s = multiply(s,m,0);
  
  sweep(s,0,s->m,rank,r);

  free(s);

}

MATRIX *
independent(MATRIX *m)
{

  int *rank, r;
  MATRIX *im;

  void mat_rank(MATRIX *m, int **rank, int *r); 

  rank = ARRAY_ALLOC(m->m, int);

  mat_rank(m,&rank,&r);

  im = remove_cols(m,rank,r);

  return im;

}

void
dump_matrix(MATRIX *mat)
{

  print_mat_d(mat->data, mat->n, mat->m);

}

