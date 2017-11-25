/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net 

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: metrics.c,v 1.9 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <math.h>
#include <stdlib.h>
#define EXTERN
#include "metrics.h"

double
euclidean(VECTOR *v1,
	  VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double d = *v1p++ - *v2p++;
    distance += d * d;
  }
  return sqrt(distance / v1->n);

}

double
ieuclidean(VECTOR *v1,
	   VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double d = (*v1p++ * -1) - *v2p++;
    distance += d * d;
  }
  return sqrt(distance / v1->n);

}

double
minkowski(VECTOR *v1,
	  VECTOR *v2,
	  double power)
{
  double distance = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double d = *v1p++ - *v2p++;
    distance += pow(fabs(d), power);
  }
  return (pow(distance, 1/power)) / v1->n;

}
 
double
iminkowski(VECTOR *v1,
	   VECTOR *v2,
	   double power)
{
  double distance = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double d = (*v1p++ * -1) - *v2p++;
    distance += pow(fabs(d), power);
  }
  return (pow(distance, 1/power)) / v1->n;

}

double
manhattan(VECTOR *v1,
	  VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    distance += fabs(*v1p++ - *v2p++);
  }
  return  distance / v1->n;

}
 
double
imanhattan(VECTOR *v1,
	   VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    distance += fabs((*v1p++ * -1) - *v2p++);
  }
  return  distance / v1->n;

}
 
double
maximum(VECTOR *v1,
	VECTOR *v2)
{
  double distance = 0.0;
  double value   = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    value = fabs(*v1p++ - *v2p++);
    if (value > distance) {
      distance = value;  
    }
  }
  return  distance;

}
 
double
imaximum(VECTOR *v1,
	 VECTOR *v2)
{
  double distance = 0.0;
  double value   = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    value = fabs((*v1p++ * -1) - *v2p++);
    if (value > distance) {
      distance = value;  
    }
  }
  return  distance;

}
 
double
pearson(VECTOR *v1,
	VECTOR *v2)
{

  double sumX    = 0.0;
  double sumY    = 0.0;
  double sumX2   = 0.0;
  double sumY2   = 0.0;
  double sumXY   = 0.0;
  double SXX     = 0.0;
  double SYY     = 0.0;
  double SXY     = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    sumX    += *v1p;
    sumY    += *v2p;
    sumX2   += *v1p * *v1p;
    sumY2   += *v2p * *v2p;
    sumXY   += *v1p++ * *v2p++;
  }
  SXX  = sumX2 - ((sumX * sumX) / v1->n);
  SYY  = sumY2 - ((sumY * sumY) / v2->n);
  SXY  = sumXY - (sumX * sumY) / v1->n;
  
  if (SXX * SYY == 0) {
    return 0.0;
  } else {
    return (SXY / (sqrt (SXX * SYY)));
  }

}
 
double
ipearson(VECTOR *v1,
	 VECTOR *v2)
{

  double sumX    = 0.0;
  double sumY    = 0.0;
  double sumX2   = 0.0;
  double sumY2   = 0.0;
  double sumXY   = 0.0;
  double SXX     = 0.0;
  double SYY     = 0.0;
  double SXY     = 0.0;
  int i;
  double *v1p = v1->data;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double invp = (*v1p++ * -1);
    sumX           += invp;
    sumY           += *v2p;
    sumX2          += invp * invp;
    sumY2          += *v2p * *v2p;
    sumXY          += invp * *v2p++;
  }
  SXX  = sumX2 - ((sumX * sumX) / v1->n);
  SYY  = sumY2 - ((sumY * sumY) / v2->n);
  SXY  = sumXY - (sumX * sumY) / v1->n;
  
  if (SXX * SYY == 0) {
    return 0.0;
  } else {
    return (SXY / (sqrt (SXX * SYY)));
  }

}

double
spearman(VECTOR *v1,
	 VECTOR *v2)
{

  int i = v1->n;
  int l = 0;
  double *v1p = v1->data;
  double *v2p = v2->data;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  spearman_corr_nna(v1p, v2p, i, l, &d, &zd, &probd, &rs, &probrs);

  return rs;

}
 
double
ispearman(VECTOR *v1,
	  VECTOR *v2)
{

  int i = v1->n;
  int l = 0;
  double *v1p = v1->data;
  double *v2p = v2->data;
  double *v3p;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  v3p = ARRAY_ALLOC(v1->n, double);

  for (i = 0; i < v1->n; i++) {
    v3p[i] = *v1p++ * -1;
  }
  
  spearman_corr_nna(v3p, v2p, i, l, &d, &zd, &probd, &rs, &probrs);
  
  return rs;

}

double
lspearman(VECTOR *v1,
	  VECTOR *v2)
{

  int i = v1->n;
  int l = 1;
  double *v1p = v1->data;
  double *v2p = v2->data;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  spearman_corr_nna(v1p, v2p, i, l, &d, &zd, &probd, &rs, &probrs);

  return rs;

}
 
double
ilspearman(VECTOR *v1,
	   VECTOR *v2)
{

  int i = v1->n;
  int l = 1;
  double *v1p = v1->data;
  double *v2p = v2->data;
  double *v3p;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  v3p = ARRAY_ALLOC(v1->n, double);

  for (i = 0; i < v1->n; i++) {
    v3p[i] = *v1p++ * -1;
  }
  
  spearman_corr_nna(v3p, v2p, i, l, &d, &zd, &probd, &rs, &probrs);
  
  return rs;

}
