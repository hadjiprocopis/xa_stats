/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: utils.c,v 1.8 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define EXTERN
#include "utils.h"

char
*mem_alloc(int n)
{   
  char *p;

  if ((p = (char *) malloc(n)) == NULL) 
    exit_on_error("Error on memory allocation");
  
  return p;
  
}

void
exit_on_error(char *str)
{
    
  fprintf(stderr,"%s\n", str);
  
  exit(1);
  
}

void
print_array_d(double data[], int n)
{
  int i;
  
  printf("[ ");
  for (i = 0; i < n; i++) printf("%10.6g  ", data[i]);
  printf(" ]\n");
  
}

void
print_array_i(int data[], int n)
{
  int i;
  
  printf("[ ");
  for (i = 0; i < n; i++) printf("%3i ", data[i]);
  printf(" ]\n");
  
}

void
print_mat_d(double **data, int n, int m)
{
  int i, j;
  double val;
  
  for (i = 0; i < n; i++) {
    printf("[ ");
    for (j = 0; j < m; j++) {
      if (fabs(data[i][j]) < EPS) {
	val = 0.0;
      } else {
	val = data[i][j];
      }
      printf("%10.6g  ", val);
    }
    printf(" ]\n");
  }
  
}

void
print_mat_i(int **data, int n, int m)
{
  int i, j;
  
  for (i = 0; i < n; i++) {
    printf("[ ");
    for (j = 0; j < m; j++) {
      printf("%3i ", data[i][j]);
    }
    printf(" ]\n");
  }
  
}

void
print_anova_data(double data[], int n, int **factor, int c)
{

  int i, j;

  printf ("   S  ");
  for (i = 0; i < c; i++) {
    printf (" F%-3i", i);
  }
  printf (" val\n");
  for (i = 0; i < n; i++) {    
    printf (" %3i ", i+1);
    for (j = 0; j < c; j++) {
      printf ("%4i ", factor[j][i]);
    }
    printf("  %-.3g\n", data[i]);
  }  

}

