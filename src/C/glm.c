/*********************************************************************
 Copyright 2004 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

 $Id: glm.c,v 1.25 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define EXTERN
#include "glm.h"

void
parse_arguments(int argc,
		char **argv,
		char **data_file,
		char **factor_file,
		char **covariate_file,
		char **output_file,
		int *samples,
		int *variables,
		char **variable_name,
		char ***fixed_effs,
		int *n_fixed_effs,
		char ***interacts,
		int *n_interacts,
		char ***random_facs,
		int *n_random_facs,
		int *ssq,
		double *Q,
		int *permuts)
{
  
  int c, arg_len, max_objs;
  extern char *optarg;
  int errflg = 0;
  
  *data_file = NULL;
  
  *factor_file = NULL;
  
  *covariate_file = NULL;

  *output_file = NULL;

  *variable_name = NULL;
  
  *samples = 0;
  
  *variables = 0;
  
  *n_fixed_effs = 0;

  *n_interacts = 0;

  *n_random_facs = 0;

  *ssq = 3;

  max_objs = 10;

  *fixed_effs = ARRAY_ALLOC(max_objs, char *);
  
  *interacts = ARRAY_ALLOC(max_objs, char *);

  *random_facs = ARRAY_ALLOC(max_objs, char *);
  
  *Q = 0.05;

  *permuts = 10;
  
  while ((c = getopt(argc, argv, "d:f:c:o:s:g:n:m:i:r:t:p:q:vD")) != EOF)
    switch (c) {
    case 'd':
      *data_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*data_file, optarg);
      break;
    case 'f':
      *factor_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*factor_file, optarg);
      break;
    case 'c':
      *covariate_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*covariate_file, optarg);
      break;
    case 'o':
      *output_file = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*output_file, optarg);
      break;
    case 's':
      *samples = atoi(optarg);
      break;
    case 'g':
      *variables = atoi(optarg);
      break;
    case 'n':
      *variable_name = ARRAY_ALLOC((strlen(optarg) + 1), char);
      strcpy(*variable_name, optarg);
      break;
    case 'm':
      if (*n_fixed_effs <= max_objs) {
	arg_len = strlen(optarg);
	(*fixed_effs)[*n_fixed_effs] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*fixed_effs)[*n_fixed_effs], optarg);
	*n_fixed_effs += 1;
      }
      break;
    case 'i':
      if (*n_interacts <= max_objs) {
	arg_len = strlen(optarg);
	(*interacts)[*n_interacts] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*interacts)[*n_interacts], optarg);
	*n_interacts += 1;
      }
      break;
    case 'r':
      if (*n_random_facs <= max_objs) {
	arg_len = strlen(optarg);
	(*random_facs)[*n_random_facs] = ARRAY_ALLOC((arg_len+1), char);
	strcpy((*random_facs)[*n_random_facs], optarg);
	*n_random_facs += 1;
      }
      break;
    case 't':
      *ssq = atoi(optarg);
      if (*ssq < 2 || *ssq > 3) {
	errflg++;
      }
      break;
    case 'q':
      *Q = atof(optarg);
      break;
    case 'p':
      *permuts = atoi(optarg);
      break;
    case 'h':
      errflg++;
      break;
    case 'v':
      verbose=1;
      break;
    case 'D':
      debug=1;
      break;
    case '?':
      errflg++;
      break;
    }

  if (*data_file == NULL || (*factor_file == NULL && *covariate_file == NULL)) { 
    errflg++;
  }
  if (*n_fixed_effs+*n_interacts==0) {
    errflg++;
  }

  if (errflg) {
    printf("\nGeneral Linear Model: A simple program to run AN[C]OVAs. Supports all sort of\n");
    printf("experimental designs (at least those that I can imagine, balanced, unbalanced,\n");
    printf("nested, empty cells, missing treatment combinations, etc). The program automatically\n");
    printf("detects nested factors so there is no need to explicitly specify them in the model.\n");
    printf("It can handle up to nine factors and nine covariates.\n");
    printf("\nUsage: glm <options> (options in square brackets are optional)\n\n");
    printf(" -d  <string>    data file. A tab delimeted file with the following format:\n");
    printf("                 A header line with the sample names followed by one or more lines\n");
    printf("                 with the data. The first column of these lines must start with the\n");
    printf("                 descriptor name. Missing data must be entered as NA.\n");
    printf(" -f  <string>    factor file. A tab delimited file containing the classes for\n");
    printf("                 the samples.\n");
    printf(" -c  [<string>]  covariate file. A tab delimited file containing the covariates for\n");
    printf("                 the samples.\n");
    printf(" -o  [<string>]  output file. If no file is specified it defaults to stdout.\n");
    printf(" -s  [<integer>] number of samples.\n");
    printf(" -g  [<integer>] number of variables.\n");
    printf(" -n  <string>    variable name. Run glm only for the specified variable name\n");
    printf(" -m  <string>    fixed effect(s) name (must match one name in the factor file).\n");
    printf("                 Multiple fixed effects are entered with multiple -m in the command\n");
    printf("                 line. Do NOT enter interactions here.\n");
    printf(" -i  [<string>]  interaction(s) between factors (must match the names in the factor\n");
    printf("                 file). Multiple interactions are entered with multiple -i in the command\n");
    printf("                 line. Each interaction must be joined by an asterix '*', for example,\n");
    printf("                 subject*treatment.\n");
    printf(" -r  [<string>]  random factor name (must match one of the names in the factor file).\n");
    printf("                 Multiple random factors are entered with multiple -r in the command\n");
    printf("                 line.\n");
    printf(" -t  [<integer>] type of sum of squares [2|3]. Type III is implemented only\n");
    printf(" -q  [<float>]   q value to use in false discovery rate. (default is 0.05).\n\n");
    printf(" -p  [<integer>] number of random permutations (not implemented yet).\n\n");
    printf("Examples\n");
    printf("\tglm -d ./t/3f.2l.2l.2l.b.dat -f ./t/3f.2l.2l.2l.b.fac -m Grp1 -m Grp2 -m Grp3\n");
    printf("\tglm -d ./t/2f.6l.2l.w.dat -f ./t/2f.6l.2l.w.fac -m Grp -m Sub\n");
    printf("\tglm -d ./t/emptycells.dat -f ./t/emptycells.fac -i 'R*T'\n");
    printf("\tglm -d ./t/glm.dat -c ./t/glm.cov -f ./t/glm.fac -m CD -m FA -v -n Var1\n");
    printf("\n");
    exit(0);
  }
   
}

RESULT *
set_result(int ndata,
	   FACTOR *factor,
	   COVARIATE *covariate,
	   MODEL *model,
	   MATRIX *d,
	   double Q)
{

  int i, j, f, t, n, idx, len;
  char str[1000], stra[1000], strn[100];
  RESULT *res;
  VECTOR *vr;

  res = TYPE_ALLOC(RESULT);
  res->n = ndata;
  res->a = model->f + model->i + 5;
  res->attributes = ARRAY_ALLOC(res->a, char *);
  res->variables = ARRAY_ALLOC(res->n, char *);
  res->vectors = ARRAY_ALLOC(res->n, VECTOR *);
  for (i=0;i<ndata;i++) {
    vr = res->vectors[i] = TYPE_ALLOC(VECTOR);
    vr->data = ARRAY_ALLOC(res->a, double);
    vr->n = res->a;
  }

  idx=0;
  for (i=0;i<model->f;i++) {
    f = model->fixed[i];
    t = model->tyfixed[i];
    if (t==0) {
      // factor
      strcpy(str, factor->fnames[f]);
      set_attribute_result(idx,str,res);
    } else {
      // covariate
      strcpy(str, covariate->cnames[f]);
      set_attribute_result(idx,str,res);
    }
    idx++;
  }

  for (i=0;i<model->i;i++) {
    n=0;
    strcpy(str, "");
    for (j=0;j<model->numinter[i];j++) {
      f = model->inter[i][j];
      t = model->tyinter[i][j];
      if (t==0 && model->infactornest[f]>0) {
	strcpy(strn, factor->fnames[f]);
	n++;
      } else {
	if (t==1) {
	  // factor
	  strcat(str, covariate->cnames[f]); 
	} else {
	  // covariate
	  strcat(str, factor->fnames[f]);
	}
	strcat(str, "*");
      }
    }

    len = strlen(str)-1;
    strncpy(stra, str, len);
    stra[len]='\0';
    if (n>0) {
      strcat(strn, "(");
      strcat(strn, stra);
      strcat(strn, ")");
      set_attribute_result(idx,strn,res);
    } else {
      set_attribute_result(idx,stra,res);
    }
    idx++;
  }

  strcpy(str, "Model");
  set_attribute_result(idx,str,res);
  idx++;
  strcpy(str, "R-square");
  set_attribute_result(idx,str,res);
  idx++;
  strcpy(str, "Adj R-sq");
  set_attribute_result(idx,str,res);
  idx++;
  strcpy(str, "AIC");
  set_attribute_result(idx,str,res);
  idx++;
  strcpy(str, "FDR (");
  sprintf(strn, "%5.3f", Q);
  strcat(str, strn);
  strcat(str, ")");
  set_attribute_result(idx,str,res);
  idx++;

  return res;

}

void
set_attribute_result(int i,
		     char *str,
		     RESULT *res)
{

  int len;

  len = strlen(str);
  res->attributes[i] = ARRAY_ALLOC((len+1), char);  
  strncpy(res->attributes[i], str, len);
  res->attributes[i][len] = '\0';

}

void
mult_test_correction(RESULT *res,
		     MODEL *model,
		     double Q)
{

  int i, nsig;
  VECTOR *v;
  double data[res->n], qdata[res->n];  
  double cutoff;

  for (i=0;i<res->n;i++) {
    v = res->vectors[i];
    data[i] = v->data[model->t];
  }

  fdr(res->n,data,qdata,Q,1,&cutoff,&nsig);

  for (i=0;i<res->n;i++) {
    v = res->vectors[i];
    v->data[res->a-1] = qdata[i];
  }

}

RESULT *
calculate_glm(XDATA *xdata,
	      FACTOR *factor,
	      COVARIATE *covariate,
	      char *variable_name,
	      MODEL *model,
	      int ssq,
	      double Q,
	      int permuts)
{

    int i, j, k, len, n, g, s, na;
  int *variables, **indeces, *rank, *rna;
  RESULT *res;
  MATRIX *des, *gef, *X;
  MATRIX *desna, *gefna, *Xna;
  NMATRIX *hyp, *con;
  NMATRIX *hypna, *conna;
  VECTOR *v, *v1;

  // Allocate memory for an array with the indeces for the data
  // vectors. If thre is one variable the lenght will be one otherwise the
  // lenght will be the size of the number of variables in the xdata
  // object.
  if (variable_name != NULL) {
    n = 1;   
  } else {
    n = xdata->n;   
  }
  variables = ARRAY_ALLOC(n, int);  
      
  // The number of samples. If there is are factors we get the length
  // from the factor object otherwise we get it from the covariate
  // object
  s = factor->s > 0 ? factor->s : covariate->s;

  // An array to store the data with NA values so we can clean the
  // design matrix if necesary
  rna = ARRAY_ALLOC(s, int);

  // Generate some matrices used in all all calculations to speed up
  // the glm implementation (for data without NA values)
  des = set_design_matrix(factor, covariate, model, &indeces);

  gef = set_gef_matrix(des, &rank);

  if (ssq == 2) {
    hyp = set_ef3_matrix(model, gef, rank, indeces);
  } else {
    hyp = set_ef3_matrix(model, gef, rank, indeces);
  }

  // inv(X'X) * X'
  X = set_X_matrix(des);
 
  // C * inv(X'X) * C'
  con = set_contrast_matrices(des, hyp);

  // Load the variables array with the indeces of the data vectors
  if (variable_name != NULL) {
    g = search_variable(xdata, variable_name);
    if (g < 0) {
      fprintf(stderr, "Unable to find %s in data\n", variable_name);
      exit(1);
    }
    variables[0] = g;
  } else {
    for (i=0;i<n;i++) {
      variables[i]=i;
    }      
  }

  // Set the result object
  res = set_result(n,
		   factor,
		   covariate,
		   model,
		   des,
		   Q);
  
  // We used this v vector to copy the data
  v = TYPE_ALLOC(VECTOR);
  v->data = ARRAY_ALLOC(des->n, double);

  // Iterate over the variables (rows)
  for (i=0;i<n;i++) {

    // Set the index for the variable
    g = variables[i];

    // The size of the vector defaults to the size of the design matrix
    v->n = des->n;

    // Get the variable name
    len = strlen(xdata->variables[g]);
    res->variables[i] = ARRAY_ALLOC((len+1), char);
    strncpy(res->variables[i], xdata->variables[g], len);
    res->variables[i][len] = '\0';

    if (verbose || debug) printf("%s\n", xdata->variables[g]);

    // We clean the data if there are skipped samples or NA data
    v1 = xdata->vectors[g];
    k=0;
    na=0;
    for (j=0;j<s;j++) {
      if (model->include[j]>0) {
	if (! isnan(v1->data[j])) {
	  v->data[k] = v1->data[j];
	  rna[k]=j;
	  k++; 
	} else {
	  na++;
	}
      }
    }

    // Resize the vector if there is NA data
    v->n -= na;

    if (na > 0) {

      // We need to create new matrices if there is NA data
      if (debug) printf("Original Matrices Modified\n");
      desna =clean_design_matrix(rna, des, na); 
      gefna = set_gef_matrix(desna, &rank);
      if (ssq == 2) {
        hypna = set_ef3_matrix(model, gefna, rank, indeces);
      } else {
	hypna = set_ef3_matrix(model, gefna, rank, indeces);
      }
      Xna = set_X_matrix(desna);
      conna = set_contrast_matrices(desna, hypna);
      glm(v,desna,Xna,hypna,conna,res,model,i);
      free(desna);
      free(gefna);
      free(Xna);
      free(hypna);
      free(conna);

    } else {

      glm(v,des,X,hyp,con,res,model,i);

    } 

    free(v1);

  }

  mult_test_correction(res,model,Q);

  free(v);
  free(rna);
  free(rank);
  free(indeces);
  free(des);
  free(gef);
  free(X);
  free(hyp);
  free(con);

  return res;

}

int
search_variable(XDATA *xdata,
		char *variable_name)
{

  int i;
  
  for (i = 0; i < xdata->n; i++) {
    if (strcmp(variable_name, xdata->variables[i]) == 0) {
      return i;
    }
  }
  
  return (int) -1;

}

MATRIX *
set_X_matrix(MATRIX *des)
{

  int *rank, r;
  MATRIX *Xt, *XtX, *iXtX, *X;
  
  Xt = transpose(des,0);
  XtX = multiply(Xt,des,0);
  rank = ARRAY_ALLOC(XtX->m, int);
  iXtX = g2invert(XtX,
		  XtX->m,
		  &rank,
		  &r,
		  1);
  X = multiply(iXtX,Xt,1);
  free(rank);
  
  return X;

}

NMATRIX *
set_contrast_matrices(MATRIX *des,
		      NMATRIX *hyp)
{

  int i;
  int *rank, r;
  MATRIX *Xt, *XtX, *iXtX;
  MATRIX *Ct, *CiXtX, *CiXtXCt;
  NMATRIX *con;  

  // Memory allocation
  con = TYPE_ALLOC(NMATRIX);

  con->n = hyp->n;
  con->m = ARRAY_ALLOC(con->n, MATRIX *);

  // The design matrix
  Xt = transpose(des,0);
  XtX = multiply(Xt,des,0);
  rank = ARRAY_ALLOC(XtX->m, int);
  iXtX = g2invert(XtX,
		  XtX->m,
		  &rank,
		  &r,
		  1);
  free(Xt);
  free(rank);

  // The contrast matrices
  for (i=0;i<con->n;i++) {
    Ct = transpose(hyp->m[i],0);
    CiXtX = multiply(hyp->m[i],iXtX,0);
    CiXtXCt = multiply(CiXtX, Ct, 1);
    rank = ARRAY_ALLOC(CiXtXCt->m, int);
    con->m[i] = g2invert(CiXtXCt,
			 CiXtXCt->m,
			 &rank,
			 &r,
			 1);
    free(rank);
  }

  free(iXtX);

  return con;

}

MATRIX *
sas_ortogonalize(MATRIX *ef,
		 MATRIX *in)
{

  int i, j, ii;
  double *w, **v, ab, bb, r;
  VECTOR *wv;
  MATRIX *m, *wm;

  if (ef->n != in->n) {
    exit_on_error("not the same number of columns.");
  }

  // Memory Allocation and initialization
  // of the orthonormal matrix
  m = TYPE_ALLOC(MATRIX); 
  m->n = ef->n;
  m->m = ef->m;
  m->data = ARRAY_ALLOC(ef->n, double *);
  for (i=0;i<ef->n;i++) {
    m->data[i] = ARRAY_ALLOC(ef->m, double);
    for (j=0;j<ef->m;j++) {
      m->data[i][j] = ef->data[i][j];
    }
  }

  // Memory allocation for the SVD method
  w = ARRAY_ALLOC(in->n, double);
  v = ARRAY_ALLOC(in->m, double *);
  for (i=0;i<in->m;i++) {
    v[i] = ARRAY_ALLOC(in->m, double);
  }

  // SVD
  if (debug) {
    printf ("SAS ORTHOGONALIZATION\n\n");
  }
  svdcmp(in->data, in->n, in->m, w, v);
  if (debug) {
    printf ("\n");
  }

  // Convert arrays to objects
  wv = TYPE_ALLOC(VECTOR); 
  wv->n = in->m;
  wv->data = ARRAY_ALLOC(in->m, double);
  for (i=0;i<in->m;i++) {
    wv->data[i] = w[i];
  }

  wm = TYPE_ALLOC(MATRIX); 
  wm = diag(wv);
  in = multiply(in,wm,1);

  for (i=0;i<ef->m;i++) {
    for (ii=0;ii<in->m;ii++) {
      ab = dot(ef,in,i,ii);
      bb = dot(in,in,ii,ii);
      if (bb!=0) {
	r = ab/bb;
      } else {
	r = 0;
      }
      for (j=0;j<ef->n;j++) {
	m->data[j][i] -= (r * in->data[j][ii]);
	if ((double)(fabs(m->data[j][i]))<EPS) {
	  m->data[j][i] = 0;
	}
      }
    }
  }

  free(v);
  free(w);
  free(in);

  return m;

}

NMATRIX *
set_ef3_matrix(MODEL *model,
	       MATRIX *gef,
               int *rank,
               int **indeces)
{
	       
  int i, j, k, l, m, n, d, dof;
  MATRIX *ef, *in;
  NMATRIX *ef3;

  // Memory allocation
  ef3 = TYPE_ALLOC(NMATRIX);

  ef3->n = model->t+1;
  ef3->m = ARRAY_ALLOC(ef3->n, MATRIX *);
  
  // Set the matrices for the terms
  for (i=0;i<model->t;i++) {
    dof = 0;
    n = gef->m;
    m = indeces[i][1] - indeces[i][0];
    ef = mem_allocate_matrix(n, m);
    for (j=indeces[i][0];j<indeces[i][1];j++) {
      if (rank[j]>0) {
	for (k=0;k<ef->n;k++) {
	  ef->data[k][dof] = gef->data[k][j];
	}
	dof++;
      }
    }
    ef->m = dof;
    dof = 0;
    if (model->numdep[i]>0) {
      m = gef->n - 1;
      in = mem_allocate_matrix(n, m);
      for (j=0;j<model->numdep[i];j++) {
	d = model->dep[i][j];
	for (k=indeces[d][0];k<indeces[d][1];k++) {
	  if (rank[k]>0) {
	    for (l=0;l<in->n;l++) {
	      in->data[l][dof] = gef->data[l][k];
	    }
	    dof++;
	  }
	}
      }
      in->m = dof;
      ef = sas_ortogonalize(ef,in);
    }
    // Copy the data to the object
    ef3->m[i] = mem_allocate_matrix(ef->m, ef->n);
    for(j=0;j<ef->m;j++) {
      for(k=0;k<ef->n;k++) {
	ef3->m[i]->data[j][k] = ef->data[k][j];
      }
    }
  }
  // Set the matrix for the model
  ef3->m[model->t] = mem_allocate_matrix(ef->n-1, ef->n);
  for(i=0;i<ef->n-1;i++) {
    for(j=0;j<ef->n;j++) {
      if (i==j-1) {
	ef3->m[model->t]->data[i][j] = 1;
      } else {
	ef3->m[model->t]->data[i][j] = 0;
      }
    }
  }
  free(ef);

  if (debug) {
    printf ("TYPE III ESTIMABLE FUNCTIONS MATRICES\n\n");
    for (i=0;i<ef3->n;i++) {
      if (i == ef3->n-1) {
	printf ("MODEL MATRIX %iX%i\n\n", ef3->m[i]->n, ef3->m[i]->m);
      } else {
	printf ("TERM %i MATRIX %iX%i\n\n", i, ef3->m[i]->n, ef3->m[i]->m);
      }
      dump_matrix(ef3->m[i]);
      printf ("\n");    
    }
    printf("\n\n");
  }
  
  return ef3;
  
}

MATRIX *
set_gef_matrix(MATRIX *x,
	       int **rank)
{

    int r, i;
  MATRIX *m, *xt, *xtx, *ixtx;

  // General Form of Estimable
  // Functions according to SAS

  xt = transpose(x,0);
  xtx = multiply(xt,x,0);
  *rank = ARRAY_ALLOC(xtx->m, int);
  ixtx = g2invert(xtx,
		  xtx->m,
		  rank,
		  &r,
		  0);
  m = multiply(ixtx,xtx,0);
  m = transpose(m,1);

  if (debug) {
    printf ("GEN. ESTIMABLE FUNCTIONS (GEF) MATRIX %iX%i\n\n", m->n, m->m);
    dump_matrix(m);
    printf ("\n\n");    
    printf("RANK FOR GEF MATRIX:\n[ ");
    for (i=0;i<xtx->m;i++) {
      printf("%10i", (*rank)[i]);
    }
    printf(" ]\n\n");
  }

  free(xt);
  free(xtx);
  free(ixtx);
  
  return m;

}

MATRIX *
set_design_matrix(FACTOR *factor,
		  COVARIATE *covariate,
		  MODEL *model,
		  int ***indeces)
{

  int i, j, idx;
  int id, ty;
  MATRIX *m, *m1, *m2;

  j = model->f + model->i;
  *indeces = ARRAY_ALLOC(j, int *);
  for (i=0;i<j;i++) {
    (*indeces)[i] = ARRAY_ALLOC(2, int); 
  }
   
  // The first the column is the means (intercept)
  m = fill(model->s-model->k, 1, 1);

  // idx will have the locations of the effects
  // and interactions
  idx = 1;

  // We start with the fixed effects
  for (i=0;i<model->f;i++) {
    (*indeces)[i][0] = idx;
    id = model->fixed[i];
    ty = model->tyfixed[i];
    m1 = design(factor, covariate, model, id, ty);
    m = augment_cols(m, m1, 1);
    idx = m->m;
    (*indeces)[i][1] = idx;
  }  

  // We continue with the interactions
  for (i=0;i<model->i;i++) {
    m1 = fill(model->s-model->k, 1, 1);
    for (j=0;j<model->numinter[i];j++) {
      id = model->inter[i][j];
      ty = model->tyinter[i][j];
      m2 = design(factor, covariate, model, id, ty);
      m1 = cross(m1, m2, 1);
    }    
    (*indeces)[model->f+i][0] = idx;
    m = augment_cols(m, m1, 1);
    idx = m->m;
    (*indeces)[model->f+i][1] = idx;
  }
   
  if (debug) {
    printf ("DESIGN MATRIX %iX%i\n\n", m->n, m->m);
    dump_matrix(m);
    printf ("\n\n");    
    printf(" Indeces for DESIGN MATRIX:\n");
    for (i=0;i<model->f+model->i;i++) {
      printf(" T%i => [ %2i - %2i ]\n", i, (*indeces)[i][0], (*indeces)[i][1]-1);
    }
    printf("\n\n");
  }
  
  return m;

}

MATRIX *
clean_design_matrix(int *rna,
		    MATRIX *des,
		    int na)
{

  int i, j;
  MATRIX *desna;

  desna = mem_allocate_matrix(des->n - na, des->m);

  for (i=0;i<des->n-na;i++) {
    for (j=0;j<des->m;j++) {
	desna->data[i][j] = des->data[rna[i]][j];
    }
  }

  if (debug) {
    printf ("DESIGN MATRIX (no NA) %iX%i\n\n", desna->n, desna->m);
    dump_matrix(desna);
    printf ("\n\n");    
  }

  return desna;

}

MATRIX *
design(FACTOR *factor,
       COVARIATE *covariate,
       MODEL *model,
       int o,
       int t)
{

  MATRIX *d;

  int i, j, c=0;

  if (t>0) {
    d = fill(model->s, 1, 0);
    for (i=0;i<d->n;i++) {
      if (model->include[i]>0) {
	d->data[c][0] = covariate->values[o][i];
	c++;
      }
    }
  } else {
    d = fill(model->s, factor->classes[o], 0);
    for (i=0;i<d->m;i++) {
      c=0;
      for (j=0;j<d->n;j++) {
	if (model->include[j]>0) {
	  if (factor->factors[o][j]==i) {
	    d->data[c][i] = 1;
	  }
	  c++;
	}
      }
    }
  }
  // We need to adjust the number of columns
  d->n = c;

  return d;

}

int
is_nested(FACTOR *factor,
	  COVARIATE *covariate,
	  MODEL *model,
	  int f)
{

  MATRIX *f1, *f2;
  MATRIX *f1t, *f2t;
  MATRIX *f2tf1, *f1tf2;
  MATRIX *res;
  int i, j, k;
  int n=0, c=0, z=0;

  f1 = design(factor, covariate, model, f, 0);
  f1t = transpose(f1,0); 

  for (i=0;i<factor->n;i++) {
    c=0;
    if (model->infactor[i]>0 && i != f && factor->classes[i]<f1->m) {
      c++;
      f2 = design(factor, covariate, model, i, 0);
      f2t = transpose(f2,0); 
    }
    if (c>0) {
      f2tf1 = multiply(f2t,f1,0);
      f1tf2 = multiply(f1t,f2,0);
      res = multiply(f2tf1, f1tf2,1);
      z=0;
      for (j=0;j<res->n;j++) {
	for (k=0;k<res->m;k++) {
	  if (res->data[j][k]==0) {
	    z++;
	  }
	}
      }
      if (z>=(res->n*res->m/2)) {
	model->nested[f][model->numnest[f]]=i;
	model->numnest[f]++;
	model->denomnest[f]/=factor->classes[i];
	n=1;
      }
      free(f2);
      free(f2t);
      free(res);
    }
  }
  free(f1);
  free(f1t);

  return n;

}

MODEL *
set_model(FACTOR *factor,
	  COVARIATE *covariate,
	  char **fixed_effs,
	  int n_fixed_effs,
	  char **interacts,
	  int n_interacts,
	  char **random_facs,
	  int n_random_facs)
{

  int i, j, k, l, si, cnt, ok;
  int *tmpt, *tmpty, *tmpsrt;
  int *srtint, *lenint;
  char *tok, *ctok;
  MODEL *m;
  
  m = TYPE_ALLOC(MODEL);

  // A. We set the fixed effects. We (1)
  // check that the name of each fixed effect
  // is either a valid factor name or a valid 
  // covariate name and (2) find all the
  // factors and covariates included in the
  // model.  We also check the names of the
  // random factors


  // 1. Allocate memory
  m->infactor = ARRAY_ALLOC(factor->n, int);
  m->infactorrand = ARRAY_ALLOC(factor->n, int);
  m->infactornest = ARRAY_ALLOC(factor->n, int);


  // 2. Initialize
  for (i=0;i<factor->n;i++) {
    // No factors are included
    m->infactor[i]=0;
    m->infactorrand[i]=0;
    m->infactornest[i]=0;
  }
  if (covariate->n > 0) {
    m->incovariate = ARRAY_ALLOC(covariate->n, int);
    for (i=0;i<covariate->n;i++) {
      // No covariates are included
      m->incovariate[i]=0;
    }
  }

  // 3. For every fixed effect name we have
  // to check all the factor and all the
  // covariate names
  m->f = n_fixed_effs;
  tmpt = ARRAY_ALLOC(m->f, int);
  tmpty = ARRAY_ALLOC(m->f, int);
  for (i=0;i<n_fixed_effs;i++) {
    ok=0;
    // check the factors first
    for (j=0;j<factor->n;j++) {
      if (!strcmp(factor->fnames[j], fixed_effs[i])) {
	// .. get the factor Id
	tmpt[i]=j;
	// .. type of fixed effect
	tmpty[i]=0;
	// .. include factor in the model
	m->infactor[j]++;
	ok++;
	break;
      }
    }
    if (ok==0) {
      // then check the covariates
      if (covariate->n > 0) {
	for (j=0;j<covariate->n;j++) {
	  if (!strcmp(covariate->cnames[j], fixed_effs[i])) {
	    // .. get the covariate id
	    tmpt[i]=j;
	    // .. type of fixed effect
	    tmpty[i]=1;
	    // include covariate in the model
	    m->incovariate[j]++;
	    ok++;
	    break;
	  }
	}
      }
      if (ok==0) {
	fprintf(stderr, "Fixed effect %s not found in factor file\n", fixed_effs[i]);
	exit(1);
      }
    }
  }

  // 4. We then sort the fixed effects
  // and set them in the model
  m->fixed = ARRAY_ALLOC(m->f, int);
  m->tyfixed = ARRAY_ALLOC(m->f, int);
  tmpsrt = ARRAY_ALLOC(m->f, int);
  sort_by_indexi(m->f, tmpt, tmpsrt);
  for (i=0;i<m->f;i++) {
    si=tmpsrt[i];
    m->fixed[i]=tmpt[si];
    m->tyfixed[i]=tmpty[si];
  }
  free( tmpt);
  free( tmpty);
  free( tmpsrt);


  // 5. We now ckeck that (1) the names of
  // the random effects are valid factor names
  // and (2) they are included in the model
  m->r = n_random_facs;
  m->rand = ARRAY_ALLOC(m->r, int);
  for (i=0;i<n_random_facs;i++) {
    ok=0;
    for (j=0;j<factor->n;j++) {
      // check the factor name
      if (!strcmp(factor->fnames[j], random_facs[i])) {
	m->rand[i] = j;
	m->infactorrand[j]++;
	ok++;
	break;
      }
    }
    if (ok == 0) {
      fprintf(stderr, "Random factor %s not found in factor file\n", random_facs[i]);
      exit(1);
    }
  }


  // B. Next, we set the interactions. We need to
  // check (1) if the terms in each interaction
  // has a valid name, either a factor or a 
  // covariate valid name, (2) include the
  // term in the model and (3) find if we have
  // enough DOF to calculate the interaction 

  // 1. We start by sorting the interactions
  // by the number of terms 
  srtint = ARRAY_ALLOC(n_interacts, int);
  lenint = ARRAY_ALLOC(n_interacts, int);
  for (i=0;i<n_interacts;i++) {
    // find number the terms in the interaction
    cnt=0;
    ctok = strdup(interacts[i]);
    tok = strtok (ctok, "*");
    while (tok != NULL) {     
      cnt++;
      tok = strtok (NULL, "*");
    }
    free(ctok);
    lenint[i]=cnt;
  }
  sort_by_indexi(n_interacts, lenint, srtint);

  // 2. Allocate memory
  m->i = n_interacts;
  m->numinter = ARRAY_ALLOC(m->i, int);
  m->inter = ARRAY_ALLOC(m->i, int *);
  m->tyinter = ARRAY_ALLOC(m->i, int *);
  m->randinter = ARRAY_ALLOC(m->i, int);

  // 3. Check the interactions
  for (i=0;i<n_interacts;i++) {
    si=srtint[i];
    m->randinter[i] = 0;
    m->numinter[i] = lenint[si];
    tmpt = ARRAY_ALLOC(1000, int);
    tmpty = ARRAY_ALLOC(1000, int);
    // check the names of all term
    cnt=0; // term number
    tok = strtok (interacts[si], "*");
    while (tok != NULL) {     
      ok=0;
      // check factor names first ..
      for (j=0;j<factor->n;j++) {
	if (!strcmp(factor->fnames[j], tok)) {
	  // .. get the factor id
	  tmpt[cnt]=j;
	  // .. type of term
	  tmpty[cnt]=0;
	  // .. include factor in the model
	  m->infactor[j]++;
	  // .. check for random factors
	  if (m->infactorrand[j]>0) {
	    m->randinter[i]++;
	  }
	  ok++;
	  break;
	}
      }
      if (ok == 0) {
	// now the covariate names
	if (covariate->n > 0) {
	  for (j=0;j<covariate->n;j++) {
	    if (!strcmp(covariate->cnames[j], tok)) {
	      // get the covariate id
	      tmpt[cnt]=j;
	      // type of term 
	      tmpty[cnt]=1;
	      // include the covariate in the model
	      m->incovariate[j]++;
	      ok++;
	      break;
	    }
	  }
	}
	if (ok == 0) {
	  fprintf(stderr, "Term %s not found in factor or covariate names\n", tok);
	  exit(1);
	}
      }
      cnt++;
      tok = strtok (NULL, "*");
    }
    // We then sort the interaction
    // and set them in the model
    m->inter[i] = ARRAY_ALLOC(cnt, int);
    m->tyinter[i] = ARRAY_ALLOC(cnt, int);
    tmpsrt = ARRAY_ALLOC(cnt, int);
    sort_by_indexi(cnt, tmpt, tmpsrt);
    for (j=0;j<cnt;j++) {
      si=tmpsrt[j];
      m->inter[i][j]=tmpt[si];
      m->tyinter[i][j]=tmpty[si];
    }
    free(tmpt);
    free(tmpty);
    free(tmpsrt);
  }
  free(lenint);
  free(srtint);

  // C. We define the number of samples included in the
  // model based on the included factors and covariates
  // (those that are not 'Unassigned')


  // 1.  Allocate memory and set number of samples.
  // We start by inclding all samples
  m->s = factor->n>0 ? factor->s : covariate->s;
  m->include = ARRAY_ALLOC(m->s, int);
  for (i=0;i<m->s;i++) {
    m->include[i] = 1;
  }


  // 2. We begin with the covariates. We exclude any
  // samples with missing data for a covariate.
  for (i=0;i<covariate->n;i++) {
    if (m->incovariate[i]>0) {
      for (j=0;j<m->s;j++) {
	if (isnan(covariate->values[i][j])) {
	  m->include[j] = 0;
	}
      }
    }
  }


  // 3. We continue with the factors.
  for (i=0;i<factor->n;i++) {
    if (m->infactor[i]>0) {
      for (j=0;j<m->s;j++) {
	if (factor->factors[i][j]<0) {
	  m->include[j] = 0;
	}
      }
    }
  }


  // 4. We count the number of excluded samples
  m->k = 0;
  for (i=0;i<m->s;i++) {
    if (m->include[i] < 1) {
      m->k++;
    }
  }


  // E. We check if we have nested factors so
  // we can enable back some of the required
  // interactions 
  m->n=0;
  m->infactornest = ARRAY_ALLOC(factor->n, int);
  m->numnest = ARRAY_ALLOC(factor->n, int);
  m->denomnest = ARRAY_ALLOC(factor->n, int);
  m->nested = ARRAY_ALLOC(factor->n, int *);
  for (i=0;i<factor->n;i++) {
    m->infactornest[i]=0;
    m->numnest[i]=0;
    m->denomnest[i]=factor->classes[i];
    m->nested[i] = ARRAY_ALLOC(factor->n, int);
    if (m->infactor[i]) {
      if (is_nested(factor, covariate, m, i)) {
	m->infactornest[i]++;
	m->n++;
      }
    }
  }


  // F. We finally check for dependencies
  // for the fixed effects and the 
  // interactions.


  // 1. Total number of terms in the model
  m->t = m->f + m->i;


  // 2. Allocate memory
  m->numdep = ARRAY_ALLOC(m->t, int);
  m->randdep = ARRAY_ALLOC(m->t, int);
  m->dep = ARRAY_ALLOC(m->t, int *);
  for (i=0;i<m->t;i++) {
    m->dep[i] = ARRAY_ALLOC(m->i, int);
  }


  // 3. Start with the fixed effects
  for (i=0;i<m->f;i++) {
    m->numdep[i]=0;
    for (j=0;j<m->i;j++) {
      for (k=0;k<m->numinter[j];k++) {
	if (m->fixed[i] == m->inter[j][k] &&
	    m->tyfixed[i] == m->tyinter[j][k]) {
	  m->dep[i][m->numdep[i]]=m->f+j;
	  m->numdep[i]++;
	}
      }
    }
    if (m->infactorrand[m->fixed[i]]>0) {
      m->randdep[i] = 1;
    } else {
      m->randdep[i] = 0;
    }
  }

  // 4. We continue with the interactions
  for (i=0;i<m->i;i++) {
    m->numdep[m->f+i]=0;
    for (j=i+1;j<m->i;j++) {
      if (m->numinter[j]>m->numinter[i]) {
	ok=0;
	for (k=0;k<m->numinter[i];k++) {
	  for (l=0;l<m->numinter[j];l++) {
	    if (m->inter[i][k] == m->inter[j][l] &&
		m->tyinter[i][k] == m->tyinter[j][l]) {
	      ok++;
	    }
	  }
	}
	if (ok==m->numinter[i]) {
	  m->dep[m->f+i][m->numdep[m->f+i]]=m->f+j;
	  m->numdep[m->f+i]++;
	}
      }
    }
    if (m->randinter[i]>0) {
      m->randdep[m->f+i] = 1;
    } else {
      m->randdep[m->f+i] = 0;
    }
  }

  // We are done!
  if (debug) {
    printf ("MODEL\n\n");
    printf (" Factors included in the model\n [ ");
    for (i=0;i<factor->n;i++) {
      if (m->infactor[i]>0) {
	printf("%3i ", i);
      }
    }
    printf (" ]\n\n");
    printf (" Factor names\n");
    for (i=0;i<factor->n;i++) {
      if (m->infactor[i]>0) {
	printf("%3i: %s\n", i, factor->fnames[i]);
      }
    }
    printf ("\n");
    if (covariate->n > 0) {
      printf (" Covariates included in the model\n [ ");
      for (i=0;i<covariate->n;i++) {
	if (m->incovariate[i]>0) {
	  printf("%3i ", i);
	}
      }
      printf (" ]\n\n");
      printf (" Covariate names\n");
      for (i=0;i<covariate->n;i++) {
	if (m->incovariate[i]>0) {
	  printf("%3i: %s\n", i, covariate->cnames[i]);
	}
      }
      printf ("\n");
    }
    printf (" Fixed Effects\n");
    printf (" Id(s) for fixed effects => ");
    print_array_i(m->fixed, m->f);
    printf (" Type of fixed effects   => ");
    print_array_i(m->tyfixed, m->f);
    printf ("\n");
    if (m->r > 0) {
      printf (" Random Effects\n");
      printf (" Id(s) for random effects => ");
      print_array_i(m->rand, m->r);
      printf ("\n");
    }
    if (m->i > 0) {
      printf (" Interactions\n");
      for (i=0;i<m->i;i++) {
	printf (" Ids for terms in I%-2i => ", i);
	printf (" ");
	print_array_i(m->inter[i], m->numinter[i]);	
	printf (" Type of terms in I%-2i => ", i);
	printf (" ");
	print_array_i(m->tyinter[i], m->numinter[i]);	
      }
      printf ("\n");
      printf (" Dependencies on terms\n");
      for (i=0;i<m->t;i++) {
	printf (" T%-2i => ", i); 
	if (m->numdep[i] > 0) {
	  printf (" ");
	  print_array_i(m->dep[i], m->numdep[i]);	
	} else {
	  printf (" [ None ]\n");	  
	}
      }
      printf ("\n");
    }
    if (m->n > 0) {
      printf (" Nested Factors\n");
      for (i=0;i<factor->n;i++) {
	if (m->infactornest[i]>0) {
	  printf (" F%-2i nested in Factor(s) ", i); 
	  print_array_i(m->nested[i], m->numnest[i]);	
	  printf (" F%-2i denominator %i\n", i, m->denomnest[i]); 
	}
      }
      printf ("\n");
    }
    printf (" Number of Samples skipped: %i\n", m->k);
    printf ("\n");
    printf (" Samples included\n");
    printf (" ");
    print_array_i(m->include, m->s);
    printf ("\n");
    printf ("\n");
  }

  return m;

}

double
get_sum_of_squares(MATRIX *H,
		   MATRIX *C,
		   MATRIX *b)
{

  double ss;
  MATRIX *Hb, *Hbt, *HbtC, *HbtCHb;


  Hb = multiply(H,b,0);
  Hbt = transpose(Hb,0);
  HbtC = multiply(Hbt, C, 0);
  HbtCHb = multiply(HbtC, Hb, 1);
  ss = HbtCHb->data[0][0];
  free(Hbt);
  free(HbtCHb);

  return ss;

}

void
parameter_estimates(VECTOR *Y,
		    MATRIX *des,
		    MATRIX *X,
		    MATRIX **b,
		    double *sse)
{

  MATRIX *y, *yt, *yty, *yhat, *ytyhat;

  // X is inv(X'X) * X'

  y = vector_to_matrix(Y,0);
  yt = transpose(y,0);
  yty = multiply(yt,y,0);
  *b = multiply(X,y,0);
  // the expected values
  yhat = multiply(des,*b,0);
  ytyhat = multiply(yt,yhat,1);
  *sse = yty->data[0][0]-ytyhat->data[0][0];
  free(y);
  free(yty);
  free(ytyhat);

}

void
satterthwaite_approx (MODEL *model,
		      double **temp,
		      double *df,
		      double *f,
		      double *ms,
		      int idx)
{

  int i, ii;
  double num=0, denom=0;

  // df are in temp[x][1]
  // ms are in temp[x][2]

  for (i=0;i<model->numdep[idx];i++) {
    ii = model->dep[idx][i];
    num += temp[ii][2];
    denom += pow(temp[ii][2],2)/temp[ii][1];
  }
  num = pow(num, 2);
  *df = num/denom;
  *f = temp[idx][2]/num;
  *ms = temp[idx][2];

  //ii = model->dep[idx][model->numdep[idx]-1];
  //*df = temp[ii][1];
  //*f = temp[idx][2]/temp[ii][2];
  //*ms = temp[ii][2];

}

void
glm(VECTOR *Y,
    MATRIX *design,
    MATRIX *X,
    NMATRIX *H,
    NMATRIX *C,
    RESULT *res,
    MODEL *model,
    int idx)
{

  int i, n, index;
  double ss, tss, ms, msm, df, sse, mse, dfe, f, pval, dft;
  double r2, ar2, aic;
  double **temp;
  VECTOR *vr;
  MATRIX *b;

  // To store results
  vr = res->vectors[idx];

  // the parameter estimates
  n = Y->n;
  parameter_estimates(Y, design, X, &b, &sse);

  // We iterate over the hypothesis
  // to find the total DOFand the error
  // DOF
  df = 0;
  for (i=0;i<H->n-1;i++) {
    df += H->m[i]->n;
  }
  dft = n - 1;
  dfe = dft - df;

  // Calculate glm for the model
  index = H->n-1;
  df = dft - dfe;
  ss = get_sum_of_squares(H->m[H->n-1],
			  C->m[C->n-1],
			  b);

  // The error terms
  tss = sse + ss;
  mse = sse / dfe;
  ms = msm = ss/df;
  f = fabs(ms/mse);
  pval = betai(0.5*dfe,
	       0.5*df,
	       dfe/(dfe+df*f));
  vr->data[index++] = pval;
  r2 = ss/(tss);
  vr->data[index++] = r2;
  ar2 = MAX(0,1 - ((dft/(dft-df))*(1-r2)));
  vr->data[index++] = ar2;
  aic = n * log(sse/n) + (2*(b->n+1));
  vr->data[index] = aic;

  // Calculate glm for the terms in the model
  // and store the results in a temporary array
  index = 0;
  temp = ARRAY_ALLOC(H->n-1, double *);
  for (i=0;i<H->n-1;i++) {
    temp[i] = ARRAY_ALLOC(5, double);
    temp[i][0] = get_sum_of_squares(H->m[i],
				    C->m[i],
				    b);
    temp[i][1] = H->m[i]->n;
    temp[i][2] = temp[i][0]/temp[i][1];
    temp[i][3] = fabs(temp[i][2]/mse);
    temp[i][4] = betai(0.5*dfe,
		       0.5*temp[i][1],
		       dfe/(dfe+temp[i][1]*temp[i][3]));
    vr->data[index++] = temp[i][4];
  }

  if (verbose || debug) {
    printf("\nANALYSIS OF VARIANCE                 SS     df           MS         F    Pr(>F)\n");
    printf ("%-24s %14.2f %6.0f %12.2f %9.3f  %7.4g", res->attributes[H->n-1], ss, df, ms, f, pval);      
    if (pval < 0.05) {
      printf (" *\n");
    } else {
      printf ("\n");
    }
    printf("Within (error)           %14.2f %6.0f %12.2f\n", sse, dfe, mse);
    printf("Total                    %14.2f %6.0f\n", tss, dft);
    printf("\nR2 = %6.4f     Adj.R2 = %6.4f   AIC = %6.4f\n", r2, ar2, aic);
    printf ("\nSOURCE OF VARIATION                  SS     df           MS         F    Pr(>F)\n");
    for (i=0;i<H->n-1;i++) {
      printf ("%-24s %14.2f %6.0f %12.2f %9.3f  %7.4g", res->attributes[i], temp[i][0], temp[i][1], temp[i][2], temp[i][3], temp[i][4]);      
      if (temp[i][4] < 0.05) {
	printf (" *\n");
      } else {
	printf ("\n");
      }
    }
    printf ("\n");
  }

  if (model->r > 0) {
    index=0;
    for (i=0;i<H->n-1;i++) {
      if (model->randdep[i]>0) {
	if (model->numdep[i]==0) {
	  temp[i][1] = df = dfe;
	  temp[i][3] = f = temp[i][2]/mse;
	} else {
	  satterthwaite_approx(model, temp, &df, &f, &ms, i);
	  temp[i][1] = df;
	  temp[i][2] = ms;
	  temp[i][3] = f;
	  temp[i][4] = pval = betai(0.5*dfe,
				    0.5*df,
				    dfe/(dfe+df*f));
	  vr->data[index++] = temp[i][4];
	}
      }
    }
    if (verbose || debug) {
      printf ("\nMODIFIED SOURCE OF VARIATION         ERR df     ERR MS      F    Pr(>F)\n");
      for (i=0;i<H->n-1;i++) {
	printf ("%-34s %6.0f %12.2f %9.3f  %7.4g", res->attributes[i], temp[i][1], temp[i][2], temp[i][3], temp[i][4]);	
	if (temp[i][4] < 0.05) {
	  printf (" *\n");
	} else {
	  printf ("\n");
	}
      }
      printf ("\n");
    }
  }

  free(temp);
  free(b);

}
	
int 
main(int argc,
     char **argv)

{

  char *data_file, *factor_file, *covariate_file, *output_file, *variable_name;
  int  samples, variables, n_fixed_effs, n_interacts, n_random_facs, ssq, permuts, lines;
  char **fixed_effs, **interacts, **random_facs; 
  double Q;
  XDATA *xdata;
  FACTOR *factor;
  COVARIATE *covariate;
  MODEL *model;
  RESULT *res;
  
  verbose = 0;

  debug = 0;

  parse_arguments(argc, argv,
		  &data_file,
		  &factor_file, &covariate_file,
		  &output_file,
		  &samples, &variables,
		  &variable_name,
		  &fixed_effs, &n_fixed_effs,
		  &interacts, &n_interacts,
		  &random_facs, &n_random_facs,
		  &ssq, &Q, &permuts);
    
  if (factor_file != NULL) {

    if (samples == 0) {
      samples = count_fields(factor_file);
    }

    lines = count_lines(factor_file);
    lines--;
  
  } else {

    lines = 0;

  }

  factor = read_factor(factor_file, lines, samples);

  if (covariate_file != NULL) {

    if (samples == 0) {
      samples = count_fields(covariate_file);
    }

    lines = count_lines(covariate_file);
    lines--;
  
  } else {

    lines = 0;

  }

  covariate = read_covariate(covariate_file, lines, samples);

  model = set_model(factor, covariate,
		    fixed_effs, n_fixed_effs,
		    interacts, n_interacts,
		    random_facs, n_random_facs);

  if (variables == 0) {
    variables = count_lines(data_file);
    variables--;
  }
  
  xdata = read_xdata(data_file, samples, variables);

  res = calculate_glm(xdata, factor, covariate, variable_name,
        	      model, ssq, Q, permuts);

  print_results(output_file, res, 1, model->t);
  
  free(data_file);
  free(factor_file);
  free(covariate_file);
  free(output_file);
  free(variable_name);
  free(fixed_effs);
  free(interacts);
  free(random_facs); 
  free(xdata);
  free(factor);
  free(covariate);
  free(model);
  free(res);

  exit(0);

}
