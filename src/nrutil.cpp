#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <vector> // clangini
#include <iostream> // clangini
#define NR_END 1
#define FREE_ARG char*

void nrerror(char const *error_text)
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"\nMemory allocation error in SEED numerical recipes \"nrutil.c\": not enough memory?\n");
  /* fprintf(stderr,"Numerical Recipes run-time error...\n"); */
  fprintf(stderr,"%s\n",error_text);

  fprintf(stderr,"\nNumber of fragments that passed the energy cut-off might be too large.\n");
  fprintf(stderr,"Try reducing the number of residues(i10) or the energy cut-off(i13).\n");

  fprintf(stderr,"\n...now exiting to system...\n");

  exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  if(nh>=nl){
    float *v;

    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl+NR_END;
  }

  return NULL;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  if(nh>=nl){
    int *v;

    v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl+NR_END;
  }

  return NULL;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  if(nh>=nl){
    unsigned char *v;

    v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
    if (!v) nrerror("allocation failure in cvector()");
    return v-nl+NR_END;
  }
  return NULL;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  if(nh>=nl){
    unsigned long *v;

    v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
    if (!v) nrerror("allocation failure in lvector()");
    return v-nl+NR_END;
  }
  return NULL;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  if(nh>=nl){
    double *v;

    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl+NR_END;
  }
  return NULL;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  if(nrh>=nrl && nch>=ncl)
    {
      long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
      float **m;

      /* allocate pointers to rows */
      m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
      if (!m) nrerror("allocation failure 1 in matrix()");
      m += NR_END;
      m -= nrl;

      /* allocate rows and set pointers to them */
      m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
      if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
      m[nrl] += NR_END;
      m[nrl] -= ncl;

      for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

      /* return pointer to array of pointers to rows */
      return m;
    }
  return NULL;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  if(nrh>=nrl && nch>=ncl)
    {
      long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
      double **m;

      /* allocate pointers to rows */
      m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
      if (!m) nrerror("allocation failure 1 in matrix()");
      m += NR_END;
      m -= nrl;

      /* allocate rows and set pointers to them */
      m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
      if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
      m[nrl] += NR_END;
      m[nrl] -= ncl;

      for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

      /* return pointer to array of pointers to rows */
      return m;
    }
  return NULL;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  if(nrh>=nrl && nch>=ncl)
    {
      long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
      int **m;

      /* allocate pointers to rows */
      m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
      if (!m) nrerror("allocation failure 1 in matrix()");
      m += NR_END;
      m -= nrl;


      /* allocate rows and set pointers to them */
      m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
      if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
      m[nrl] += NR_END;
      m[nrl] -= ncl;

      for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

      /* return pointer to array of pointers to rows */
      return m;
    }
  return NULL;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
		  long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;

  /* allocate array of pointers to rows */
  m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/*clangini start*/
double **subdmatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
		  long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  double **m;

  /* allocate array of pointers to rows */
  m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}
/*clangini end*/

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
   declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
   and ncol=nch-ncl+1. The routine should be called with the address
   &a[0][0] as the first argument. */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;

  /* allocate pointers to pointers to rows */
  t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  if(v!=NULL)
    {
      free((FREE_ARG) (v+nl-NR_END));
    }
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  if(v!=NULL)
    {
      free((FREE_ARG) (v+nl-NR_END));
    }
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
  if(v!=NULL)
    {
      free((FREE_ARG) (v+nl-NR_END));
    }
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
  if(v!=NULL)
    {
      free((FREE_ARG) (v+nl-NR_END));
    }
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  if(v!=NULL)
    {
      free((FREE_ARG) (v+nl-NR_END));
    }
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  if(m!=NULL)
    {
      free((FREE_ARG) (m[nrl]+ncl-NR_END));
      free((FREE_ARG) (m+nrl-NR_END));
    }
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  if(m!=NULL)
    {
      free((FREE_ARG) (m[nrl]+ncl-NR_END));
      free((FREE_ARG) (m+nrl-NR_END));
    }
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  if(m!=NULL)
    {
      free((FREE_ARG) (m[nrl]+ncl-NR_END));
      free((FREE_ARG) (m+nrl-NR_END));
    }
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free a float 3tensor allocated by f3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}

/* ------ added by NM ------ */
float ***ppfvector(long nl, long nh)
/* allocate a vector of pointer to pointer of float with subscript
   range v[nl..nh] */
{
  float ***v;

  v=(float ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float**)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

/* clangini start */
double ***ppdvector(long nl, long nh)
/* allocate a vector of pointer to pointer of double with subscript
   range v[nl..nh] */
{
  double ***v;

  v=(double ***)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double**)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}
/* clangini end */
char **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a character matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  char **m;

  /* allocate pointers to rows */
  m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char*)));
  if (!m) nrerror("allocation failure 1 in cmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
  if (!m[nrl]) nrerror("allocation failure 2 in cmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate an integer 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  int ***t;

  /* allocate pointers to pointers to rows */
  t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int**)));
  if (!t) nrerror("allocation failure 1 in i3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
  if (!t[nrl]) nrerror("allocation failure 2 in i3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in i3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_ppfvector(float ***v, long nl, long nh)
/* free a vector of pointer to pointer of float allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

/* clangini start */
void free_ppdvector(double ***v, long nl, long nh)
/* free a vector of pointer to pointer of double allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}
/* clangini end */

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free a character matrix allocated by cmatrix() */
{
  if(m!=NULL)
    {
      free((FREE_ARG) (m[nrl]+ncl-NR_END));
      free((FREE_ARG) (m+nrl-NR_END));
    }
}

void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free an integer 3tensor allocated by i3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
/* ------------------------- */

/* ------ added by MS ------ */
struct point {
  double x;
  double y;
  double z;
};
struct point *structpointvect(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  struct point *v;

  v=(struct point *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(struct point)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

char ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate an integer 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  char ***t;

  /* allocate pointers to pointers to rows */
  t=(char ***) malloc((size_t)((nrow+NR_END)*sizeof(char**)));
  if (!t) nrerror("allocation failure 1 in c3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(char **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char*)));
  if (!t[nrl]) nrerror("allocation failure 2 in c3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(char *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(char)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in c3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_structpointvect(struct point *v, long nl, long nh)
/* free a struct point allocated with structpointvect() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_c3tensor(char ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free an integer 3tensor allocated by i3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free a double 3tensor allocated by d3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
/* ------------------------- */

float * fvecresize(float * vec, int newsize)
{
  float * tmp;
  /* fprintf(stderr,"FVEC-RESIZE %d\n",newsize); */
  /*tmp = realloc(vec,(newsize * sizeof(float)) ); clangini */
  tmp = (float *) realloc(vec,(newsize * sizeof(float)) ); /* clangini */

  if(tmp != NULL)
    { /* vec = tmp;  */
      return tmp;
    }

  printf("Memory allocation error while increasing memory usage\n");
  nrerror("allocation failure in vector()");
  exit(13);
  return NULL;
}
/*clangini START*/
double * dvecresize(double * vec, int newsize)
{
  double * tmp;
  /* fprintf(stderr,"FVEC-RESIZE %d\n",newsize); */
  /*tmp = realloc(vec,(newsize * sizeof(double)) ); clangini */
  tmp = (double *) realloc(vec,(newsize * sizeof(double)) ); /* clangini */

  if(tmp != NULL)
    { /* vec = tmp;  */
      return tmp;
    }

  printf("Memory allocation error while increasing memory usage\n");
  nrerror("allocation failure in dvector()");
  exit(13);
  return NULL;
}
/*clangini END*/
int * ivecresize(int * vec, int newsize)
{
  int * tmp;
  // tmp = realloc(vec,(newsize * sizeof(int)) ); clangini
  tmp = (int *) realloc(vec,(newsize * sizeof(int)) ); // clangini
  /* fprintf(stderr,"IVEC-RESIZE %d\n",newsize); */
  if(tmp != NULL)
    {
      /* vec = tmp;   */
      return tmp;
    }


  printf("Memory allocation error while increasing memory usage\n");
  nrerror("allocation failure in vector()");
  exit(13);

  return NULL;
}

/* Added by clangini */
double **zero_dmatrix(int rl,int rh,int cl,int ch){
/* Create matrix with range [rl..rh][cl..ch] whose entries are all zero. */
  double **A;
  int i,j;

  A = dmatrix(rl,rh,cl,ch);
  for (i=rl;i<=rh;i++)
    for (j=cl;j<=ch;j++)
      A[i][j] = 0.0;
  return A;
}

void dmm_prod(int m1,int n1,int m2,int n2,double **A,double **B,double **P){
/*
* Matrix product.  A is a m1 X n1 matrix with range [1..m1][1..n1]
*  and B is a m2 X n2 matrix with range [1..m2][1..n2].  n1 = m2.
*  P = A * B.
*/
  int i,j,k;

  for (i=1;i<=m1;i++){
    for (j=1;j<=n2;j++){
      P[i][j] = 0.0;
      for (k=1;k<=n1;k++){
        P[i][j] = P[i][j] + A[i][k] * B[k][j];
      }
    }
  }
}

void dm_transpose(int m, int n, double **M,double **T){
  /*
  * Create a Matrix T transpose of matrix M.
  * M is a m X n matrix with range [1..m][1..n]
  */
  int i,j;

  for (i=1;i<=m;i++){
    for (j=1;j<=n;j++){
      T[j][i] = M[i][j];
    }
  }

}

std::vector<double> lineseq(double start, double end, double step){
  int nstep = (int) (end - start)/step + 1;
  std::vector<double> sequence(nstep);

  for (int i = 0; i < nstep; i++)
    sequence[i] = start + i * step;
  return sequence;
}
/*---------------------*/

#else /* ANSI */
/* traditional - K&R */

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
     char error_text[];
     /* Numerical Recipes standard error handler */
{
  void exit();

  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

float *vector(nl,nh)
     long nh,nl;
     /* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

int *ivector(nl,nh)
     long nh,nl;
     /* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
     long nh,nl;
     /* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;

  v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
     long nh,nl;
     /* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;

  v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v) nrerror("allocation failure in lvector()");
  return v-nl+NR_END;
}

double *dvector(nl,nh)
     long nh,nl;
     /* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
     long nch,ncl,nrh,nrl;
     /* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((unsigned int)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
     long nch,ncl,nrh,nrl;
     /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix(nrl,nrh,ncl,nch)
     long nch,ncl,nrh,nrl;
     /* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
     float **a;
     long newcl,newrl,oldch,oldcl,oldrh,oldrl;
     /* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;

  /* allocate array of pointers to rows */
  m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
     float *a;
     long nch,ncl,nrh,nrl;
     /* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
	declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
	and ncol=nch-ncl+1. The routine should be called with the address
	&a[0][0] as the first argument. */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
  if (!m)	nrerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
     long nch,ncl,ndh,ndl,nrh,nrl;
     /* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;

  /* allocate pointers to pointers to rows */
  t=(float ***) malloc((unsigned int)((nrow+NR_END)*sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_vector(v,nl,nh)
     float *v;
     long nh,nl;
     /* free a float vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
     int *v;
     long nh,nl;
     /* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
     long nh,nl;
     unsigned char *v;
     /* free an unsigned char vector allocated with cvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
     long nh,nl;
     unsigned long *v;
     /* free an unsigned long vector allocated with lvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
     double *v;
     long nh,nl;
     /* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
     float **m;
     long nch,ncl,nrh,nrl;
     /* free a float matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
     double **m;
     long nch,ncl,nrh,nrl;
     /* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
     int **m;
     long nch,ncl,nrh,nrl;
     /* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
     float **b;
     long nch,ncl,nrh,nrl;
     /* free a submatrix allocated by submatrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
     float **b;
     long nch,ncl,nrh,nrl;
     /* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
     float ***t;
     long nch,ncl,ndh,ndl,nrh,nrl;
     /* free a float f3tensor allocated by f3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}


#endif /* ANSI */
