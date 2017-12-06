/*
*    This file is part of SEED.
*
*    Copyright (C) 2017, Caflisch Lab, University of Zurich
*
*    SEED is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    SEED is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

/* -->> unused variables !
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
*/

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */


void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
/*clangini start*/
double **subdmatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
		  long newrl, long newcl);
/*clangini end*/
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

/* ------ added by NM ------ */
float ***ppfvector(long nl, long nh);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ppfvector(float ***v, long nl, long nh);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
	           long ndl, long ndh);
/* ------------------------- */

/*clangini start*/
double ***ppdvector(long nl, long nh);
void free_ppdvector(double ***v, long nl, long nh);
/*clangini end*/

/* ------ added by MS ------ */
struct point *structpointvect(long nl, long nh);
char ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_structpointvect(struct point *v, long nl, long nh);
void free_c3tensor(char ***t, long nrl, long nrh, long ncl, long nch,
	           long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	           long ndl, long ndh);
/* ------------------------- */

/* added by FD .... test */
float * fvecresize(float *,int newsize);
int   * ivecresize(int   *,int newsize);
/* ... */
/*clangini START*/
double * dvecresize(double *,int newsize);
/*clangini END*/
/* ------ added by clangini ------ */
double **zero_dmatrix(int rl,int rh,int cl,int ch);
void dmm_prod(int m1,int n1,int m2,int n2,double **A,double **B,double **P);
void dm_transpose(int m, int n, double **M, double **T);
/* ------------------------- */

#else /* ANSI */
/* traditional - K&R */

void nrerror();
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();

/* ------ added by NM ------ */
float ***ppfvector();
char **cmatrix();
int ***i3tensor();
void free_ppfvector();
void free_cmatrix();
void free_i3tensor();
/* ------------------------- */

/* ------ added by MS ------ */
struct point *structpointvect();
char ***c3tensor();
double ***d3tensor();
void free_structpointvect();
void free_c3tensor();
void free_d3tensor();
/* ------------------------- */

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */
