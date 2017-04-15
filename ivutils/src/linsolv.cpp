#include"linsolv.hpp"     
#include<complex>

using namespace std;

template int linsys(float *matr,float *vect,int num,float *x,int lda);
template int linsys(double *matr,double *vect,int num,double *x,int lda);
//template int linsys(complex<float> *matr,complex<float> *vect,int num,complex<float> *x,int lda);
//template int linsys(complex<double> *matr,complex<double> *vect,int num,complex<double> *x,int lda);

#ifdef USE_LAPACK
int inverse(complex<double>* A, int N){
  int *ipiv = new int[N+1];
  int lwork = N*N;
  complex<double> *work = new complex<double>[lwork];
  int info;
#ifndef USE_MKL
  zgetrf_(&N,&N,A,&N,ipiv,&info);
  zgetri_(&N,A,&N,ipiv,work,&lwork,&info);
#else
  zgetrf_(&N,&N,(MKL_Complex16 *)A,&N,ipiv,&info);
  zgetri_(&N,(MKL_Complex16 *)A,&N,ipiv,(MKL_Complex16 *)work,&lwork,&info);
#endif
  delete ipiv;
  delete work;

  return info ? -1 : 1;
}
#endif

template class linear_solver<double>;

template class newton_raphson<double, linear_solver<double> >;

#ifdef USE_LAPACK
template<>
void tgesv_(int *N, int *nrhs, double *A, int *lda,int *ipiv, double *b, int *ldb, int *info){
  dgesv_(N,nrhs,A,lda,ipiv,b,ldb,info);
}

template<>
void theev_(char *jobz,char *uplo,int *N,complex<double> *matr,int *lda,double *x,complex<double> *work,int *lwork,
  double *rwork,int *info){
#ifndef USE_MKL
  zheev_(jobz,uplo,N,matr,lda,x,work,lwork,rwork,info);
#else
  zheev_(jobz,uplo,N,(MKL_Complex16 *)matr,lda,x,(MKL_Complex16 *)work,lwork,rwork,info);
#endif
}


#endif

#ifdef USE_PARDISO
int pardiso(int n, double *a, int *ia, int *ja, double *vect, double *x){

  int mtype = 11; /* Real unsymmetric matrix */
  int nrhs = 1; /* Number of right hand sides. */

  /* Internal solver memory pointer pt, */
  void *pt[64];
  /* Pardiso control parameters. */
  int iparm[64];
  int maxfct, mnum, phase, error, msglvl;

  double ddum; /* Double dummy */
  int idum; /* Integer dummy. */

  for (int i = 0; i < 64; i++)
    iparm[i] = 0;

  iparm[0] = 1; /* No solver default */
  iparm[1] = 2; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm[2] = 1;
  iparm[3] = 0; /* No iterative-direct algorithm */
  iparm[4] = 0; /* No user fill-in reducing permutation */
  iparm[5] = 0; /* Write solution into x */
  iparm[6] = 0; /* Not in use */
  iparm[7] = 2; /* Max numbers of iterative refinement steps */
  iparm[8] = 0; /* Not in use */
  iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0; /* Not in use */
  iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  iparm[13] = 0; /* Output: Number of perturbed pivots */
  iparm[14] = 0; /* Not in use */
  iparm[15] = 0; /* Not in use */
  iparm[16] = 0; /* Not in use */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1; /* Output: Mflops for LU factorization */
  iparm[19] = 0; /* Output: Numbers of CG Iterations */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 0; /* Print statistical information in file */
  error = 0; /* Initialize error flag */

  for (int i = 0; i < 64; i++)
    pt[i] = 0;

  phase = 11;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, a, ia, ja, &idum, &nrhs,iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during symbolic factorization: %d", error);
    return -1;
  }

  phase = 22;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, a, ia, ja, &idum, &nrhs,iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
	return -1;
  }

  phase = 33;
  iparm[7] = 2; /* Max numbers of iterative refinement steps. */
  /* Set right hand side to one. */

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, a, ia, ja, &idum, &nrhs,iparm, &msglvl, vect, x, &error);
  if (error != 0) {
    printf("\nERROR during solution: %d", error);
    return -1;
  }

  phase = -1; /* Release internal memory. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,&n, &ddum, ia, ja, &idum, &nrhs,iparm, &msglvl, &ddum, &ddum, &error);

  return 1;
}

#endif

#ifdef USE_PARDISO
template class pardiso_solver<double>;
#endif

#ifdef USE_LAPACK
template class eigen_solver<complex<double> >;
#endif
