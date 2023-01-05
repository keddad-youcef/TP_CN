/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include "../include/blaslapack_headers.h"
#include <time.h>


int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X, *RHS2;
  double **AAB;
  double *AB, *AB_LU;

  double temp, relres;

  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  RHS2=(double *) malloc(sizeof(double)*la);


  
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);
  AB_LU = (double *) malloc(sizeof(double)*lab*la);   // to store AB after LU factorisation

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_GB_operator_colMajor_poisson1D(AB_LU, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  


  // CBLAS DGBMV
  cblas_dgbmv(CblasColMajor, CblasConjNoTrans, la, la, kl, ku, 1.0, AB+1, lab, EX_SOL, 1, 0.0, RHS2, 1);

  // Validation
  cblas_daxpy(la, -1, RHS2, 1, RHS, 1);
  double norm_mv = cblas_dnrm2(la, RHS, 1);

  
  
  //******** A*X = b PROBLEM **********//

  printf("Solution with LAPACK\n");
  /* LU Factorization */
  info=0;
  ipiv = (int *) calloc(la, sizeof(int));
  


  //******LAPACKE_dgbtrs**********// 
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  clock_t dgbtrs_begin = clock();
  info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, la, la, kl, ku, AB, lab, ipiv);   // LAPACKE CALL
  info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
  clock_t dgbtrs_end = clock();
  write_vec(RHS, &la, "dgbtrs_sol.dat");
  
  
  // Computes execution time
  double dgbtrs_delta = (double) (dgbtrs_end - dgbtrs_begin) / CLOCKS_PER_SEC;
  printf("Execution time DGBTRS = %f sec\n", dgbtrs_delta);

  

  //******LAPACKE_dgbsv**********// 

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  clock_t dgbsv_begin = clock();
  info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
  write_vec(RHS, &la, "dgbsv_sol.dat");
  clock_t dgbsv_end = clock();

  // Computes execution time
  double dgbsv_delta = (double) (dgbsv_end - dgbsv_begin) / CLOCKS_PER_SEC;
  printf("Execution time DGBSV = %f sec\n", dgbsv_delta);


  


  
  //*********** LU FACTORISATION ************//
  

  // FUNCTION CALL
  LU_Facto(AB_LU, &lab, &la, &kv);     // Factorisation Implementation                         
  write_GB_operator_colMajor_poisson1D(AB_LU, &lab, &la, "AB_LU_Facto.dat");


  // LAPACK dgbtrf
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  LAPACKE_dgbtrf(LAPACK_COL_MAJOR, la, la, kl, ku, AB, lab, ipiv);
  write_GB_operator_colMajor_poisson1D(AB_LU, &lab, &la, "LAPACKE_dgbtrf.dat");


  // Validation
  cblas_daxpy(la, -1, AB, 1, AB_LU, 1);
  double norm_lu = cblas_dnrm2(la, AB_LU, 1);

  // Print norms
  printf("\n\n dgbmv() est valide et l'erreur = %e\n", norm_mv);
  printf("\n\n LU_FACTO() est valide et l'erreur = %e\n", norm_lu);



  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */

  // Calcul de || X ||
  temp = cblas_dnrm2(la, EX_SOL, 1);

  // Calcul de || X - X'||
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_dnrm2(la, EX_SOL, 1);

  // || X - X'|| / || X ||
  relres = relres / temp;
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(RHS2);
  free(EX_SOL);
  free(X);
  free(AB);
  free(AB_LU);
  printf("\n\n--------- End -----------\n");
}
