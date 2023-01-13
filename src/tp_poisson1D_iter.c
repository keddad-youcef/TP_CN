/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include "../include/blaslapack_headers.h"
#include <time.h>

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;
  double *MB;
  
  double temp, relres;

  double opt_alpha;

  /* Size of the problem */
  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double)); 
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  /********************************************/
  /* Solution (Richardson with optimal alpha) */
  printf("*********Solution Richardson***********\n\n\n\n\n");
  /* Computation of optimum alpha */
  opt_alpha = richardson_alpha_opt(&la);
  printf("Optimal alpha for simple Richardson iteration is : %lf",opt_alpha); 

  /* Solve */
  double tol=1e-3;
  int maxit=10000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  clock_t richardson_begin = clock();
  richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  clock_t richardson_end = clock();


  // Computes execution time
  double richardson_delta = (double) (richardson_end - richardson_begin) / CLOCKS_PER_SEC;
  printf("Execution time RICHARDSON = %f sec\n", richardson_delta);


  /* Write solution */
  write_vec(SOL, &la, "SOL_RICHARDSON.dat");

  /* Write convergence history */
  write_vec(resvec, &nbite, "RESVEC_RICHARDSON.dat");


  /* Relative forward error */

  // Calcul de || X ||
  temp = cblas_dnrm2(la, EX_SOL, 1);
  // Calcul de || X - X'||
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_dnrm2(la, EX_SOL, 1);
  // || X - X'|| / || X ||
  relres = relres / temp;

  printf("\nThe relative forward error for RICHARDSON is relres = %e\n\n",relres);
  
  
  free(SOL);
  free(resvec);
  nbite = 0;

  /* Richardson General Tridiag */

  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  kv = 1;
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(la)*(la));


  // Jakobi Call
  printf("*********Solution JAKOBI***********\n\n\n\n\n");
  SOL=(double *) calloc(la, sizeof(double));
  resvec=(double *) calloc(maxit, sizeof(double));


  // Extract Jaccobi MATRIX
  clock_t jacobi_begin = clock();
  extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  
  /* Solve with General Richardson */
  richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite, 0);
  clock_t jacobi_end = clock();

  // Computes execution time
  double jacobi_delta = (double) (jacobi_end - jacobi_begin) / CLOCKS_PER_SEC;
  printf("Execution time JAKOBI = %f sec\n", jacobi_delta);


  /* Write solution */
  write_vec(SOL, &la, "SOL_JAKOBI.dat");

  /* Write convergence history */
  write_vec(resvec, &nbite, "RESVEC_JAKOBI.dat");


  /* Relative forward error */
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  // Calcul de || X ||
  temp = cblas_dnrm2(la, EX_SOL, 1);

  // Calcul de || X - X'||
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_dnrm2(la, EX_SOL, 1);

  // || X - X'|| / || X ||
  relres = relres / temp;

  printf("\nThe relative forward error for JACOBI is relres = %e\n\n",relres);


  free(SOL);
  free(MB);
  free(resvec);
  nbite = 0;
  



  printf("*********Solution GAUSS SEIDEL***********\n\n\n\n\n");
  // Gauss SEIDEL Call
  SOL=(double *) calloc(la, sizeof(double));
  MB = (double *) malloc(sizeof(double)*(la)*(la));
  resvec=(double *) calloc(maxit, sizeof(double));

  // EXTRACT GAUSS SEIDEL MATRIX
  clock_t gs_begin = clock();
  extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  
  /* Solve with General Richardson */
  richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite, 1);
  clock_t gs_end = clock();

  

  // Computes execution time
  double gs_delta = (double) (gs_end - gs_begin) / CLOCKS_PER_SEC;
  printf("Execution time GAUSS SEIDEL = %f sec\n", gs_delta);


  /* Write solution */
  write_vec(SOL, &la, "SOL_GAUSS_SEIDEL.dat");

  /* Write convergence history */
  write_vec(resvec, &nbite, "RESVEC_GAUSS_SEIDEL.dat");


  /* Relative forward error */
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  // Calcul de || X ||
  temp = cblas_dnrm2(la, EX_SOL, 1);
  // Calcul de || X - X'||
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_dnrm2(la, EX_SOL, 1);
  // || X - X'|| / || X ||
  relres = relres / temp;

  printf("\nThe relative forward error for GAUSS SEIDEL is relres = %e\n",relres);

  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End -----------\n");
}
