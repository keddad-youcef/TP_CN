/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"
#include "../include/blaslapack_headers.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  for (int j=0;j<(*la);j++){
    int k = j * (*lab);
    if (*kv>=0){
      for (int i=0;i< *kv;i++){
	      AB[k + i] = 0.0;
      }
    }
    // Set up the tridiagonals values to -1 & 2
    AB[k + *kv] = -1.0;
    AB[k + *kv + 1] = 2.0;
    AB[k + *kv + 2] = -1.0;
  }
  // First & Last element = 0
  AB[0] = 0.0;
  if (*kv == 1) {
    AB[1]=0;
  }
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  for (int j=0;j<(*la);j++){
    int k = j * (*lab);
    if (*kv>=0){
      for (int i=0;i< *kv;i++){
	      AB[k + i] = 0.0;
      }
    }
    // Set up the tridiagonals values to -1 & 2
    AB[k + *kv]=0.0;
    AB[k + *kv + 1] = 1.0;
    AB[k + *kv + 2] = 0.0;
  }
  // First & Last element = 0
  AB[1] = 0.0;
  AB[(*lab)*(*la)-1] = 0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // Set up the First & Last element of the Vector to -5 et +5
  RHS[0] = *BC0;
  RHS[*la-1] = *BC1;

  // the remaining elements are equal to 0
  for (int i=1;i<(*la-1);i++) {
    RHS[i] = 0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  const double DELTA_T = *BC1 - *BC0;   // DELTA TEMPERATURE

  for (int i=0;i<(*la);i++) {
    EX_SOL[i] = *BC0 + X[i] * DELTA_T;      // Analytical Sol Equation
  }
}  

void set_grid_points_1D(double* x, int* la){
  double h = 1.0 /((*la) + 1);           // le pas h
  for (int i = 0; i < (*la); i++) {
    x[i] = (i + 1) * h;
  }
}

// LU FACTORISATION FUNCTION
void LU_Facto(double* AB, int *lab, int *la, int *kv){
  int i, j, k, diags = 3;

  if (*kv>=0){
    diags = 4;
    for (i=0;i< *kv;i++){
      AB[i]=0.0;        // Set the first element to zero 
    }
  }
  AB[*kv+2]/=AB[*kv+1];
    

  for (j=1;j<(*la);j++){
    k = j*(*lab);
    if (*kv>=0){
      for (i=0;i< *kv;i++){
        AB[k+i]=0.0;      // Set First Column to zeros
      }
    }

      AB[k+ *kv+1]-=AB[k+ *kv]*AB[(k-diags)+ *kv+2];
      AB[k+ *kv+2]/=AB[k+ *kv+1];
  }

}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab){
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}



void eig_poisson1D(double* eigval, int *la){
  int i;
  double val;
  for (i = 0; i < *la; i++) {
    val = (1.0 * i + 1.0) * M_PI_2 * (1.0 / (*la + 1));
    eigval[i] = sin(val);
    eigval[i] = 4 * eigval[i] * eigval[i];
  }
}

double eigmax_poisson1D(int *la){
  return 4 * pow(sin(M_PI * (*la) * (1.0 /((*la) + 1)) / 2), 2);
}

double eigmin_poisson1D(int *la){
  return 4 * pow(sin(M_PI * (1.0 /((*la) + 1)) / 2), 2);
}

double richardson_alpha_opt(int *la){
  
  return 2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  
  cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);    // r = A * X
  cblas_dscal(*la, -1.0, resvec, 1);              // r = -A * x
  cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);       // r = b -A * X


  double releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);      // erreur residuelle

  while (releres > (*tol) && *nbite < *maxit)
  {
    cblas_daxpy(*la, *alpha_rich, resvec, 1, X, 1);           // iteration de Richardson


    // Calcul de residu r
    cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);
    cblas_dscal(*la, -1.0, resvec, 1);
    cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);


    //Erreur residuelle
    releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);

    // increment iteration
    *nbite += 1;
  }
  printf("\n\n\n iter richardson = %d\n\n\n", *nbite);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for (int j=0;j<(*la);j++){

    int k = j * (*lab);
    MB[k + *kv] = 1 / AB[k + *kv];    // M = inv(Diag(AB))
    
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  double *ONES;
  int NRHS = 1;
  int *info;
  double *M = malloc(sizeof(double) * (*la) * (*la));


  // M = D - E
  for (int i = 0; i < (*la) * (*la); i++) { M[i] = 0.0; }

  for (int i = 0; i < (*la); i++) { M[i * (*la) + i] = AB[(*lab) * i + 1]; }    // extraction du diagonale

  for (int i = 0; i < (*la) - 1; i++) {
    M[(i + 1) * (*la) + i] = AB[(*lab) * i + 2];      // extraction de la sous diagonale
  }


  // Calcul de l'inverse de M
  int *ipiv = malloc(sizeof(int) * (*la));
  LAPACKE_dgetrf(CblasRowMajor, (*la), (*la), M, (*la), ipiv);
  LAPACKE_dgetri(CblasRowMajor, (*la), M, (*la), ipiv);
  free(ipiv);
  

  
  // Copie dans MB
  cblas_dcopy((*la)*(*la),M, 1, MB, 1);
  

  free(M);
  
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite, int resol_type){

  
  // Calcul de residu r
  cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);
  cblas_dscal(*la, -1.0, resvec, 1);
  cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);


  //Erreur residuelle
  double releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);
  


  while (releres > (*tol) && *nbite < *maxit)
  {
    if (resol_type == 1)
    {
      cblas_dgemv(CblasColMajor, CblasConjNoTrans, *la, *la, 1.0, MB, *la, resvec, 1, 1.0, X, 1);   // Xk = Xk + M(gauss) * r
    }
    else {cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, 1.0, MB, *lab, resvec, 1, 1.0, X, 1);}   // Xk = Xk + M(jacobi) * r
    
    
    
    // Calcul de residu r
    cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, 1.0, AB, *lab, X, 1, 0.0, resvec, 1);
    cblas_dscal(*la, -1.0, resvec, 1);
    cblas_daxpy(*la, 1.0, RHS, 1, resvec, 1);

    //Erreur residuelle
    releres = cblas_dnrm2(*la, resvec, 1) / cblas_dnrm2(*la, RHS, 1);

    // increment iteration
    *nbite += 1;
 
  }
  printf("\n\n\n iter = %d\n\n\n", *nbite);
}

