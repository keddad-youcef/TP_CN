/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

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
    int i, j, k, k1 = 3;


      if (*kv>=0){
        k1 = 4;
        for (i=0;i< *kv;i++){
            AB[i]=0.0;
        }
      }
      AB[*kv+2]/=AB[*kv+1];
    

    for (j=1;j<(*la);j++){
      k = j*(*lab);
      if (*kv>=0){
        for (i=0;i< *kv;i++){
            AB[k+i]=0.0;
        }
      }

      AB[k+ *kv+1]-=AB[k+ *kv]*AB[(k-k1)+ *kv+2];
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
