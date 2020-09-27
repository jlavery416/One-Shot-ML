#include <stdio.h>
#include <stdlib.h>

double * invert(double *x, int N) {
  int c,i,j,k;
  double *aug = (double *) malloc(N * N * sizeof(double));
  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++) {
      if(i == j) {
	*(aug + (i * N + j)) = 1.0;
      } else {
	*(aug + (i * N + j)) = 0.0;
      }
    }
  }
  for(c = 0; c < N; c++) {
    double divisor = *(x + (c * N + c)); 
    for(j = 0; j < N; j++) {
      *(x + (c * N + j)) /= divisor;
      *(aug + (c * N + j)) /= divisor;
    }
    for(i = 0; i < N; i++) {
      if(i == c) {
	continue;
      }
      double scalar = (-1.0) * (*(x + (i * N + c)));
      if(scalar != 0.0) {
	for(k = 0; k < N; k++) { 
	  *(x + (i * N + k)) += *(x + (c * N + k)) * scalar;
	  *(aug + (i * N + k)) += *(aug + (c * N + k))* scalar;
	}
      }
    }
  }
  return aug;
}

double * multi(double *x, double *y, int m, int n, int p) {
  int i, j, k;
  double *product = (double *) malloc(m * p * sizeof(double));
  for(i = 0; i < m; i++) {
    for(j = 0; j < p; j++) {
      *(product + (i * p + j)) = 0.0;
    }
  }
  for(i = 0; i < m; i++) {
    for(j = 0; j < p; j++) {
      for(k = 0; k < n; k++) {
	*(product + (i * p + j)) += (*(x + (i * n + k))) * (*(y + (k * p + j)));
      }
    }
  }
  return product;
}

double * trans(double *x, int m, int n) {
  int i, j;
  double *xtran = (double *) malloc(n * m * sizeof(double));
  for(i = 0; i < n; i++) {
    for(j = 0; j < m; j++) {
      *(xtran + (i * m + j)) = 0.0;
    }
  }

  for(i = 0; i < n; i++) {
    for(j = 0; j < m; j++) {
      *(xtran + (i * m + j)) = *(x + (j * n + i));
    }
  }
  return xtran;
}

int main(int argc, char *argv[]) {

  //------SET-UP--------------------//
  int K, N, M, i, j;

  FILE *trainptr = fopen(argv[1], "r");                      
  FILE *testptr = fopen(argv[2], "r");

  fscanf(trainptr,"%d", &K);                            
  fscanf(trainptr,"%d", &N);
  fscanf(testptr, "%d", &M);

  double X[N][K+1];
  double Y[N][1];
  double Alpha[M][K+1];
  double *xptr = &X[0][0];
  double *yptr = &Y[0][0];
  double *alphaPtr = &Alpha[0][0];

  //----INITIALIZE ARRAYS------------//    

  for (i = 0; i < N; i++) {
    for (j = 0; j < (K+1); j++) {
      if(j == 0) {
	X[i][j] = 1.0;
      } else {
	X[i][j] = 0.0;
      }
    }
  }                             
  

  for (i = 0; i < N; i++) {
    for (j = 0; j < 1; j++) {
      Y[i][j] = 0.0;
    }
  }

  for(i = 0; i < M; i++) {
    for(j = 0; j < (K+1); j++) {
      if(j == 0) {
	Alpha[i][j] = 1.0;
      } else {
	Alpha[i][j] = 0.0;
      }
    }
  }


  //------READ IN DATA--------------//

  for(i = 0; i < N; i++) {                                    
    for(j = 0; j <= K; j++) {
      if(j == 0) {
        fscanf(trainptr,"%lf,", &Y[i][j]);
      } else {
	fscanf(trainptr,"%lf,", &X[i][j]);
      }
    }
  }



  for(i = 0; i < M; i++) {
    for(j = 0; j < K; j++) {
      fscanf(testptr, "%lf,", &Alpha[i][j+1]);
    }
  }

  //------PERFORM COMPUTATIONS-----//

  double *xT = trans(xptr, N, (K+1));
  double *A = multi(xT, xptr, (K+1), N, (K+1));
  double *in = invert(A, (K+1));
  double *B = multi(xT, yptr, (K+1), N, 1);
  double *w = multi(in, B, (K+1), (K+1), 1);
  double *res = multi(alphaPtr, w, M, (K+1), 1);

  //-------PRODUCE OUTPUT----------//

  for(i = 0; i < M; i++) {
    for(j = 0; j < 1; j++) {
      printf("%0.0lf\n", *(res + i + j));
    }
  }
  //------FREE UP MEMORY----------// 
  free(xT);
  free(A);
  free(in);
  free(B);
  free(w);
  free(res);
  fclose(trainptr);
  fclose(testptr);
  return 0;
}
