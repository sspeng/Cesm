#include <stdio.h>

int main() {
  double a[10][10][10][10];
  double *b = (double *)(a);
  int i, j, k, m;
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 10; j++) {
      for (k = 0; k < 10; k++) {
        for (m = 0; m < 10; m++) {
          a[i][j][k][m] = i * 1000 + j * 100 + k * 10 + m;
        }
      }
    }
  }
  int count = 0;
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 10; j++) {
      for (k = 0; k < 10; k++) {
        for (m = 0; m < 10; m++) {
          printf("b:%lf, *B:%lf \n", b[i * 1000 + j * 100 + k * 10 + m], *(b + i * 1000 + j * 100 + k * 10 + m));
        }
        //printf("\n");
      }
    }
  }

  return 0;
}
