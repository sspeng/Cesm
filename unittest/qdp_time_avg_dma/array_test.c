#include <stdio.h>

typedef struct {
  int a1, a2, a3, a4, a5;
  double a[10][10][10][10];
} data;

int main() {
  data st[10];
  double *b = (double *)(st[0].a);
  long c = (long)(b);
  double *d = (double *)(c);
  int i, j, k, m, dd;

  for (dd = 0; dd < 10; dd++) {
    for (i = 0; i < 10; i++) {
      for (j = 0; j < 10; j++) {
        for (k = 0; k < 10; k++) {
          for (m = 0; m < 10; m++) {
            st[dd].a[i][j][k][m] = dd * 10000 + i * 1000 + j * 100 + k * 10 + m;
          }
        }
      }
    }
  }
  int count = 0;
  double *b1, b2, b3, d1, d2, d3;
  long c1, c2, c3, c_delta, cs, ce;
  int n_d;
  cs = (long)(st[0].a);
  ce = (long)(st[1].a);
  c1 = (long)(st[0].a);
  c2 = (long)(st[2].a);
  c_delta = ce - cs;
  n_d = (int)(c_delta / sizeof(double));
  c3 = c1 + 2 * n_d * sizeof(double);
  printf("c1:%ld, c2:%ld, c_delta:%ld\n", c1, c2, c3);
  //for (i = 0; i < 10; i++) {
  //  for (j = 0; j < 10; j++) {
  //    for (k = 0; k < 10; k++) {
  //      for (m = 0; m < 10; m++) {
  //        printf("b:%lf, *B:%lf \n", b[i * 1000 + j * 100 + k * 10 + m], *(d + i * 1000 + j * 100 + k * 10 + m));
  //      }
  //      //printf("\n");
  //    }
  //  }
  //}

  return 0;
}
