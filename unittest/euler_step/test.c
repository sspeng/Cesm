#include <stdio.h>
#include <stdlib.h>
#define max(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define min(val1, val2) ((val1) < (val2) ? (val1) : (val2))
#define  abs(value) (__asm__ __volatile__ ("fcpys $31 , %0 , %0 " :  : "r"(va) ))
#define sum(sum, array, len) {   \
  int i;   \
  for (i = 0; i < len; i++)     \
    sum = sum + array[i];  \
}

int main() {
  //double a, b, max, min;
  //a = 7.0f; b = 9.0f;
  //max = max(a, b);
  //min = min(a, b);
  //printf("max:%lf, min:%lf\n", max, min);
  //return 0;
  //int a[10][10];
  //int i, j, k;
  //for (i = 0; i < 10; i++) {
  //  for (j = 0; j < 10; j++) {
  //    a[i][j] = 10*i + j;
  //  }
  //}
  //for (i = 0; i < 10; i++) {
  //  for (j = 0; j < 10; j++) {
  //    if (j > 3 && j < 5)
  //      continue;
  //    printf("%d, ", a[i][j]);
  //    if (j > 6)
  //      break;
  //  }
  //  printf("\n");
  //}
  //return 0;
  int sum_a = 0, sum_b = 0, sum_c = 0;
  int a[16], b[16], c[16];
  int test = -6;
  int i = 0;
  for (i = 0; i < 16; i++) {
    a[i] = 1;
    b[i] = 1;
  }
  for (i = 0; i < 16; i++) {
    sum_b = sum_b + b[i];
  }
  test = abs(test);
  sum(sum_a, a, 16);
  printf("sum_a:%d, sum_b:%d, test_abs:%d\n", sum_a, sum_b, test);
  return 0;
}
