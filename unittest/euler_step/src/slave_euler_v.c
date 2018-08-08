#include <slave.h>
#include "dma_macros.h"
//#include "ldm_alloc.h"
#define NP 4
#define NLEV 30
#define UC 4         // a unit that divides nete by column direction
#define UR 3         // a unit that divides qsize by row direcion
#define NC 4
#define NR 16
#define stripe_qdp (NLEV*NP*NP)
#define qstep_Qten (NLEV*NP*NP)       // stripe in q axis of Qtens_biharmonic array
#define qstep_qmax  NLEV
#define block  (UC*NLEV*NP*NP)


#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0":"=r"(var))
#define max(val1, val2) ((val1) > (val2) ? (val1) : (val2))
#define min(val1, val2) ((val1) < (val2) ? (val1) : (val2))


typedef struct {
  double *qdp_s_ptr, *qdp_leap_ptr, *divdp_proj, *vn0, *dp_temp, *gradQ_temp;
  double dt;
  int  nets, nete, rhs_multiplier, qsize;
} param_t;

void slave_euler_v_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
	dma_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

}
