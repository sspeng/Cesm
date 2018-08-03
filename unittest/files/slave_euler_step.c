#include <slave.h>
#include <dma.h>
#define NP 4
#define NLEV 30
#define UC 3         // a unit that divides nete by column direction
#define UR 2         // a unit that divides qsize by row direcion
#define NC 8
#define NR 8
#define stripe_qdp (NLEV*NP*NP)

#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0":"=r"(var))

#define athread_get_(src, dst, size_block, n_reply) {\
	get_reply = 0;\
	athread_get(PE_MODE, src, dst, size_block, &get_reply, 0, 0, 0);\
	while(get_reply!=n_reply);\
	asm volatile("memb\n\t");\
}

#define athread_put_(src, dst, size_block, n_reply) {\
	put_reply = 0;\
	athread_put(PE_MODE, src, dst, size_block, &put_reply, 0, 0);\
	while(put_reply!=n_reply);\
	asm volatile("memb\n\t");\
}

#define regcomm_put(RC, src_buf, dst_id) {\
  int i, j;\
  doublev4 *src_p = (doublev4*)(src_buf);\
  register doublev4 d0;\
  for (i = 30; i < NLEV; i++) {\
    for (j = 0; j < NP; j++) {\
      simd_load(d0, src_p);\
      REG_PUT##RC(d0, dst_id);\
      src_p ++;\
    }\
    athread_syn(ROW_SCOPE, 0xff);\
  }\
}

#define regcomm_get(RC, dst_buf) {\
	int i, j;\
	doublev4 *dst_p = (doublev4 *)(dst_buf);\
	register doublev4 d0;\
	for(i = 0; i < NLEV; i++) {\
    for (j = 0; j < NP; j++) {\
      REG_GET##RC(d0);\
      simd_store(d0, dst_p);\
      dst_p ++;\
    }\
    athread_syn(ROW_SCOPE, 0xff);\
	}\
}

typedef struct {
  long qdp_s_ptr, qdp_leap_ptr,dp_s_ptr, dp_leap_ptr, divdp_proj_s_ptr, divdp_proj_leap_ptr, qdp_test_ptr;
  double dt;
  int  nets, nete, np1_qdp, n0_qdp, DSSopt, rhs_multiplier, qsize;
} param_t;

void slave_euler_step_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  get_row_id(rid);
  get_col_id(cid);

  long qdp_s_ptr = param_s->qdp_s_ptr;
  long qdp_leap_ptr = param_s->qdp_leap_ptr;
  long dp_s_ptr = param_s->dp_s_ptr;
  long dp_leap_ptr = param_s->dp_leap_ptr;
  long divdp_proj_s_ptr = param_s->divdp_proj_s_ptr;
  long divdp_proj_leap_ptr = param_s->divdp_proj_leap_ptr;
  long qdp_test_ptr = param_s->qdp_test_ptr;
  double dt = param_s->dt;
  int nets = param_s->nets;
  int nete = param_s->nete;
  int np1_qdp = param_s->np1_qdp;
  int n0_qdp = param_s->n0_qdp;
  int DSSopt = param_s->DSSopt;
  int rhs_multiplier = param_s->rhs_multiplier;
  int qsize = param_s->qsize;

  int block = UR*NLEV*NP*NP;
  double Qdp[UC*block];
  double dp[NLEV*NP*NP];
  int slice_qdp = (int)((qdp_leap_ptr - qdp_s_ptr)/sizeof(double));
  int slice_dp = (int)((dp_leap_ptr - dp_s_ptr)/sizeof(double));
  double *src_np1_qdp = (double *)(qdp_s_ptr) + (np1_qdp - 1)*qsize*stripe_qdp;
  double *src_n0_qdp = (double *)(qdp_s_ptr) + (n0_qdp - 1)*qsize*stripe_qdp;
  double *dst_qdp = (double *)(Qdp);
  double *src_dp = (double *)(dp_s_ptr);
  double *dst_dp = (double *)(dp);

  double *src_qdp_ptr, *dst_qdp_ptr, *src_dp_ptr, *dst_dp_ptr;

  //int rid = id % NR;
  //int cid = id / NC;
  int loop_c = (qsize + UC*NC - 1)/(UC*NC);
  int loop_r = ((nete - nets) + UR*NR - 1)/(UR*NR);
  int c, r, i, j, k, q, ie, cbeg, cend, rbeg, rend, cn, rn;

  if (id == 88)
    printf("%d\n", block);
  if (rid == 0)
    printf("%d", id);
  for (r = 0; r < loop_r; r++) {
    rbeg = r*NR*UR + rid*UR;
    rend = r*NR*UR + (rid + 1)*UR;
    rend = rend < (nete - nets) ? rend : (nete - nets);
    rn = rend - rbeg;
    for (ie = 0; ie < rn; ie++) {
      if (cid == 0) {
        src_dp_ptr = src_dp + (rbeg + ie)*slice_dp;
        dst_dp_ptr = dst_dp;
        athread_get_(src_qdp_ptr, dst_qdp_ptr, (NLEV*NP*NP*sizeof(double)), 1);
        int rank;
        for (rank = 1; rank < 8; rank++) {
          regcomm_put(R, dst_dp_ptr, rank);
          athread_syn(COL_SCOPE, 0xff);
        }
      } else {
        regcomm_get(R, dst_dp_ptr);
      }
      for (c = 0; c < loop_c; c++) {
        cbeg = c*NC*UC + cid*UC;
        cend = c*NC*UC + (cid + 1)*UC;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        src_qdp_ptr = src_n0_qdp + (rbeg + ie)*slice_qdp + cbeg*stripe_qdp;
        dst_qdp_ptr = dst_qdp + ie*block;
        if (cn > 0) {
          athread_get_(src_qdp_ptr, dst_qdp_ptr, (block*sizeof(double)), 1);
          athread_put_(dst_qdp_ptr, src_qdp_ptr, (rn*stripe_qdp*sizeof(double)), 1);
        }
      }
    }
  }

  if (id ==0) {
    printf("%lf,%d,%d,%d,%d,%d,%d,%d\n", dt, nets, nete, np1_qdp, n0_qdp, DSSopt, rhs_multiplier, qsize);
  }

  if (id == 66) {
    athread_get_(src_n0_qdp, dst_qdp, (block*sizeof(double)), 1);
    int i;
    for(i = 0; i < block; i++) {
      printf("%lf\n", Qdp[i]);
    }
    athread_put_(dst_qdp, src_n0_qdp, (block*sizeof(double)), 1);
  }
}
