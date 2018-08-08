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
  long qdp_s_ptr, qdp_leap_ptr, dp_s_ptr, divdp_proj_s_ptr \
	    , Qtens_biharmonic, qmax, qmin;
  double dt;
  int  nets, nete, np1_qdp, n0_qdp, DSSopt, rhs_multiplier, qsize;
} param_t;

void slave_euler_step_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
  //get_row_id(rid);
  //get_col_id(cid);
	dma_init();
	//ldm_alloc_init();

	/*
  long qdp_s_ptr = param_s->qdp_s_ptr;
  long qdp_leap_ptr = param_s->qdp_leap_ptr;
  long dp_s_ptr = param_s->dp_s_ptr;
  long dp_leap_ptr = param_s->dp_leap_ptr;
  long divdp_proj_s_ptr = param_s->divdp_proj_s_ptr;
  long divdp_proj_leap_ptr = param_s->divdp_proj_leap_ptr;
	long Qtens_biharmonic_ptr = param_s->Qtens_biharmonic;
  double dt = param_s->dt;
  int nets = param_s->nets;
  int nete = param_s->nete;
  int np1_qdp = param_s->np1_qdp;
  int n0_qdp = param_s->n0_qdp;
  int DSSopt = param_s->DSSopt;
  int rhs_multiplier = param_s->rhs_multiplier;
  int qsize = param_s->qsize;
	*/

	param_t param_d;
	pe_get(param_s, &param_d, sizeof(param_t));
	dma_syn();
	long qdp_s_ptr = param_d.qdp_s_ptr;
	long qdp_leap_ptr = param_d.qdp_leap_ptr;
	long dp_s_ptr = param_d.dp_s_ptr;
	long divdp_proj_s_ptr = param_d.divdp_proj_s_ptr;
	long Qtens_biharmonic_ptr = param_d.Qtens_biharmonic;
	long qmax_ptr = param_d.qmax;
	long qmin_ptr = param_d.qmin;
	double dt = param_d.dt;
	int nets = param_d.nets;
	int nete = param_d.nete;
	int np1_qdp = param_d.np1_qdp;
	int n0_qdp = param_d.n0_qdp;
	int DSSopt = param_d.DSSopt;
	int rhs_multiplier = param_d.rhs_multiplier;
	int qsize = param_d.qsize;

  int istep_Qten = qsize*NLEV*NP*NP; // stripe in ie axis of Qtens_biharmonic array
	int istep_qmax = qsize*NLEV;
  double Qdp[block];
  double dp[NLEV*NP*NP];
	double divdp_proj[NLEV*NP*NP];
	double dp_l[NLEV*NP*NP];
	double Qtens_biharmonic_l[block];
	double qmax_val[NLEV];
	double qmin_val[NLEV];
	double qmax[UC*NLEV];
	double qmin[UC*NLEV];
	//double *qmax, *qmin;
	//ldm_alloc(qmax, (nete - nets + 1)*qsize*NLEV);
	//ldm_alloc(qmin, (nete - nets + 1)*qsize*NLEV);

  int slice_qdp = (int)((qdp_leap_ptr - qdp_s_ptr)/sizeof(double));  //stripe in ie axis of elem.state.Qdp
  double *src_np1_qdp = (double *)(qdp_s_ptr) + (np1_qdp - 1)*qsize*stripe_qdp;
  double *src_n0_qdp = (double *)(qdp_s_ptr) + (n0_qdp - 1)*qsize*stripe_qdp;
  double *src_dp = (double *)(dp_s_ptr);
	double *src_divdp_proj = (double *)(divdp_proj_s_ptr);
	double *src_Qtens_bi = (double *)(Qtens_biharmonic_ptr);
	double *src_qmax = (double *)(qmax_ptr);
	double *src_qmin = (double *)(qmin_ptr);

  double *src_qdp_ptr, *src_dp_ptr,  *src_divdp_proj_ptr, *src_Qtens_bi_ptr  \
			, *src_qmax_ptr, *src_qmin_ptr;

  rid = id / NC;
  cid = id % NC;
  int loop_r = ((nete - nets + 1) + UR*NR - 1)/(UR*NR);
	int loop_c = (qsize + UC*NC - 1)/(UC*NC);
  int c, r, i, j, k, q, ie, cbeg, cend, rbeg, rend, cn, rn, pos_dp, pos_qdp, pos_Qtens_bi, pos_qmax;

  for (r = 0; r < loop_r; r++) {
    rbeg = r*NR*UR + rid*UR;
    rend = r*NR*UR + (rid + 1)*UR;
    rend = rend < (nete - nets + 1) ? rend : (nete - nets + 1);
    rn = rend - rbeg;
		rn = rn < 0 ? 0 : rn;
    for (ie = 0; ie < rn; ie++) {
			src_dp_ptr = src_dp + (rbeg + ie)*slice_qdp;
			src_divdp_proj_ptr = src_divdp_proj + (rbeg + ie)*slice_qdp;
      pe_get(src_dp_ptr, dp, (NLEV*NP*NP*sizeof(double)));
			pe_get(src_divdp_proj_ptr, divdp_proj, (NLEV*NP*NP*sizeof(double)));
			dma_syn();
			for (k = 0; k < NLEV; k++) {
				for (j = 0; j < NP; j++) {
					for (i = 0; i < NP; i++) {
						pos_dp = k*NP*NP + j*NP + i;
						dp_l[pos_dp] = dp[pos_dp] - rhs_multiplier*dt*divdp_proj[pos_dp];
					}
				}
			}
      for (c = 0; c < loop_c; c++) {
        cbeg = c*NC*UC + cid*UC;
        cend = c*NC*UC + (cid + 1)*UC;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        src_qdp_ptr = src_n0_qdp + (rbeg + ie)*slice_qdp + cbeg*stripe_qdp;
				src_Qtens_bi_ptr = src_Qtens_bi + (rbeg + ie)*istep_Qten + cbeg*qstep_Qten;
				src_qmax_ptr = src_qmax + (rbeg + ie)*istep_qmax + cbeg*qstep_qmax;
				src_qmin_ptr = src_qmin + (rbeg + ie)*istep_qmax + cbeg*qstep_qmax;
        if (cn > 0) {
          pe_get(src_qdp_ptr, Qdp, (block*sizeof(double)));
					pe_get(src_Qtens_bi_ptr, Qtens_biharmonic_l, (block*sizeof(double)));
					dma_syn();
					for (q = 0; q < cn; q++) {
						for (k = 0; k < NLEV; k++) {
							qmin_val[k] = +1.0e+24;
							qmax_val[k] = -1.0e+24;
							for (j = 0; j < NP; j++) {
								for (i = 0; i < NP; i++) {
									pos_Qtens_bi = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
									pos_dp = k*NP*NP + j*NP + i;
									Qtens_biharmonic_l[pos_Qtens_bi] = Qdp[pos_Qtens_bi]/dp_l[pos_dp];
									qmin_val[k] = min(qmin_val[k], Qtens_biharmonic_l[pos_Qtens_bi]);
									qmax_val[k] = max(qmax_val[k], Qtens_biharmonic_l[pos_Qtens_bi]);
								}
							}
							if (rhs_multiplier == 1) {
								pos_qmax = q*NLEV + k;
								qmin[pos_qmax] = min(qmin[pos_qmax], qmin_val[k]);
								qmin[pos_qmax] = max(qmin[pos_qmax], 0.0);
								qmax[pos_qmax] = max(qmax[pos_qmax], qmax_val[k]);
							} else {
								pos_qmax = q*NLEV + k;
								qmin[pos_qmax] = max(qmin_val[k], 0.0);
								qmax[pos_qmax] = qmax_val[k];
							}
						}
					}
					pe_put(src_Qtens_bi_ptr, Qtens_biharmonic_l, (cn*NLEV*NP*NP*sizeof(double)));
					pe_put(src_qmax_ptr, qmax, cn*NLEV*sizeof(double));
					pe_put(src_qmin_ptr, qmin, cn*NLEV*sizeof(double));
					dma_syn();
        }
      }
    }
  }

	//ldm_dealloc(qmax);
	//ldm_dealloc(qmin);
}
