#include <slave.h>
#include "dma_macros.h"
//#include "ldm_alloc.h"
#define NP    4
#define NLEV  30
#define UC    1         // a unit that divides nete by column direction
#define UR    3         // a unit that divides qsize by row direcion
#define NC    4
#define NR    16
#define stripe_qdp        (NLEV*NP*NP)
#define istep_dp          (NLEV*NP*NP)
#define qstep_Qten        (NLEV*NP*NP)       // stripe in q axis of Qtens_biharmonic array
#define qstep_dp_star     (NLEV*NP*NP)
#define qstep_qmax        NLEV
#define qstep_Qtens       (NLEV*NP*NP)
#define block             (UC*NLEV*NP*NP)
#define block_dp          (NLEV*NP*NP)
#define block_vn0         (NLEV*2*NP*NP)
#define block_gradQ       (2*NP*NP)
#define block_dp_star     (UC*NLEV*NP*NP)
#define block_Dinv        (2*2*NP*NP)
#define block_det         (NP*NP)

typedef struct {
  double *qdp_s_ptr, *qdp_leap_ptr, *divdp_proj, *dp, *vn0, *dp_temp,          \
      *dp_star_temp, *Dvv, *Dinv, *metdet, *rmetdet, *Qtens_temp,              \
      *Qtens_biharmonic, *divdp, *dpdiss_biharmonic, *spheremp;
  double dt, rrearth, nu_p, nu_q;
  int nets, nete, rhs_multiplier, qsize, n0_qdp, np1_qdp, limiter_option       \
      , rhs_viss;
} param_t;

void slave_euler_v_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
	dma_init();

  param_t param_d;
  pe_get(param_s, &param_d, sizeof(param_t));
  dma_syn();

  double *qdp_s_ptr = param_d.qdp_s_ptr;
  double *qdp_leap_ptr = param_d.qdp_leap_ptr;
  double *divdp_proj_ptr = param_d.divdp_proj;
  double *dp_ptr = param_d.dp;
  double *vn0_ptr = param_d.vn0;
  double *dp_temp_ptr = param_d.dp_temp;
  double *dp_star_ptr = param_d.dp_star_temp;
  double *Dvv_ptr = param_d.Dvv;
  double *Dinv_ptr = param_d.Dinv;
  double *metdet_ptr = param_d.metdet;
  double *rmetdet_ptr = param_d.rmetdet;
  double *Qtens_temp_ptr = param_d.Qtens_temp;
  double *Qtens_biharmonic_ptr = param_d.Qtens_biharmonic;
  double *divdp_ptr = param_d.divdp;
  double *dpdiss_biharmonic_ptr = param_d.dpdiss_biharmonic;
  double *spheremp_ptr = param_d.spheremp;
  double dt = param_d.dt;
  double rrearth = param_d.rrearth;
  double nu_p = param_d.nu_p;
  double nu_q = param_d.nu_q;
  int nets = param_d.nets;
  int nete = param_d.nete;
  int rhs_multiplier = param_d.rhs_multiplier;
  int qsize = param_d.qsize;
  int n0_qdp = param_d.n0_qdp;
  int np1_qdp = param_d.np1_qdp;
  int limiter_option = param_d.limiter_option;
  int rhs_viss = param_d.rhs_viss;

  rid = id / NC;
  cid = id % NC;
  int loop_r = ((nete - nets + 1) + UR*NR - 1)/(UR*NR);
  int loop_c = (qsize + UC*NC - 1)/(UC*NC);
  int slice_qdp = (int)(qdp_leap_ptr - qdp_s_ptr);
  int istep_dp_star = qsize*NLEV*NP*NP;
  int istep_Qtens = qsize*NLEV*NP*NP;
  int c, r, i, j, l, k, q, ie, cbeg, cend, rbeg, rend, cn, rn;

  double Qdp[block];
  double dp[block_dp];
  double dp_temp[block_dp];
  double dp_star[block_dp_star];
  double vn0[block_vn0];
  double Vstar[block_vn0];                 // same size as vn0
  double divdp_proj[block_dp];             // same size as dp_tmp
  double gradQ[block_gradQ];
  double Dvv[NP*NP];
  double Dinv[2*2*NP*NP];
  double metdet[NP*NP];
  double rmetdet[NP*NP];
  double Qtens_temp[block];
  double Qtens_biharmonic[block];
  double dpdiss_biharmonic[block_dp];     // same size as dp
  double divdp[block_dp];                 // same size as dp
  double dpdiss[NP*NP];
  double spheremp[NP*NP];
  double xx[NP*NP];
  double cc[NP*NP];
  double qmax[UC*NLEV];
  double qmin[UC*NLEV];

  if (id == 66) {
    int size_tol = sizeof(double)*(block + block_dp + block_dp_star + block_vn0 \
        + block_vn0 + block_dp + block_gradQ + NP*NP + 4*NP*NP + NP*NP + NP*NP  \
        + block + block + block_dp + block_dp + NP*NP + NP*NP + NP*NP + NP*NP   \
        + UC*NLEV + UC*NLEV);
    printf("size_tol:%dk\n", size_tol/1024);
  }

  pe_get(Dvv_ptr, Dvv, NP*NP*sizeof(double));
  dma_syn();

  // local Variables
  double dudx00, dvdy00;
  double gv[2*NP*NP], vvtemp[NP*NP];

  double *src_qdp, *src_dp, *src_dp_temp, *src_vn0, *src_divdp_proj,   \
      *src_Dinv, *src_metdet, *src_rmetdet, *src_dp_star, *src_Qtens_temp,     \
      *src_Qtens_biharmonic, *src_divdp, *src_dpdiss_biharmonic, *src_spheremp;
  double *src_n0_qdp = qdp_s_ptr + (n0_qdp - 1)*qsize*stripe_qdp;
  double *src_np1_qdp = qdp_s_ptr + (np1_qdp - 1)*qsize*stripe_qdp;

  for (r = 0; r < loop_r; r++) {
    rbeg = r*NR*UR + rid*UR;
    rend = r*NR*UR + (rid + 1)*UR;
    rend = rend < (nete - nets + 1) ? rend : (nete - nets + 1);
    rn = rend - rbeg;
    rn = rn < 0 ? 0 : rn;
    for (ie = 0; ie < rn; ie++) {
      src_dp = dp_ptr + (rbeg + ie)*slice_qdp;
      src_dp_temp = dp_temp_ptr + (rbeg + ie)*istep_dp;
      src_vn0 = vn0_ptr + (rbeg + ie)*slice_qdp;
      src_divdp_proj = divdp_proj_ptr + (rbeg + ie)*slice_qdp;
      src_divdp = divdp_ptr + (rbeg + ie)*slice_qdp;
      src_dpdiss_biharmonic = dpdiss_biharmonic_ptr + (rbeg + ie)*slice_qdp;
      src_Dinv = Dinv_ptr + (rbeg + ie)*slice_qdp;
      src_metdet = metdet_ptr + (rbeg + ie)*slice_qdp;
      src_rmetdet = rmetdet_ptr + (rbeg + ie)*slice_qdp;
      src_spheremp = spheremp_ptr + (rbeg + ie)*slice_qdp;
      pe_get(src_dp, dp, block_dp*sizeof(double));
      pe_get(src_dp_temp, dp_temp, block_dp*sizeof(double));
      pe_get(src_vn0, vn0, block_vn0*sizeof(double));
      pe_get(src_divdp_proj, divdp_proj, block_dp*sizeof(double));
      pe_get(src_divdp, divdp, block_dp*sizeof(double));
      pe_get(src_dpdiss_biharmonic, dpdiss_biharmonic, block_dp*sizeof(double));
      pe_get(src_Dinv, Dinv, block_Dinv*sizeof(double));
      pe_get(src_metdet, metdet, block_det*sizeof(double));
      pe_get(src_rmetdet, rmetdet, block_det*sizeof(double));
      pe_get(src_spheremp, spheremp, NP*NP*sizeof(double));
      dma_syn();

      for (k = 0; k < NLEV; k++) {
        for (j = 0; j < NP; j++) {
          for (i = 0; i < NP; i++) {
            int pos_dp = k*NP*NP + j*NP + i;
            int pos_vn0_1 = k*2*NP*NP + j*NP + i;
            int pos_vn0_2 = k*2*NP*NP + NP*NP + j*NP + i;
            dp_temp[pos_dp] = dp[pos_dp] - rhs_multiplier*dt*divdp_proj[pos_dp];
            Vstar[pos_vn0_1] = vn0[pos_vn0_1]/dp_temp[pos_dp];
            Vstar[pos_vn0_2] = vn0[pos_vn0_2]/dp_temp[pos_dp];
          }
        }
      }
      for (c = 0; c < loop_c; c++) {
        cbeg = c*NC*UC + cid*UC;
        cend = c*NC*UC + (cid + 1)*UC;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        if (cn > 0) {
          src_qdp = src_n0_qdp + (rbeg + ie)*slice_qdp + cbeg*stripe_qdp;
          src_Qtens_temp = Qtens_temp_ptr + (rbeg + ie)*istep_Qtens + cbeg*qstep_Qtens;
          src_Qtens_biharmonic = Qtens_biharmonic_ptr + (rbeg + ie)*istep_Qtens   \
              + cbeg*qstep_Qtens;
          src_dp_star = dp_star_ptr + (rbeg + ie)*istep_dp_star + cbeg*qstep_dp_star;
          pe_get(src_qdp, Qdp, block*sizeof(double));
          pe_get(src_Qtens_temp, Qtens_temp, block*sizeof(double));
          pe_get(src_Qtens_biharmonic, Qtens_biharmonic, block*sizeof(double));
          pe_get(src_dp_star, dp_star, block_dp_star*sizeof(double));
          dma_syn();
          for (q = 0; q < cn; q++) {
            for (k = 0; k < NLEV; k++) {
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_qdp = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  int pos_vn0_1 = k*2*NP*NP + j*NP + i;
                  int pos_vn0_2 = k*2*NP*NP + NP*NP + j*NP + i;
                  int pos_gradQ_1 = j*NP + i;
                  int pos_gradQ_2 = NP*NP + j*NP + i;
                  int pos_Dinv_1 = j*NP + i;
                  int pos_Dinv_2 = NP*NP + j*NP + i;
                  int pos_Dinv_3 = 2*NP*NP + j*NP + i;
                  int pos_Dinv_4 = 3*NP*NP + j*NP + i;
                  int pos_gv_1 = j*NP + i;
                  int pos_gv_2 = NP*NP + j*NP + i;
                  int pos_det = j*NP + i;
                  gradQ[pos_gradQ_1] = Vstar[pos_vn0_1]*Qdp[pos_qdp];
                  gradQ[pos_gradQ_2] = Vstar[pos_vn0_2]*Qdp[pos_qdp];
                  gv[pos_gv_1] = metdet[pos_det]*(Dinv[pos_Dinv_1]*gradQ[pos_gradQ_1] \
                      + Dinv[pos_Dinv_3]*gradQ[pos_gradQ_2]);
                  gv[pos_gv_2] = metdet[pos_det]*(Dinv[pos_Dinv_2]*gradQ[pos_gradQ_1] \
                      + Dinv[pos_Dinv_4]*gradQ[pos_gradQ_2]);
                }
              }
              for (j = 0; j < NP; j++) {
                for (l = 0; l < NP; l++) {
                  dudx00 = 0.0;
                  dvdy00 = 0.0;
                  for (i = 0; i < NP; i++) {
                    int pos_Dvv = l*NP + i;
                    int pos_gv_1 = j*NP + i;
                    int pos_gv_2 = NP*NP + i*NP + j;
                    dudx00 = dudx00 + Dvv[pos_Dvv]*gv[pos_gv_1];
                    dvdy00 = dvdy00 + Dvv[pos_Dvv]*gv[pos_gv_2];
                  }
                  int pos_dp_star = q*NLEV*NP*NP + k*NP*NP + j*NP + l;
                  int pos_vvtemp = l*NP + j;
                  dp_star[pos_dp_star] = dudx00;
                  vvtemp[pos_vvtemp] = dvdy00;
                }
              }
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_dp_star = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  int pos_vvtemp = j*NP + i;
                  int pos_det = j*NP + i;
                  dp_star[pos_dp_star] = (dp_star[pos_dp_star] + vvtemp[pos_vvtemp])   \
                      *(rmetdet[pos_det]*rrearth);
                }
              }   // end of divergence_sphere
              for (j = 0; j < NP; j++) {
                for (i = 0; i < NP; i++) {
                  int pos_Qtens = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  int pos_dp_star = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  int pos_qdp = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                  Qtens_temp[pos_Qtens] = Qdp[pos_qdp] - dt*dp_star[pos_dp_star];
                }
              }
              if (rhs_viss != 0) {
                for (j = 0; j < NP; j++) {
                  for (i = 0; i < NP; i++) {
                    int pos_Qtens = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                    Qtens_temp[pos_Qtens] = Qtens_temp[pos_Qtens] + Qtens_biharmonic[pos_Qtens];
                  }
                }
              }
            }   // end loop k
            if (limiter_option == 8) {
              for (k = 0; k < NLEV; k++) {
                for (j = 0; j < NP; j++) {
                  for (i = 0; i < NP; i++) {
                    int pos_dp_star = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                    int pos_dp = k*NP*NP + j*NP + i;
                    dp_star[pos_dp_star] = dp_temp[pos_dp] - dt*divdp[pos_dp];
                  }
                }
                if (nu_p > 0 && rhs_viss != 0) {
                  for (j = 0; j < NP; j++) {
                    for (i = 0; i < NP; i++) {
                      int pos_dpdiss = j*NP + i;
                      int pos_spheremp = j*NP + i;
                      int pos_dpdiss_bi = k*NP*NP + j*NP + i;
                      int pos_dp_star = q*NLEV*NP*NP + k*NP*NP + j*NP + i;
                      dpdiss[pos_dpdiss] = dpdiss_biharmonic[pos_dpdiss_bi];
                      dp_star[pos_dp_star] = dp_star[pos_dp_star]              \
                          - rhs_viss*dt*nu_q*dpdiss[pos_dpdiss]/spheremp[pos_spheremp];
                    }
                  }
                }
              }
            }   // end loop k
          }     // end loop q
          pe_put(src_dp_star, dp_star, cn*NLEV*NP*NP*sizeof(double));
          pe_put(src_Qtens_temp, Qtens_temp, cn*NLEV*NP*NP*sizeof(double));
          dma_syn();
        }
      }
      pe_put(src_dp_temp, dp_temp, block_dp*sizeof(double));
      dma_syn();
    }
  }

  if (id == 66) {
    src_Dinv = Dinv_ptr + 23*slice_qdp;
    for (i = 0; i < block_Dinv; i++) {
      printf("%lf\n", Dinv[i]);
    }
  }
}
