#include <slave.h>
#include <dma.h>


#define np_qdp 2
#define get_slv_id(tid) asm volatile ("rcsr %0, 0" : "=r"(tid))
#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))



typedef struct {
	int nets, nete, np, nlev, n0_qdp, np1_qdp, rkstage, limiter_opiton, qsize;
	long qdp_ptr1;
	long qdp_ptr2;
} param;

inline void dma_get_set(dma_desc *p_dma_get, volatile int *p_replyget, int size_block) {
	dma_set_op(p_dma_get, DMA_GET);
	dma_set_mode(p_dma_get, PE_MODE);
	dma_set_reply(p_dma_get, p_replyget);

	dma_set_size(p_dma_get, size_block);
	dma_set_bsize(p_dma_get, size_block);
	dma_set_stepsize(p_dma_get, 0);
}

inline void dma_put_set(dma_desc *p_dma_put, volatile int *p_replyput, int size_block) {
	dma_set_op(p_dma_put, DMA_PUT);
	dma_set_mode(p_dma_put, PE_MODE);
	dma_set_reply(p_dma_put, p_replyput);

	dma_set_size(p_dma_put, size_block);
	dma_set_bsize(p_dma_put, size_block);
	dma_set_stepsize(p_dma_put, 0);
}

#define dma_get(src, dst) {\
	dma(get_desc, (long)(src), (long)(dst));\
	dma_wait(&get_reply, 1);\
	get_reply = 0;\
}

#define dma_put(src, dst) {\
	dma(put_desc, (long)(src), (long)(dst));\
	dma_wait(&put_reply, 1);\
	put_reply = 0;\
}

void slave_qdp_time_avg_(param *param_s) {
	int id = athread_get_id(-1);
	int rid, cid;
	get_row_id(rid);
	get_col_id(cid);
	dma_desc get_desc, put_desc;
	volatile int get_reply = 0, put_reply = 0;

	param param_d;
	athread_get(PE_MODE, param_s, &param_d, sizeof(param), &get_reply, 0, 0, 0);
	int nets = param_d.nets;
	int nete = param_d.nete;
	int np = param_d.np;
	int nlev = param_d.nlev;
	int n0_qdp = param_d.n0_qdp;
	int np1_qdp = param_d.np1_qdp;
	int rkstage = param_d.rkstage;
	int limiter_opiton = param_d.limiter_opiton;
	int qsize = param_d.qsize;
	double *qdp_ptr1 = (double *)(param_d.qdp_ptr1);
	double *qdp_ptr2 = (double *)(param_d.qdp_ptr2);
	//long p_delta = (long)(qdp_ptr2) - (long)(qdp_ptr1);
	//int n_step = (int)(p_delta / sizeof(double));
	int size_blk = nlev * np * np;
	//double data[size_blk];
	//double *src, *dst;


	if (id == 0) {
		printf("%d, %d, %d, %d, %d, %d, %d, %d %d\n", nets, nete, np, nlev, n0_qdp, np1_qdp, rkstage, limiter_opiton, qsize);
	}


	//dma_get_set(&get_desc, &get_reply, size_blk * sizeof(double));
	//dma_put_set(&put_desc, &put_reply, size_blk * sizeof(double));
	////src = (double *)(qdp_ptr1 + n_step);
	//src = qdp_ptr1;
	//dst = (double *)(data);
	//if (id == 0) {
	//	int i;
	//	printf("%ld, %ld, %ld\n", (long)(qdp_ptr1), (long)(qdp_ptr2), (long)(src));
	//	dma_get(src, dst);
	//	for (i = 0; i < size_blk; i++) {
	//    printf("%lf\n", dst[i]);
	//	}
	//}



}
