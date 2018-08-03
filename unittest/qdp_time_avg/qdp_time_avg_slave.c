#include <slave.h>

typedef struct {
	long qdp_ptr1, qdp_ptr2;
	int nets, nete, np, nlev, n0_qdp, np1_qdp, rkstage, limiter_option, qsize;
} param;

#define athread_get_(src, dst, size_block, n_reply) {\
	get_reply = 0;\
	athread_get(PE_MODE, src, dst, size_block, &get_reply, 0, 0, 0);\
	while(get_reply!=n_reply);\
	asm volatile("memb\n\t");\
}

void slave_qdp_time_avg_(param *param_s) {
	volatile int id = athread_get_id(-1);
	volatile unsigned long get_reply, put_reply;

	long qdp_ptr1 = param_s->qdp_ptr1;
	long qdp_ptr2 = param_s->qdp_ptr2;
	int nets = param_s->nets;
	int nete = param_s->nete;
	int np = param_s->np;
	int nlev = param_s->nlev;
	int n0_qdp = param_s->n0_qdp;
	int np1_qdp = param_s->np1_qdp;
	int rkstage = param_s->rkstage;
	int limiter_option = param_s->limiter_option;
	int qsize = param_s->qsize;

	int size_blk = nlev * np * np;
	double data[size_blk];
	double *src, *dst;
	long ptr_delta = qdp_ptr2 - qdp_ptr1;
	int n_step = (int)(ptr_delta/sizeof(double));
	src = (double *)(qdp_ptr1) + n_step;
	dst = (double *)(data);
	if (id == 0)
		printf("%d, %ld, %ld\n",n_step, qdp_ptr2, src);
	if (id == 0) {
		printf("%d,%d,%d,%d,%d,%d,%d,%d,%d\n", nets, nete, np, nlev, n0_qdp, np1_qdp, rkstage, limiter_option, qsize);
		athread_get_(src, dst, (size_blk * sizeof(double)), 1);
		int i;
		for (i = 0; i < size_blk; i++)
			printf("%lf\n", data[i]);
	}
}
