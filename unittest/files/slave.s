# 1 "src/slave_euler_v.c"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "src/slave_euler_v.c"
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave.h" 1 3


# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 1 3 4
# 29 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/features.h" 1 3 4
# 35 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/features.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_config.h" 1 3 4
# 36 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/features.h" 2 3 4

# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_arch_features.h" 1 3 4
# 38 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/features.h" 2 3 4
# 356 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/features.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/sys/cdefs.h" 1 3 4
# 357 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/features.h" 2 3 4
# 30 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4



# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigset.h" 1 3 4
# 23 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigset.h" 3 4
typedef int __sig_atomic_t;




typedef struct
  {
    unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
  } __sigset_t;
# 103 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigset.h" 3 4
extern int __sigismember (__const __sigset_t *, int);
extern int __sigaddset (__sigset_t *, int);
extern int __sigdelset (__sigset_t *, int);
# 34 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4







typedef __sig_atomic_t sig_atomic_t;








typedef __sigset_t sigset_t;






# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/types.h" 1 3 4
# 28 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/types.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/wordsize.h" 1 3 4
# 29 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/types.h" 2 3 4


# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stddef.h" 1 3 4
# 214 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stddef.h" 3 4
typedef long unsigned int size_t;
# 32 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/types.h" 2 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/kernel_types.h" 1 3 4
# 10 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/kernel_types.h" 3 4
typedef unsigned int __kernel_dev_t;
typedef unsigned int __kernel_ino_t;
typedef unsigned int __kernel_mode_t;
typedef unsigned int __kernel_nlink_t;
typedef long __kernel_off_t;
typedef long __kernel_loff_t;
typedef int __kernel_pid_t;
typedef int __kernel_ipc_pid_t;
typedef unsigned int __kernel_uid_t;
typedef unsigned int __kernel_gid_t;
typedef unsigned long __kernel_size_t;
typedef long __kernel_ssize_t;
typedef long __kernel_ptrdiff_t;
typedef long __kernel_time_t;
typedef long __kernel_suseconds_t;
typedef long __kernel_clock_t;
typedef int __kernel_daddr_t;
typedef char * __kernel_caddr_t;
typedef unsigned long __kernel_sigset_t;
typedef unsigned short __kernel_uid16_t;
typedef unsigned short __kernel_gid16_t;
typedef __kernel_uid_t __kernel_old_uid_t;
typedef __kernel_gid_t __kernel_old_gid_t;
typedef __kernel_uid_t __kernel_uid32_t;
typedef __kernel_gid_t __kernel_gid32_t;
typedef __kernel_dev_t __kernel_old_dev_t;

typedef struct {
 int val[2];
} __kernel_fsid_t;
# 33 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/types.h" 2 3 4


typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;


typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;

typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;







typedef long int __quad_t;
typedef unsigned long int __u_quad_t;
# 135 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/types.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/typesizes.h" 1 3 4
# 136 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/types.h" 2 3 4


typedef unsigned long int __dev_t;
typedef unsigned int __uid_t;
typedef unsigned int __gid_t;
typedef unsigned int __ino_t;
typedef unsigned long int __ino64_t;
typedef unsigned int __mode_t;
typedef unsigned int __nlink_t;
typedef long int __off_t;
typedef long int __off64_t;
typedef int __pid_t;
typedef struct { int __val[2]; } __fsid_t;
typedef long int __clock_t;
typedef unsigned long int __rlim_t;
typedef unsigned long int __rlim64_t;
typedef unsigned int __id_t;
typedef long int __time_t;
typedef unsigned int __useconds_t;
typedef long int __suseconds_t;

typedef int __daddr_t;
typedef long int __swblk_t;
typedef int __key_t;


typedef int __clockid_t;


typedef void * __timer_t;


typedef unsigned int __blksize_t;




typedef unsigned int __blkcnt_t;
typedef unsigned long int __blkcnt64_t;


typedef int __fsblkcnt_t;
typedef long int __fsblkcnt64_t;


typedef unsigned int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;

typedef long int __ssize_t;



typedef __off64_t __loff_t;
typedef __quad_t *__qaddr_t;
typedef char *__caddr_t;


typedef long int __intptr_t;


typedef unsigned int __socklen_t;





typedef __kernel_ipc_pid_t __ipc_pid_t;
# 58 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/signum.h" 1 3 4
# 59 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4
# 75 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
typedef void (*__sighandler_t) (int);
# 91 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4


extern __sighandler_t signal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__));
# 105 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4

# 118 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
extern int kill (__pid_t __pid, int __sig) __attribute__ ((__nothrow__));






extern int killpg (__pid_t __pgrp, int __sig) __attribute__ ((__nothrow__));




extern int raise (int __sig) __attribute__ ((__nothrow__));




extern __sighandler_t ssignal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__));
extern int gsignal (int __sig) __attribute__ ((__nothrow__));




extern void psignal (int __sig, __const char *__s);
# 154 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
extern int __sigpause (int __sig_or_mask, int __is_sig);
# 179 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
extern int sigblock (int __mask) __attribute__ ((__nothrow__)) __attribute__ ((__deprecated__));


extern int sigsetmask (int __mask) __attribute__ ((__nothrow__)) __attribute__ ((__deprecated__));


extern int siggetmask (void) __attribute__ ((__nothrow__)) __attribute__ ((__deprecated__));
# 199 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
typedef __sighandler_t sig_t;







# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/time.h" 1 3 4
# 121 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/time.h" 3 4
struct timespec
  {
    __time_t tv_sec;
    long int tv_nsec;
  };
# 208 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4


# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/siginfo.h" 1 3 4
# 31 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/siginfo.h" 3 4
typedef union sigval
  {
    int sival_int;
    void *sival_ptr;
  } sigval_t;
# 45 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/siginfo.h" 3 4
typedef struct siginfo
  {
    int si_signo;
    int si_errno;

    int si_code;

    union
      {
 int _pad[((128 / sizeof (int)) - 4)];


 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
   } _kill;


 struct
   {
     int si_tid;
     int si_overrun;
     sigval_t si_sigval;
   } _timer;


 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
     sigval_t si_sigval;
   } _rt;


 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
     int si_status;
     __clock_t si_utime;
     __clock_t si_stime;
   } _sigchld;


 struct
   {
     void *si_addr;
   } _sigfault;


 struct
   {
     int si_band;
     int si_fd;
   } _sigpoll;
      } _sifields;
  } siginfo_t;
# 123 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/siginfo.h" 3 4
enum
{
  SI_ASYNCNL = -60,

  SI_TKILL = -6,

  SI_SIGIO,

  SI_ASYNCIO,

  SI_MESGQ,

  SI_TIMER,

  SI_QUEUE,

  SI_USER,

  SI_KERNEL = 0x80

};



enum
{
  ILL_ILLOPC = 1,

  ILL_ILLOPN,

  ILL_ILLADR,

  ILL_ILLTRP,

  ILL_PRVOPC,

  ILL_PRVREG,

  ILL_COPROC,

  ILL_BADSTK

};


enum
{
  FPE_INTDIV = 1,

  FPE_INTOVF,

  FPE_FLTDIV,

  FPE_FLTOVF,

  FPE_FLTUND,

  FPE_FLTRES,

  FPE_FLTINV,

  FPE_FLTSUB

};


enum
{
  SEGV_MAPERR = 1,

  SEGV_ACCERR

};


enum
{
  BUS_ADRALN = 1,

  BUS_ADRERR,

  BUS_OBJERR

};


enum
{
  TRAP_BRKPT = 1,

  TRAP_TRACE

};


enum
{
  CLD_EXITED = 1,

  CLD_KILLED,

  CLD_DUMPED,

  CLD_TRAPPED,

  CLD_STOPPED,

  CLD_CONTINUED

};


enum
{
  POLL_IN = 1,

  POLL_OUT,

  POLL_MSG,

  POLL_ERR,

  POLL_PRI,

  POLL_HUP

};
# 263 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/siginfo.h" 3 4
typedef struct sigevent
  {
    sigval_t sigev_value;
    int sigev_signo;
    int sigev_notify;

    union
      {
 int _pad[((64 / sizeof (int)) - 4)];



 __pid_t _tid;

 struct
   {
     void (*_function) (sigval_t);
     void *_attribute;
   } _sigev_thread;
      } _sigev_un;
  } sigevent_t;






enum
{
  SIGEV_SIGNAL = 0,

  SIGEV_NONE,

  SIGEV_THREAD,


  SIGEV_THREAD_ID = 4

};
# 211 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4



extern int sigemptyset (sigset_t *__set) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int sigfillset (sigset_t *__set) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int sigaddset (sigset_t *__set, int __signo) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int sigdelset (sigset_t *__set, int __signo) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int sigismember (__const sigset_t *__set, int __signo)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
# 244 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigaction.h" 1 3 4
# 25 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigaction.h" 3 4
struct sigaction
  {


    union
      {

 __sighandler_t sa_handler;

 void (*sa_sigaction) (int, siginfo_t *, void *);
      }
    __sigaction_handler;







    __sigset_t sa_mask;


    unsigned int sa_flags;
  };
# 245 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4


extern int sigprocmask (int __how, __const sigset_t *__restrict __set,
   sigset_t *__restrict __oset) __attribute__ ((__nothrow__));






extern int sigsuspend (__const sigset_t *__set) __attribute__ ((__nonnull__ (1)));


extern int sigaction (int __sig, __const struct sigaction *__restrict __act,
        struct sigaction *__restrict __oact) __attribute__ ((__nothrow__));


extern int sigpending (sigset_t *__set) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));






extern int sigwait (__const sigset_t *__restrict __set, int *__restrict __sig)
     __attribute__ ((__nonnull__ (1, 2)));
# 308 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
struct sigvec
  {
    __sighandler_t sv_handler;
    int sv_mask;

    int sv_flags;

  };
# 328 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
extern int sigvec (int __sig, __const struct sigvec *__vec,
     struct sigvec *__ovec) __attribute__ ((__nothrow__));



# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigcontext.h" 1 3 4
# 28 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigcontext.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/asm/sigcontext.h" 1 3 4



struct sigcontext {
# 13 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/asm/sigcontext.h" 3 4
  long sc_onstack;
  long sc_mask;
  long sc_pc;
  long sc_ps;
  long sc_regs[32];
  long sc_ownedfp;
  long sc_fpregs[32];
  unsigned long sc_fpcr;
  unsigned long sc_fp_control;
  unsigned long sc_reserved1, sc_reserved2;
  unsigned long sc_ssize;
  char * sc_sbase;
  unsigned long sc_traparg_a0;
  unsigned long sc_traparg_a1;
  unsigned long sc_traparg_a2;
  unsigned long sc_fp_trap_pc;
  unsigned long sc_fp_trigger_sum;
  unsigned long sc_fp_trigger_inst;
};
# 29 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigcontext.h" 2 3 4
# 334 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4


extern int sigreturn (struct sigcontext *__scp) __attribute__ ((__nothrow__));
# 346 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
extern int siginterrupt (int __sig, int __interrupt) __attribute__ ((__nothrow__));

# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigstack.h" 1 3 4
# 26 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigstack.h" 3 4
struct sigstack
  {
    void * ss_sp;
    int ss_onstack;
  };



enum
{
  SS_ONSTACK = 1,

  SS_DISABLE

};
# 50 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/sigstack.h" 3 4
typedef struct sigaltstack
  {
    void * ss_sp;
    int ss_flags;
    size_t ss_size;
  } stack_t;
# 349 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 2 3 4
# 357 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
extern int sigstack (struct sigstack *__ss, struct sigstack *__oss)
     __attribute__ ((__nothrow__)) __attribute__ ((__deprecated__));



extern int sigaltstack (__const struct sigaltstack *__restrict __ss,
   struct sigaltstack *__restrict __oss) __attribute__ ((__nothrow__));
# 394 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/signal.h" 3 4
extern int __libc_current_sigrtmin (void) __attribute__ ((__nothrow__));

extern int __libc_current_sigrtmax (void) __attribute__ ((__nothrow__));




# 4 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave.h" 2 3

# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/share.h" 1 3



# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stddef.h" 1 3 4
# 152 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stddef.h" 3 4
typedef long int ptrdiff_t;
# 326 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stddef.h" 3 4
typedef int wchar_t;
# 5 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/share.h" 2 3
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 1 3 4
# 30 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4




# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stddef.h" 1 3 4
# 35 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 2 3 4
# 44 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4


typedef struct __STDIO_FILE_STRUCT FILE;





# 62 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
typedef struct __STDIO_FILE_STRUCT __FILE;
# 72 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_stdio.h" 1 3 4
# 119 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_stdio.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_mutex.h" 1 3 4
# 120 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_stdio.h" 2 3 4
# 170 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_stdio.h" 3 4
typedef struct {
 __off_t __pos;






} __STDIO_fpos_t;


typedef struct {
 __off64_t __pos;






} __STDIO_fpos64_t;




typedef __off64_t __offmax_t;
# 233 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_stdio.h" 3 4
struct __STDIO_FILE_STRUCT {
 unsigned short __modeflags;







 unsigned char __ungot[2];

 int __filedes;
# 284 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_stdio.h" 3 4
};
# 384 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/uClibc_stdio.h" 3 4
extern int __fgetc_unlocked(FILE *__stream);
extern int __fputc_unlocked(int __c, FILE *__stream);
# 73 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 2 3 4



# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stdarg.h" 1 3 4
# 43 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/include/stdarg.h" 3 4
typedef __builtin_va_list __gnuc_va_list;
# 77 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 2 3 4




typedef __STDIO_fpos_t fpos_t;




# 131 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/bits/stdio_lim.h" 1 3 4
# 132 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 2 3 4
# 144 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern FILE *_slave_stdin;
extern FILE *_slave_stdout;
extern FILE *_slave_stderr;
# 158 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4


extern int remove (__const char *__filename) __attribute__ ((__nothrow__));

extern int rename (__const char *__old, __const char *__new) __attribute__ ((__nothrow__));









extern FILE *tmpfile (void) ;
# 186 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern char *tmpnam (char *__s) __attribute__ ((__nothrow__)) ;





extern char *tmpnam_r (char *__s) __attribute__ ((__nothrow__)) ;
# 204 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern char *tempnam (__const char *__dir, __const char *__pfx)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;








extern int fclose (FILE *__stream);




extern int fflush (FILE *__stream);

# 229 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern int fflush_unlocked (FILE *__stream);
# 243 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4






extern FILE *fopen (__const char *__restrict __filename,
      __const char *__restrict __modes) ;




extern FILE *freopen (__const char *__restrict __filename,
        __const char *__restrict __modes,
        FILE *__restrict __stream) ;
# 272 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4

# 283 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern FILE *fdopen (int __fd, __const char *__modes) __attribute__ ((__nothrow__)) ;
# 307 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4



extern void setbuf (FILE *__restrict __stream, char *__restrict __buf) __attribute__ ((__nothrow__));



extern int setvbuf (FILE *__restrict __stream, char *__restrict __buf,
      int __modes, size_t __n) __attribute__ ((__nothrow__));





extern void setbuffer (FILE *__restrict __stream, char *__restrict __buf,
         size_t __size) __attribute__ ((__nothrow__));


extern void setlinebuf (FILE *__stream) __attribute__ ((__nothrow__));








extern int fprintf (FILE *__restrict __stream,
      __const char *__restrict __format, ...);




extern int printf (__const char *__restrict __format, ...);

extern int sprintf (char *__restrict __s,
      __const char *__restrict __format, ...) __attribute__ ((__nothrow__));





extern int vfprintf (FILE *__restrict __s, __const char *__restrict __format,
       __gnuc_va_list __arg);




extern int vprintf (__const char *__restrict __format, __gnuc_va_list __arg);

extern int vsprintf (char *__restrict __s, __const char *__restrict __format,
       __gnuc_va_list __arg) __attribute__ ((__nothrow__));





extern int snprintf (char *__restrict __s, size_t __maxlen,
       __const char *__restrict __format, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 4)));

extern int vsnprintf (char *__restrict __s, size_t __maxlen,
        __const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 0)));

# 403 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4





extern int fscanf (FILE *__restrict __stream,
     __const char *__restrict __format, ...) ;




extern int scanf (__const char *__restrict __format, ...) ;

extern int sscanf (__const char *__restrict __s,
     __const char *__restrict __format, ...) __attribute__ ((__nothrow__));

# 445 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4





extern int fgetc (FILE *__stream);
extern int getc (FILE *__stream);





extern int getchar (void);

# 469 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern int getc_unlocked (FILE *__stream);
extern int getchar_unlocked (void);
# 483 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern int fgetc_unlocked (FILE *__stream);











extern int fputc (int __c, FILE *__stream);
extern int putc (int __c, FILE *__stream);





extern int putchar (int __c);

# 516 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern int fputc_unlocked (int __c, FILE *__stream);







extern int putc_unlocked (int __c, FILE *__stream);
extern int putchar_unlocked (int __c);
# 535 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern int getw (FILE *__stream);


extern int putw (int __w, FILE *__stream);








extern char *fgets (char *__restrict __s, int __n, FILE *__restrict __stream)
     ;






extern char *gets (char *__s) ;

# 602 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4





extern int fputs (__const char *__restrict __s, FILE *__restrict __stream);





extern int puts (__const char *__s);






extern int ungetc (int __c, FILE *__stream);






extern size_t fread (void *__restrict __ptr, size_t __size,
       size_t __n, FILE *__restrict __stream) ;




extern size_t fwrite (__const void *__restrict __ptr, size_t __size,
        size_t __n, FILE *__restrict __s) ;

# 655 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern size_t fread_unlocked (void *__restrict __ptr, size_t __size,
         size_t __n, FILE *__restrict __stream) ;
extern size_t fwrite_unlocked (__const void *__restrict __ptr, size_t __size,
          size_t __n, FILE *__restrict __stream) ;








extern int fseek (FILE *__stream, long int __off, int __whence);




extern long int ftell (FILE *__stream) ;




extern void rewind (FILE *__stream);

# 710 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4






extern int fgetpos (FILE *__restrict __stream, fpos_t *__restrict __pos);




extern int fsetpos (FILE *__stream, __const fpos_t *__pos);
# 733 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4

# 742 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4


extern void clearerr (FILE *__stream) __attribute__ ((__nothrow__));

extern int feof (FILE *__stream) __attribute__ ((__nothrow__)) ;

extern int ferror (FILE *__stream) __attribute__ ((__nothrow__)) ;




extern void clearerr_unlocked (FILE *__stream) __attribute__ ((__nothrow__));
extern int feof_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;
extern int ferror_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;








extern void perror (__const char *__s);

# 779 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern int fileno (FILE *__stream) __attribute__ ((__nothrow__)) ;




extern int fileno_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;
# 794 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern FILE *popen (__const char *__command, __const char *__modes) ;





extern int pclose (FILE *__stream);





extern char *ctermid (char *__s) __attribute__ ((__nothrow__));
# 834 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4
extern void flockfile (FILE *__stream) __attribute__ ((__nothrow__));



extern int ftrylockfile (FILE *__stream) __attribute__ ((__nothrow__)) ;


extern void funlockfile (FILE *__stream) __attribute__ ((__nothrow__));
# 896 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/sys-include/stdio.h" 3 4

# 6 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/share.h" 2 3
# 49 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/share.h" 3
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef unsigned short uint16_t;
typedef unsigned long uint64_t;


typedef int int32_t;
typedef short int16_t;
typedef long int64_t;






enum {
        SW3_SIGSYSC = 1,
        SW3_SIGHALT,
        SW3_SIGEND,
 SW3_SIGERR,
 SW3_SIGPRINT = 11,
 SW3_SIGEXPT = 62,
        SW3_SIGMAX = 63,
};
# 171 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/share.h" 3
enum Exception_Kind {
        UNALIGN = 0,
        ILLEGAL,
        OVI,
 INE,
        UNF,
 OVF,
 DBZ,
        INV,
 DNO,
 CLASSE1,
 CLASSE2,
        SELLDWE,
 LDME,
 SDLBE,
 MABC,
 MATNO,
 RAMAR,
        IFLOWE,
        SBMDE1,
        SBMDE2,
        SYNE1,
        SYNE2,
 RCE,
 DMAE1,
 DMAE2,
 DMAE3,
 DMAE4,
 DMAE5,
        IOE1,
        IOE2,
        IOE3,
 OTHERE,
 _SYS_NSIG,
};

typedef struct slave_expt {
    int coreno;
    unsigned long expt_vector;
    unsigned long expt_pc;
    unsigned long dma_expt_type;
    unsigned long dma0;
    unsigned long dma1;
    unsigned long dma2;
    unsigned long dma3;
    unsigned long tc_sdlb_err_spot;
    unsigned long reserve_data;
} slave_expt;
# 227 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/share.h" 3
extern long *_tdata_local_start __attribute__ ((weak));
extern long *_tdata_local_end __attribute__ ((weak));
extern long *_tdata_private_start __attribute__ ((weak));
extern long *_tdata_private_end __attribute__ ((weak));
extern long *_tdata_local_fix_end __attribute__ ((weak));

typedef struct s_thread_info{
        char valid;
        int thread_id;
        int core_num;
        int state_flag;
        void * pc;
        void * arg;
        char fini_sig;
        long gid;
        int team_size;
}thread_info_t;

typedef struct s_core_type {
 unsigned long coremask_fault;
 unsigned long coremask_employ;
 unsigned long coremask_free;
 unsigned long spe_access;
}core_type_t;

typedef struct s_msg{
        unsigned long own;
        int type;
        char * msg;
        unsigned long msg_size;
}msg_t;


extern void slave_start();
# 279 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/share.h" 3
extern char __id_core_map[4][64];
extern char __core_id_map[4][64];
extern thread_info_t * __v_core_thread_map[4][64];


extern volatile unsigned long *__total_tasks;
extern int __cgid;
extern int __cgnum;
# 6 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave.h" 2 3
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave_intc.h" 1 3
# 10 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave_intc.h" 3
extern void _Waiting_For_Task();
extern void _MPE_stop();
# 7 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave.h" 2 3
# 1 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave_sig.h" 1 3
# 8 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave.h" 2 3
# 118 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave.h" 3
typedef enum {
        PE_MODE,
        BCAST_MODE,
        ROW_MODE,
        BROW_MODE,
        RANK_MODE
} dma_mode;

typedef enum {
        DMA_PUT,
        DMA_GET,
        DMA_PUT_P,
        DMA_GET_P,
        DMA_BARRIER = 5
} DMA_OP;

typedef enum {
 ROW_SCOPE,
 COL_SCOPE,
 ARRAY_SCOPE,
}scope;
# 153 "/usr/sw-mpp/swcc/sw5gcc-binary/bin/../lib/gcc/sw_64thl-unknown-linux-gnu/4.1.2/../../../../sw_64thl-unknown-linux-gnu/include/slave.h" 3
extern int athread_get(dma_mode mode, void *src, void *dest, int len, void *reply, char mask,int stride,int bsize);
extern int athread_put(dma_mode mode, void *src, void *dest, int len, void *reply, int stride,int bsize);

extern int athread_get_core( int id);
extern int athread_get_id(int core);
extern int athread_syn(scope scp,int mask);





extern __attribute__ ((section (".tdata_local_fix"))) char _CGN,_ROW,_COL,_PEN;
extern __attribute__ ((section (".tdata_local_fix"))) volatile void (*_PC)();
extern __attribute__ ((section (".tdata_local_fix"))) int _MYID;
# 2 "src/slave_euler_v.c" 2
# 1 "src/dma_macros.h" 1
# 3 "src/slave_euler_v.c" 2
# 52 "src/slave_euler_v.c"
typedef struct {
  double *gl_qdp, *gl_qdp_leap, *divdp_proj, *dp, *vn0, *dp_temp, *dp_star_temp, *Dvv, *Dinv, *metdet, *rmetdet, *Qtens_temp, *Qtens_biharmonic, *divdp, *dpdiss_biharmonic, *spheremp, *qmax, *qmin;



  double dt, rrearth, nu_p, nu_q;
  int nets, nete, rhs_multiplier, qsize, n0_qdp, np1_qdp, limiter_option , rhs_viss;

} param_t;

void slave_euler_v_(param_t *param_s) {
  volatile int id = athread_get_id(-1);
  volatile unsigned long get_reply, put_reply;
  volatile int cid, rid;
 volatile int reply_shadow = 0; int count_shadow = 0; long pe_get_stride_shadow = 0; long pe_put_stride_shadow = 0; long row_get_stride_shadow = 0; long row_put_stride_shadow = 0; long rank_get_stride_shadow = 0; long rank_put_stride_shadow = 0; long brow_get_stride_shadow = 0; long bcast_get_stride_shadow = 0; long pe_get_bsize_shadow = 0; long pe_put_bsize_shadow = 0; long row_get_bsize_shadow = 0; long row_put_bsize_shadow = 0; long rank_get_bsize_shadow = 0; long rank_put_bsize_shadow = 0; long brow_get_bsize_shadow = 0; long bcast_get_bsize_shadow = 0; int bcast_get_mask_shadow = 0xff; int brow_get_mask_shadow = 0xff;;

  param_t param_d;
  { athread_get(PE_MODE, (param_s), (&param_d), (sizeof(param_t)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
  { while (reply_shadow != count_shadow) { }; asm volatile("memb\n\t"); };

  double *gl_qdp = param_d.gl_qdp;
  double *gl_qdp_leap = param_d.gl_qdp_leap;
  double *gl_divdp_proj = param_d.divdp_proj;
  double *gl_dp = param_d.dp;
  double *gl_vn0 = param_d.vn0;
  double *gl_dp_temp = param_d.dp_temp;
  double *gl_dp_star = param_d.dp_star_temp;
  double *gl_Dvv = param_d.Dvv;
  double *gl_Dinv = param_d.Dinv;
  double *gl_metdet = param_d.metdet;
  double *gl_rmetdet = param_d.rmetdet;
  double *gl_Qtens_temp = param_d.Qtens_temp;
  double *gl_Qtens_biharmonic = param_d.Qtens_biharmonic;
  double *gl_divdp = param_d.divdp;
  double *gl_dpdiss_biharmonic = param_d.dpdiss_biharmonic;
  double *gl_spheremp = param_d.spheremp;
  double *gl_qmax = param_d.qmax;
  double *gl_qmin = param_d.qmin;
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

  rid = id / 4;
  cid = id % 4;
  int loop_r = ((nete - nets + 1) + 3*16 - 1)/(3*16);
  int loop_c = (qsize + 1*4 - 1)/(1*4);
  int slice_qdp = (int)(gl_qdp_leap - gl_qdp);
  int istep_dp_star = qsize*30*4*4;
  int istep_Qtens = qsize*30*4*4;
  int istep_qmax = qsize*30;
  int c, r, i, j, l, k, k1, iter, q, ie, cbeg, cend, rbeg, rend, cn, rn;

  double Qdp[(1*30*4*4)];
  double dp[(30*4*4)];
  double dp_temp[(30*4*4)];
  double dp_star[(1*30*4*4)];
  double vn0[(30*2*4*4)];
  double Vstar[(30*2*4*4)];
  double divdp_proj[(30*4*4)];
  double gradQ[(2*4*4)];
  double Dvv[4*4];
  double Dinv[2*2*4*4];
  double metdet[4*4];
  double rmetdet[4*4];
  double Qtens_temp[(1*30*4*4)];
  double Qtens_biharmonic[(1*30*4*4)];
  double dpdiss_biharmonic[(30*4*4)];
  double divdp[(30*4*4)];
  double dpdiss[4*4];
  double spheremp[4*4];
  double qmax[1*30];
  double qmin[1*30];
# 144 "src/slave_euler_v.c"
  { athread_get(PE_MODE, (gl_Dvv), (Dvv), (4*4*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
  { while (reply_shadow != count_shadow) { }; asm volatile("memb\n\t"); };


  double dudx00, dvdy00;
  double gv[2*4*4], vvtemp[4*4];

  double cc[4*4], xx[4*4];
  double addmass, weightssum;
  double mass = 0.0;
  double sumc = 0.0;
  double tol_limiter = 5e-14;

  double *src_qdp, *src_dp, *src_dp_temp, *src_vn0, *src_divdp_proj, *src_Dinv, *src_metdet, *src_rmetdet, *src_dp_star, *src_Qtens_temp, *src_Qtens_biharmonic, *src_divdp, *src_dpdiss_biharmonic, *src_spheremp, *src_qmax, *src_qmin;



  double *src_n0_qdp = gl_qdp + (n0_qdp - 1)*qsize*(30*4*4);
  double *src_np1_qdp = gl_qdp + (np1_qdp - 1)*qsize*(30*4*4);



  for (r = 0; r < loop_r; r++) {
    rbeg = r*16*3 + rid*3;
    rend = r*16*3 + (rid + 1)*3;
    rend = rend < (nete - nets + 1) ? rend : (nete - nets + 1);
    rn = rend - rbeg;
    rn = rn < 0 ? 0 : rn;
    for (ie = 0; ie < rn; ie++) {
      src_dp = gl_dp + (rbeg + ie)*slice_qdp;

      src_vn0 = gl_vn0 + (rbeg + ie)*slice_qdp;
      src_divdp_proj = gl_divdp_proj + (rbeg + ie)*slice_qdp;
      src_divdp = gl_divdp + (rbeg + ie)*slice_qdp;
      src_dpdiss_biharmonic = gl_dpdiss_biharmonic + (rbeg + ie)*slice_qdp;
      src_Dinv = gl_Dinv + (rbeg + ie)*slice_qdp;
      src_metdet = gl_metdet + (rbeg + ie)*slice_qdp;
      src_rmetdet = gl_rmetdet + (rbeg + ie)*slice_qdp;
      src_spheremp = gl_spheremp + (rbeg + ie)*slice_qdp;
      { athread_get(PE_MODE, (src_dp), (dp), ((30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };

      { athread_get(PE_MODE, (src_vn0), (vn0), ((30*2*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { athread_get(PE_MODE, (src_divdp_proj), (divdp_proj), ((30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { athread_get(PE_MODE, (src_divdp), (divdp), ((30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { athread_get(PE_MODE, (src_dpdiss_biharmonic), (dpdiss_biharmonic), ((30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { athread_get(PE_MODE, (src_Dinv), (Dinv), ((2*2*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { athread_get(PE_MODE, (src_metdet), (metdet), ((4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { athread_get(PE_MODE, (src_rmetdet), (rmetdet), ((4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { athread_get(PE_MODE, (src_spheremp), (spheremp), (4*4*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
      { while (reply_shadow != count_shadow) { }; asm volatile("memb\n\t"); };

      for (k = 0; k < 30; k++) {
        for (j = 0; j < 4; j++) {
          for (i = 0; i < 4; i++) {
            int pos_dp = k*4*4 + j*4 + i;
            int pos_vn0_1 = k*2*4*4 + j*4 + i;
            int pos_vn0_2 = k*2*4*4 + 4*4 + j*4 + i;
            dp_temp[pos_dp] = dp[pos_dp] - rhs_multiplier*dt*divdp_proj[pos_dp];
            Vstar[pos_vn0_1] = vn0[pos_vn0_1]/dp_temp[pos_dp];
            Vstar[pos_vn0_2] = vn0[pos_vn0_2]/dp_temp[pos_dp];
          }
        }
      }
      for (c = 0; c < loop_c; c++) {
        cbeg = c*4*1 + cid*1;
        cend = c*4*1 + (cid + 1)*1;
        cend = cend < qsize ? cend : qsize;
        cn = cend - cbeg;
        if (cn > 0) {
          src_qdp = src_n0_qdp + (rbeg + ie)*slice_qdp + cbeg*(30*4*4);
          src_Qtens_temp = gl_Qtens_temp + (rbeg + ie)*istep_Qtens + cbeg*(30*4*4);
          src_Qtens_biharmonic = gl_Qtens_biharmonic + (rbeg + ie)*istep_Qtens + cbeg*(30*4*4);

          src_dp_star = gl_dp_star + (rbeg + ie)*istep_dp_star + cbeg*(30*4*4);
          src_qmax = gl_qmax + (rbeg + ie)*istep_qmax + cbeg*30;
          src_qmin = gl_qmin + (rbeg + ie)*istep_qmax + cbeg*30;
          { athread_get(PE_MODE, (src_qdp), (Qdp), ((1*30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
          { athread_get(PE_MODE, (src_Qtens_temp), (Qtens_temp), ((1*30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
          { athread_get(PE_MODE, (src_Qtens_biharmonic), (Qtens_biharmonic), ((1*30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
          { athread_get(PE_MODE, (src_dp_star), (dp_star), ((1*30*4*4)*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
          { athread_get(PE_MODE, (src_qmax), (qmax), (1*30*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
          { athread_get(PE_MODE, (src_qmin), (qmin), (1*30*sizeof(double)), (void*)&reply_shadow , 0, pe_get_stride_shadow, pe_get_bsize_shadow); count_shadow ++; };
          { while (reply_shadow != count_shadow) { }; asm volatile("memb\n\t"); };
          for (q = 0; q < cn; q++) {
            for (k = 0; k < 30; k++) {
              for (j = 0; j < 4; j++) {
                for (i = 0; i < 4; i++) {
                  int pos_qdp = q*30*4*4 + k*4*4 + j*4 + i;
                  int pos_vn0_1 = k*2*4*4 + j*4 + i;
                  int pos_vn0_2 = k*2*4*4 + 4*4 + j*4 + i;
                  int pos_gradQ_1 = j*4 + i;
                  int pos_gradQ_2 = 4*4 + j*4 + i;
                  int pos_Dinv_1 = j*4 + i;
                  int pos_Dinv_2 = 4*4 + j*4 + i;
                  int pos_Dinv_3 = 2*4*4 + j*4 + i;
                  int pos_Dinv_4 = 3*4*4 + j*4 + i;
                  int pos_gv_1 = j*4 + i;
                  int pos_gv_2 = 4*4 + j*4 + i;
                  int pos_det = j*4 + i;
                  gradQ[pos_gradQ_1] = Vstar[pos_vn0_1]*Qdp[pos_qdp];
                  gradQ[pos_gradQ_2] = Vstar[pos_vn0_2]*Qdp[pos_qdp];
                  gv[pos_gv_1] = metdet[pos_det]*(Dinv[pos_Dinv_1]*gradQ[pos_gradQ_1] + Dinv[pos_Dinv_3]*gradQ[pos_gradQ_2]);

                  gv[pos_gv_2] = metdet[pos_det]*(Dinv[pos_Dinv_2]*gradQ[pos_gradQ_1] + Dinv[pos_Dinv_4]*gradQ[pos_gradQ_2]);

                }
              }
              for (j = 0; j < 4; j++) {
                for (l = 0; l < 4; l++) {
                  dudx00 = 0.0;
                  dvdy00 = 0.0;
                  for (i = 0; i < 4; i++) {
                    int pos_Dvv = l*4 + i;
                    int pos_gv_1 = j*4 + i;
                    int pos_gv_2 = 4*4 + i*4 + j;
                    dudx00 = dudx00 + Dvv[pos_Dvv]*gv[pos_gv_1];
                    dvdy00 = dvdy00 + Dvv[pos_Dvv]*gv[pos_gv_2];
                  }
                  int pos_dp_star = q*30*4*4 + k*4*4 + j*4 + l;
                  int pos_vvtemp = l*4 + j;
                  dp_star[pos_dp_star] = dudx00;
                  vvtemp[pos_vvtemp] = dvdy00;
                }
              }
              for (j = 0; j < 4; j++) {
                for (i = 0; i < 4; i++) {
                  int pos_dp_star = q*30*4*4 + k*4*4 + j*4 + i;
                  int pos_vvtemp = j*4 + i;
                  int pos_det = j*4 + i;
                  dp_star[pos_dp_star] = (dp_star[pos_dp_star] + vvtemp[pos_vvtemp]) *(rmetdet[pos_det]*rrearth);

                }
              }
              for (j = 0; j < 4; j++) {
                for (i = 0; i < 4; i++) {
                  int pos_Qtens = q*30*4*4 + k*4*4 + j*4 + i;
                  int pos_dp_star = q*30*4*4 + k*4*4 + j*4 + i;
                  int pos_qdp = q*30*4*4 + k*4*4 + j*4 + i;
                  Qtens_temp[pos_Qtens] = Qdp[pos_qdp] - dt*dp_star[pos_dp_star];
                }
              }
              if (rhs_viss != 0) {
                for (j = 0; j < 4; j++) {
                  for (i = 0; i < 4; i++) {
                    int pos_Qtens = q*30*4*4 + k*4*4 + j*4 + i;
                    Qtens_temp[pos_Qtens] = Qtens_temp[pos_Qtens] + Qtens_biharmonic[pos_Qtens];
                  }
                }
              }
            }
            if (limiter_option == 8) {
              for (k = 0; k < 30; k++) {
                for (j = 0; j < 4; j++) {
                  for (i = 0; i < 4; i++) {
                    int pos_dp_star = q*30*4*4 + k*4*4 + j*4 + i;
                    int pos_dp = k*4*4 + j*4 + i;
                    dp_star[pos_dp_star] = dp_temp[pos_dp] - dt*divdp[pos_dp];
                  }
                }
                if (nu_p > 0 && rhs_viss != 0) {
                  for (j = 0; j < 4; j++) {
                    for (i = 0; i < 4; i++) {
                      int pos_dpdiss = j*4 + i;
                      int pos_spheremp = j*4 + i;
                      int pos_dpdiss_bi = k*4*4 + j*4 + i;
                      int pos_dp_star = q*30*4*4 + k*4*4 + j*4 + i;
                      dpdiss[pos_dpdiss] = dpdiss_biharmonic[pos_dpdiss_bi];
                      dp_star[pos_dp_star] = dp_star[pos_dp_star] - rhs_viss*dt*nu_q*dpdiss[pos_dpdiss]/spheremp[pos_spheremp];

                    }
                  }
                }
              }

              for (k = 0; k < 30; k++) {
                for (k1 = 0; k1 < 4*4; k1++) {
                  int pos_dp_star = q*30*4*4 + k*4*4 + k1;
                  int pos_Qtens = q*30*4*4 + k*4*4 + k1;
                  cc[k1] = spheremp[k1]*dp_star[pos_dp_star];
                  xx[k1] = Qtens_temp[pos_Qtens]/dp_star[pos_dp_star];
                }



                { int i; for (i = 0; i < 4*4; i++) sumc = sumc + cc[i]; };
                if (sumc <= 0) continue;
                { int i; double _array[16]; for (i = 0; i < 16; i++) { _array[i] = cc[i] * xx[i]; mass = mass + _array[i]; } };
                int pos_qmax = q*30 + k;
                if (mass < qmin[pos_qmax]*sumc) qmin[pos_qmax] = mass/sumc;
                if (mass > qmax[pos_qmax]*sumc) qmax[pos_qmax] = mass/sumc;
                for (iter = 0; iter < (4*4 - 1); iter++) {
                  addmass = 0.0;
                  for (k1 = 0; k1 < 4*4; k1++) {
                    if(xx[k1] > qmax[pos_qmax]) {
                      addmass = addmass + (xx[k1] - qmax[pos_qmax])*cc[k1];
                      xx[k1] = qmax[pos_qmax];
                    }
                    if(xx[k1] < qmin[pos_qmax]) {
                      addmass = addmass - (qmin[pos_qmax] - xx[k1])*cc[k1];
                      xx[k1] = qmin[pos_qmax];
                    }
                  }
                  double addmass_abs;
                  double mass_abs;
                  asm volatile ("fcpys $31, %1, %0" : "=r"(mass_abs) : "r"(mass));
                  asm volatile ("fcpys $31, %1, %0" : "=r"(addmass_abs) : "r"(addmass));
                  if (addmass_abs <= tol_limiter*mass_abs) break;
                  weightssum = 0.0;
                  if (addmass > 0) {
                    for (k1 = 0; k1 < 4*4; k1++) {
                      if (xx[k1] < qmax[pos_qmax])
                        weightssum = weightssum + cc[k1];
                    }

                    for (k1 = 0; k1 < 4*4; k1++) {
                      if (xx[k1] < qmax[pos_qmax])
                        xx[k1] = xx[k1] + addmass/weightssum;
                    }

                  } else {
                    for (k1 = 0; k1 < 4*4; k1++) {
                      if (xx[k1] > qmin[pos_qmax])
                        weightssum = weightssum + cc[k1];
                    }

                    for (k1 = 0; k1 < 4*4; k1++) {
                      if (xx[k1] > qmin[pos_qmax])
                        xx[k1] = xx[k1] + addmass/weightssum;
                    }

                  }
                }
                for (k1 = 0; k1 < 4*4; k1++) {
                  int pos_Qtens = q*30*4*4 + k*4*4 + k1;
                  Qtens_temp[pos_Qtens] = xx[k1];
                }
              }
              for (k = 0; k < 30; k++) {
                for (k1 = 0; k1 < 4*4; k1++) {
                  int pos_Qtens = q*30*4*4 + k*4*4 + k1;
                  int pos_dp_star = q*30*4*4 + k*4*4 + k1;
                  Qtens_temp[pos_Qtens] = Qtens_temp[pos_Qtens]*dp_star[pos_dp_star];
                }
              }
            }
          }
          { athread_put(PE_MODE, (dp_star), (src_dp_star), (cn*30*4*4*sizeof(double)), (void*)&reply_shadow , pe_put_stride_shadow, pe_put_bsize_shadow); count_shadow ++; };
          { athread_put(PE_MODE, (Qtens_temp), (src_Qtens_temp), (cn*30*4*4*sizeof(double)), (void*)&reply_shadow , pe_put_stride_shadow, pe_put_bsize_shadow); count_shadow ++; };
          { athread_put(PE_MODE, (qmax), (src_qmax), (cn*30*sizeof(double)), (void*)&reply_shadow , pe_put_stride_shadow, pe_put_bsize_shadow); count_shadow ++; };
          { athread_put(PE_MODE, (qmin), (src_qmin), (cn*30*sizeof(double)), (void*)&reply_shadow , pe_put_stride_shadow, pe_put_bsize_shadow); count_shadow ++; };
          { while (reply_shadow != count_shadow) { }; asm volatile("memb\n\t"); };
        }
      }


    }
  }
# 415 "src/slave_euler_v.c"
}
