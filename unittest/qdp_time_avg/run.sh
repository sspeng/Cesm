#! /bin/bash
make
bsub -I -b -q q_sw_expr -host_stack 1024 -share_size 5000  -n 1 -cgsp 64  ./bin/qdp_time_avg
