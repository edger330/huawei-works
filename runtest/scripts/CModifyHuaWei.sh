#!/bin/bash
#Note:
#---------------------------
#This script is used for huawei C and Java test.
#Make your job easier.

#----params----
TEST_DIR=/home/hust/runtest
DPDK_APP=/home/hust/zzt/poc-fpga/fp1/software/app/dpdk_app

        echo "----modify the C code----"
        cd $DPDK_APP/example3/test
        echo "----check org_*.c and queue_*.c!----"
	sleep 2
	vim org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM.c

