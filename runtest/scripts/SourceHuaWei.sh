#!/bin/bash

#Note:
#---------------------------
#This script is used for huawei C and Java test.
#Make your job easier.

#----params----
TEST_DIR=/home/hust/runtest
DPDK_APP=/home/hust/zzt/poc-fpga/fp1/software/app/dpdk_app

	echo "----source the huawei script----"
	source $DPDK_APP/build_dpdk_app.sh
	sleep 2
	echo "----source done!----"

