#!/bin/bash
#Note:
#---------------------------
#This script is used for huawei C and Java test.
#Make your job easier.

#----params----
TEST_DIR=/home/hust/runtest
DPDK_APP=/home/hust/zzt/poc-fpga/fp1/software/app/dpdk_app

        echo "----make program----"
        cd $DPDK_APP
        make clean all
        sleep 2
        echo "----copy the *.so!----"
        cp $DPDK_APP/bin/libFPGAPairHMM.so $TEST_DIR/JNILib
        sleep 2
        echo "----Lib copied!----"
        echo "----check your run.sh or run_wes.sh----"
        cd $TEST_DIR

