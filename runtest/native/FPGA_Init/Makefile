.PHONY:	build clean
VPATH = .
G++ = g++
build:	liborg_broadinstitute_gatk_engine_FPGA_Init.so
liborg_broadinstitute_gatk_engine_FPGA_Init.so:  \
	org_broadinstitute_gatk_engine_FPGA_Init.o
	${G++} -o $@ -shared $^
	cp liborg_broadinstitute_gatk_engine_FPGA_Init.so ../../JNILib
org_broadinstitute_gatk_engine_FPGA_Init.o:	\
	org_broadinstitute_gatk_engine_FPGA_Init.cc org_broadinstitute_gatk_engine_FPGA_Init.h
	${G++} -c  -fPIC -I ${JAVA_HOME}/include -I ${JAVA_HOME}/include/linux $^
clean:
	rm -f *.o *.so
