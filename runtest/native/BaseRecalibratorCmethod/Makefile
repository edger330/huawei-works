#****************************************************************************
# Copyright       :
# Author          :   Shawn Williams
# Date            :   2017-03-08
# Version         :   JNI 动态链接库
# Description     :   Demo
#****************************************************************************
.PHONY:	build clean
VPATH = .
G++ = g++
build:	liborg_broadinstitute_gatk_utils_baq_BAQ.so
liborg_broadinstitute_gatk_utils_baq_BAQ.so:  \
	org_broadinstitute_gatk_utils_baq_BAQ.o BAQ.o
	${G++} -o $@ -shared -O3 $^
	cp liborg_broadinstitute_gatk_utils_baq_BAQ.so ../../JNILib
org_broadinstitute_gatk_utils_baq_BAQ.o:	\
	org_broadinstitute_gatk_utils_baq_BAQ.cpp org_broadinstitute_gatk_utils_baq_BAQ.h  BAQ.h
	${G++} -c  -fPIC -I ${JAVA_HOME}/include -I ${JAVA_HOME}/include/linux -O3 $^
BAQ.o:	\
	BAQ.cpp BAQ.h
	${G++} -c -fPIC -O3 $^
clean:
	rm -f *.o *.so
#	rm -f *.so