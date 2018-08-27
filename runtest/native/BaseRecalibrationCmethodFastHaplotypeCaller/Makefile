#****************************************************************************
# Copyright       :
# Author          :   Shawn Williams
# Date            :   2017-02-28
# Version         :   JNI 动态链接库
# Description     :   Demo
#****************************************************************************
.PHONY:	build clean
VPATH = .
G++ = g++
build:	libBaseRecalibrationFastHaplotypeCaller.so
libBaseRecalibrationFastHaplotypeCaller.so:  \
    org_broadinstitute_gatk_engine_recalibration_BaseRecalibration.o \
    BaseRecalibration.o \
    ContextSignal.o \
    CycleSignal.o \
    QualitySignal.o \
    MathUtils.o \
    SetTables.o \
    ParseGroupFile.o
	${G++} -o $@ -shared -O3 $^
	cp libBaseRecalibrationFastHaplotypeCaller.so ../../JNILib

org_broadinstitute_gatk_engine_recalibration_BaseRecalibration.o:	org_broadinstitute_gatk_engine_recalibration_BaseRecalibration.cpp	org_broadinstitute_gatk_engine_recalibration_BaseRecalibration.h  BaseRecalibration.h ParseGroupFile.h
#	g++ -c  -fPIC -I /home/shawn/jdk1.8.0_111/include -I /home/shawn/jdk1.8.0_111/include/linux -O3 $^
	${G++} -c  -fPIC -I ${JAVA_HOME}/include -I ${JAVA_HOME}/include/linux -O3 $^
BaseRecalibration.o:	BaseRecalibration.cpp BaseRecalibration.h CycleSignal.h ContextSignal.h SetTables.h
	${G++} -c  -fPIC -O3 $^
ContextSignal.o:	ContextSignal.cpp ContextSignal.h
	${G++} -c  -fPIC -O3 $^
CycleSignal.o:	CycleSignal.cpp CycleSignal.h
	${G++} -c  -fPIC -O3 $^
QualitySignal.o:	QualitySignal.cpp QualitySignal.h
	${G++} -c  -fPIC -O3 $^
MathUtils.o:	MathUtils.cpp MathUtils.h
	${G++} -c  -fPIC -O3 $^
#QualityUtils.o: QualityUtils.cpp QualityUtils.h
#	${G++} -c  -fPIC -O3 $^
# ReadGroupCovariate.o: ReadGroupCovariate.cpp ReadGroupCovariate.h
# 	g++ -c  -fPIC -O3 $^
SetTables.o:	SetTables.cpp SetTables.h
	${G++} -c  -fPIC -O3 $^
ParseGroupFile.o: ParseGroupFile.cpp ParseGroupFile.h
	${G++} -c  -fPIC -O3 $^

clean:
	rm -f *.o *.so