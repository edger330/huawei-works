#Copyright (c) 2012 The Broad Institute

#Permission is hereby granted, free of charge, to any person
#obtaining a copy of this software and associated documentation
#files (the "Software"), to deal in the Software without
#restriction, including without limitation the rights to use,
#copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the
#Software is furnished to do so, subject to the following
#conditions:

#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
#THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


OMPCFLAGS=-fopenmp
OMPLFLAGS=-fopenmp 
#-openmp-link static

#CFLAGS=-O2 -std=c++11 -W -Wall -march=corei7-avx -Wa,-q            -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas
#CFLAGS=-O2             -W -Wall -march=corei7 -mfpmath=sse -msse4.2 -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas

JRE_HOME?=${JAVA_HOME}/jre
JNI_COMPILATION_FLAGS=-D_REENTRANT -fPIC -I${JRE_HOME}/../include -I${JRE_HOME}/../include/linux

#COMMON_COMPILATION_FLAGS=$(JNI_COMPILATION_FLAGS) -O3 -W -Wall -pedantic $(OMPCFLAGS) -Wno-unknown-pragmas
COMMON_COMPILATION_FLAGS=$(JNI_COMPILATION_FLAGS) -O3 -Wall $(OMPCFLAGS) -Wno-unknown-pragmas -Wno-write-strings -Wno-unused-variable -Wno-unused-but-set-variable
ifdef DISABLE_FTZ
  COMMON_COMPILATION_FLAGS+=-DDISABLE_FTZ
endif

#define USE_GCC
#endef

#ifdef USE_GCC
  NVCC=nvcc
  C_COMPILER?=gcc
  CPP_COMPILER?=g++
  AVX_FLAGS=-mavx
  SSE41_FLAGS=-msse4.1
  COMMON_COMPILATION_FLAGS+=-Wno-char-subscripts
#else
#  C_COMPILER?=icc
#  CPP_COMPILER?=icc
#  AVX_FLAGS=-xAVX
#  SSE41_FLAGS=-xSSE4.1
#  LIBFLAGS=-static-intel
#  ifdef DISABLE_FTZ
#    COMMON_COMPILATION_FLAGS+=-no-ftz
#  endif
#endif

LDFLAGS=-lm -lrt -lgcc_s -lgcc $(OMPLDFLAGS)
CUDALDFLAGS=-L /usr/local/cuda/lib64 -lcuda -lcudart

PAPI_DIR=/home/karthikg/softwares/papi-5.3.0
ifdef USE_PAPI
  ifeq ($(USE_PAPI),1)
    COMMON_COMPILATION_FLAGS+=-I$(PAPI_DIR)/include -DUSE_PAPI
    LDFLAGS+=-L$(PAPI_DIR)/lib -lpapi
  endif
endif

#BIN=libVectorLoglessPairHMM.so checker pairhmm-template-main 
BIN=libCUDAPairHMM.so
#BIN=checker

DEPDIR=.deps
DF=$(DEPDIR)/$(*).d

#Common across libJNI and sandbox
COMMON_SOURCES=utils.cc baseline.cc LoadTimeInitializer.cc
#Part of libJNI
LIBSOURCES=org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM.cc $(COMMON_SOURCES)
CUDA_SOURCES=cudabaseline.cu
CUDA_OBJECTS=cudabaseline.o

SOURCES=$(LIBSOURCES)
LIBOBJECTS=$(LIBSOURCES:.cc=.o) $(CUDA_OBJECTS)
test:
	echo ${LIBOBJECTS}
COMMON_OBJECTS=$(COMMON_SOURCES:.cc=.o)


#No vectorization for these files
NO_VECTOR_SOURCES=utils.cc baseline.cc LoadTimeInitializer.cc org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM.cc

NO_VECTOR_OBJECTS=$(NO_VECTOR_SOURCES:.cc=.o)
$(NO_VECTOR_OBJECTS): CXXFLAGS=$(COMMON_COMPILATION_FLAGS)
OBJECTS=$(NO_VECTOR_OBJECTS)

all: $(BIN)

-include $(addprefix $(DEPDIR)/,$(SOURCES:.cc=.d))

checker: pairhmm-1-base.o $(COMMON_OBJECTS) cudabaseline.o
	$(CPP_COMPILER) $(OMPLFLAGS) -o $@ $^ $(LDFLAGS) ${CUDALDFLAGS}

pairhmm-template-main:	pairhmm-template-main.o $(COMMON_OBJECTS) cudabaseline.o
	$(CPP_COMPILER) $(OMPLFLAGS) -o $@ $^ $(LDFLAGS) ${CUDALDFLAGS}

libCUDAPairHMM.so: $(LIBOBJECTS) 
	$(CPP_COMPILER) $(OMPLFLAGS) -shared $(LIBFLAGS) -o $@ $(LIBOBJECTS) -I ./common_data_structure.h ${LDFLAGS} ${CUDALDFLAGS}


$(OBJECTS): %.o: %.cc
	@mkdir -p $(DEPDIR)
	$(CPP_COMPILER) -c -g -fopenmp -MMD -MF $(DF) $(CXXFLAGS) $(OUTPUT_OPTION) $<

$(CUDA_OBJECTS): %.o: %.cu
	$(NVCC) -c -O3 -g --compiler-options "-Wall -Wfatal-errors -Ofast -DOPENCV -DGPU -DCUDNN -fPIC" $<


clean:
	rm -rf $(BIN) *.o $(DEPDIR) *.class
