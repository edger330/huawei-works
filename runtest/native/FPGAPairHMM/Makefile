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



JRE_HOME?=${JAVA_HOME}/jre/
JNI_COMPILATION_FLAGS=-D_REENTRANT -fPIC -I${JRE_HOME}/../include -I${JRE_HOME}/../include/linux

COMMON_COMPILATION_FLAGS=$(JNI_COMPILATION_FLAGS) -O3 -Wall $(OMPCFLAGS) -Wno-unknown-pragmas -Wno-write-strings -Wno-unused-variable -Wno-unused-but-set-variable

C_COMPILER?=gcc
CPP_COMPILER?=g++
AVX_FLAGS=-mavx
SSE41_FLAGS=-msse4.1
COMMON_COMPILATION_FLAGS+=-Wno-char-subscripts

LDFLAGS=-lm -lrt 


BIN=libFPGAPairHMM.so
#BIN=checker

DEPDIR=.deps
DF=$(DEPDIR)/$(*).d

#Common across libJNI and sandbox
#Part of libJNI
LIBSOURCES=org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM.cc
SOURCES=$(LIBSOURCES)
LIBOBJECTS=$(LIBSOURCES:.cc=.o)

OBJECTS=$(LIBOBJECTS)
CXXFLAGS=$(COMMON_COMPILATION_FLAGS)

all: $(BIN) copied_lib

-include $(addprefix $(DEPDIR)/,$(SOURCES:.cc=.d))


libFPGAPairHMM.so: $(LIBOBJECTS) 
	$(CPP_COMPILER)  -shared $(LIBFLAGS) -o $@ $(LIBOBJECTS) ${LDFLAGS}


$(OBJECTS): %.o: %.cc
	@mkdir -p $(DEPDIR)
	$(CPP_COMPILER) -c -MMD -MF $(DF) $(CXXFLAGS) $(OUTPUT_OPTION) $<


copied_lib: libFPGAPairHMM.so
ifdef OUTPUT_DIR
	mkdir -p $(OUTPUT_DIR)
	rsync -a libFPGAPairHMM.so $(OUTPUT_DIR)/
endif

clean:
	rm -rf $(BIN) *.o $(DEPDIR) *.class
