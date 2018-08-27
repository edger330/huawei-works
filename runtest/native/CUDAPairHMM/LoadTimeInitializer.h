/*Copyright (c) 2012 The Broad Institute

*Permission is hereby granted, free of charge, to any person
*obtaining a copy of this software and associated documentation
*files (the "Software"), to deal in the Software without
*restriction, including without limitation the rights to use,
*copy, modify, merge, publish, distribute, sublicense, and/or sell
*copies of the Software, and to permit persons to whom the
*Software is furnished to do so, subject to the following
*conditions:

*The above copyright notice and this permission notice shall be
*included in all copies or substantial portions of the Software.

*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
*EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
*OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
*NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
*HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
*WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
*FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
*THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef LOAD_TIME_INITIALIZER_H
#define LOAD_TIME_INITIALIZER_H
#include "headers.h"
#include <jni.h>
/*#include "template.h"*/

enum LoadTimeInitializerStatsEnum
{
  NUM_REGIONS_IDX=0,
  NUM_READS_IDX,
  NUM_HAPLOTYPES_IDX,
  NUM_TESTCASES_IDX,
  NUM_DOUBLE_INVOCATIONS_IDX,
  HAPLOTYPE_LENGTH_IDX,
  READ_LENGTH_IDX,
  PRODUCT_READ_LENGTH_HAPLOTYPE_LENGTH_IDX,
  TOTAL_NUMBER_STATS
};
extern char* LoadTimeInitializerStatsNames[];

class LoadTimeInitializer
{
  public:
    LoadTimeInitializer(){};		//will be called when library is loaded
    
    jfieldID m_readBasesFID;
    jfieldID m_readQualsFID;
    jfieldID m_insertionGOPFID;
    jfieldID m_deletionGOPFID;
    jfieldID m_overallGCPFID;
    jfieldID m_haplotypeBasesFID;
    //profiling - update stats
};
extern LoadTimeInitializer g_load_time_initializer;


#endif
