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

#include <stdlib.h>
#include "headers.h"
#include "common_data_structure.h"
#include "utils.h"
#include "LoadTimeInitializer.h"
using namespace std;


template<class NUMBER>
inline void writeMatrix(FILE *fp, NUMBER ** mm,int rows,int cols){
  for (int r = 0; r < rows; r++){
    for (int c = 0; c < cols; c++){
//        fprintf(fp,"%x ", *((unsigned int *)&mm[r][c]));
          fprintf(fp,"%f ", mm[r][c]);
    }
    fprintf(fp,"\n");
  }
}


template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER *before_last_log)
{

  int r, c;
  int ROWS = tc->rslen + 1;
  int COLS = tc->haplen + 1;


  Context<NUMBER> ctx;
  //#define USE_STACK_ALLOCATION 1
#ifdef USE_STACK_ALLOCATION
  NUMBER M[ROWS][COLS];
  NUMBER X[ROWS][COLS];
  NUMBER Y[ROWS][COLS];
  NUMBER p[ROWS][6];
#else
  //allocate on heap in way that simulates a 2D array. Having a 2D array instead of
  //a straightforward array of pointers ensures that all data lies 'close' in memory, increasing
  //the chance of being stored together in the cache. Also, prefetchers can learn memory access
  //patterns for 2D arrays, not possible for array of pointers
  //NUMBER* common_buffer = 0;
  NUMBER* common_buffer = new NUMBER[3*ROWS*COLS + ROWS*6];
  //pointers to within the allocated buffer
  NUMBER** common_pointer_buffer = new NUMBER*[4*ROWS];
  NUMBER* ptr = common_buffer;
  unsigned i = 0;
  for(i=0;i<3*ROWS;++i, ptr+=COLS)
    common_pointer_buffer[i] = ptr;
  for(;i<4*ROWS;++i, ptr+=6)
    common_pointer_buffer[i] = ptr;

  NUMBER** M = common_pointer_buffer;
  NUMBER** X = M + ROWS;
  NUMBER** Y = X + ROWS;
  NUMBER** p = Y + ROWS;
#endif


  p[0][MM] = ctx._(0.0);
  p[0][GapM] = ctx._(0.0);
  p[0][MX] = ctx._(0.0);
  p[0][XX] = ctx._(0.0);
  p[0][MY] = ctx._(0.0);
  p[0][YY] = ctx._(0.0);

  for (r = 1; r < ROWS; r++)
  {
    int _i = tc->i[r-1] & 127;
    int _d = tc->d[r-1] & 127;
    int _c = tc->c[r-1] & 127;
    //p[r][MM] = ctx._(1.0) - ctx.ph2pr[(_i + _d) & 127];
    SET_MATCH_TO_MATCH_PROB(p[r][MM], _i, _d);
    p[r][GapM] = ctx._(1.0) - ctx.ph2pr[_c];
    p[r][MX] = ctx.ph2pr[_i];
    p[r][XX] = ctx.ph2pr[_c];
    p[r][MY] = ctx.ph2pr[_d];
    p[r][YY] = ctx.ph2pr[_c];
    //p[r][MY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_d];
    //p[r][YY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_c];
  }
  for (c = 0; c < COLS; c++)
  {
    M[0][c] = ctx._(0.0);
    X[0][c] = ctx._(0.0);
    Y[0][c] = ctx.INITIAL_CONSTANT / (tc->haplen);
  }

  for (r = 1; r < ROWS; r++)
  {
    M[r][0] = ctx._(0.0);
    X[r][0] = X[r-1][0] * p[r][XX];
    Y[r][0] = ctx._(0.0);
  }

  NUMBER result = ctx._(0.0);

  for (r = 1; r < ROWS; r++)
    for (c = 1; c < COLS; c++)
    {
      fexcept_t flagp;
      char _rs = tc->rs[r-1];
      char _hap = tc->hap[c-1];
      int _q = tc->q[r-1] & 127;
      NUMBER distm = ctx.ph2pr[_q];
      if (_rs == _hap || _rs == 'N' || _hap == 'N')
        distm = ctx._(1.0) - distm;
      else
        distm = distm/3;


      //feclearexcept(FE_ALL_EXCEPT);
      M[r][c] = distm * (M[r-1][c-1] * p[r][MM] + X[r-1][c-1] * p[r][GapM] + Y[r-1][c-1] * p[r][GapM]);
      //STORE_FP_EXCEPTIONS(flagp, exceptions_array);

      //feclearexcept(FE_ALL_EXCEPT);
      X[r][c] = M[r-1][c] * p[r][MX] + X[r-1][c] * p[r][XX];
      //STORE_FP_EXCEPTIONS(flagp, exceptions_array);

      //feclearexcept(FE_ALL_EXCEPT);
      Y[r][c] = M[r][c-1] * p[r][MY] + Y[r][c-1] * p[r][YY];
      //STORE_FP_EXCEPTIONS(flagp, exceptions_array);

      //CONVERT_AND_PRINT(M[r][c]);
      //CONVERT_AND_PRINT(X[r][c]);
      //CONVERT_AND_PRINT(Y[r][c]);

    }


  for (c = 0; c < COLS; c++)
  {
    result += M[ROWS-1][c] + X[ROWS-1][c];
  }

  if (before_last_log != NULL)
    *before_last_log = result;

#ifndef USE_STACK_ALLOCATION
  delete[] common_pointer_buffer;
  delete[] common_buffer;
#endif
  return result;
  //return ctx.LOG10(result) - ctx.LOG10_INITIAL_CONSTANT;
}

template double compute_full_prob<double>(testcase* tc, double* nextbuf);
template float compute_full_prob<float>(testcase* tc, float* nextbuf);

