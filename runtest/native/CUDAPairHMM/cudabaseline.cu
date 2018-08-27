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


#include "common_data_structure.h"
#include "utils.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <omp.h>

using namespace std;

template<class NUMBER>
__device__ NUMBER set_match_to_match_prob(int i, int d, NUMBER * d_jacobianLogTable, NUMBER * d_matchToMatchProb)
{
    NUMBER output;
    NUMBER result;
    int minQual = d;
    int maxQual = i;
    if (i <= d)
    {
        minQual = i;
        maxQual = d;
    }
    NUMBER small = ((NUMBER)-0.1)*minQual;
    NUMBER big =  ((NUMBER)-0.1)*maxQual;
    NUMBER diff = big - small;
    if (isinf(small) == -1 || isinf(big) == -1 || diff >= ((NUMBER)(MAX_JACOBIAN_TOLERANCE))){
        result = big;
    } else {
        NUMBER d = (NUMBER)(diff * ((NUMBER)JACOBIAN_LOG_TABLE_INV_STEP));
        int ind = (d > ((NUMBER)0.0)) ? (int) (d + ((NUMBER)0.5)) : (int) (d - ((NUMBER)0.5));
        result = big + d_jacobianLogTable[ind];
    }
    output = (MAX_QUAL < maxQual) ? ((NUMBER)1.0) - pow(((NUMBER)10), result) : d_matchToMatchProb[((maxQual * (maxQual + 1)) >> 1) + minQual];
    return output;
}

__device__ int find_index(int ROWS, int COLS, int position_x, int position_y)
{
    if (position_x + position_y < ROWS) {
        return (position_x + position_y) * (position_x + position_y + 1) / 2 + position_y;
    } else if (ROWS - 1 - position_x + COLS - 1 - position_y < ROWS) {
        return ROWS * COLS - 1 - (ROWS + COLS - 2 -position_x - position_y) * (ROWS + COLS - 1 - position_x - position_y) / 2 - (COLS - 1 - position_y);
    } else {
        return ROWS * (position_y - 1) - (ROWS - position_x - 1) * (ROWS - 1) + ROWS * (ROWS + 1) / 2;;
    }
}

template <class NUMBER>
__global__ void compute(char * d_tc_rs_all,
                        char * d_tc_hap_all,
                        char * d_tc_q_all,
                        char * d_tc_i_all,
                        char * d_tc_d_all,
                        char * d_tc_c_all,
                        NUMBER * d_ph2pr,
                        NUMBER * d_matchToMatchProb,
                        NUMBER * d_jacobianLogTable,
                        int * d_ROWS, int * d_COLS,
                        int * d_row_offset, int * d_col_offset,
                        NUMBER * d_result,
                        NUMBER INITIAL_CONSTANT)
{
    int index = blockIdx.x;
    int ROWS = d_ROWS[index];
    int COLS = d_COLS[index];
    char * d_tc_rs = d_tc_rs_all + d_row_offset[index];
    char * d_tc_hap = d_tc_hap_all + d_col_offset[index];
    char * d_tc_q = d_tc_q_all + d_row_offset[index];
    char * d_tc_i = d_tc_i_all + d_row_offset[index];
    char * d_tc_d = d_tc_d_all + d_row_offset[index];
    char * d_tc_c = d_tc_c_all + d_row_offset[index];
    int shift = threadIdx.x;

    extern __shared__ __align__(sizeof(NUMBER)) unsigned char temp[];
    NUMBER * common_buffer = reinterpret_cast<NUMBER *>(temp);

    // initialize
    NUMBER * M = common_buffer;
    NUMBER * X = common_buffer + ROWS * COLS;
    NUMBER * Y = X + ROWS * COLS;
    NUMBER * p = Y + ROWS * COLS;

    p[MM] = NUMBER(0.0);
    p[GapM] = NUMBER(0.0);
    p[MX] = NUMBER(0.0);
    p[XX] = NUMBER(0.0);
    p[MY] = NUMBER(0.0);
    p[YY] = NUMBER(0.0);


    if (shift < ROWS) {
        // init p
        int _i = d_tc_i[shift-1] & 127;
        int _d = d_tc_d[shift-1] & 127;
        int _c = d_tc_c[shift-1] & 127;
        p[shift * 6 + MM] = set_match_to_match_prob(_i, _d, d_jacobianLogTable, d_matchToMatchProb);
        p[shift * 6 + GapM] = NUMBER(1.0) - d_ph2pr[_c];
        p[shift * 6 + MX] = d_ph2pr[_i];
        p[shift * 6 + XX] = d_ph2pr[_c];
        p[shift * 6 + MY] = d_ph2pr[_d];
        p[shift * 6 + YY] = d_ph2pr[_c];

        // init row of MXY
        int index = find_index(ROWS, COLS, shift, 0);
        int last_index = find_index(ROWS, COLS, shift - 1, 0);
        M[index] = NUMBER(0.0);
        X[index] = X[last_index] * p[shift * 6 + XX];
        Y[index] = NUMBER(0.0);
    }

    if (shift > COLS) {
        // init col of MXY
        int index = find_index(ROWS, COLS, 0, shift);
        M[index] = NUMBER(0.0);
        X[index] = NUMBER(0.0);
        Y[index] = INITIAL_CONSTANT / (COLS - 1);
    }

    __syncthreads();

    if (shift < ROWS - 1) {
        int pos_x = 1;
        int pos_y = 1;
        int num_thread = 0;
        for (int i = 0; i < ROWS + COLS - 3; i++) {
            if (i < ROWS - 1) // first triangle
            {
                num_thread ++;
                if (i != 0) {
                    pos_x ++;
                }
            } else if (i < COLS - 1) // second area
            {
                pos_y ++;
            } else // third triangle
            {
                num_thread --;
                pos_y ++;
            }
            if (shift < num_thread) {
                int r = pos_x - shift;
                int c = pos_y + shift;
                char _rs = d_tc_rs[r - 1];
                char _hap = d_tc_hap[c - 1];
                int _q = d_tc_q[r - 1] & 127;
                NUMBER distm = d_ph2pr[_q];
                if (_rs == _hap || _rs == 'N' || _hap == 'N')
                    distm = NUMBER(1.0) - distm;
                else
                    distm = distm / 3;

                int position = find_index(ROWS, COLS, r, c);
                int left = find_index(ROWS, COLS, r , c - 1);
                int up = find_index(ROWS, COLS, r - 1, c);
                int leftup = find_index(ROWS, COLS, r - 1, c - 1);

                M[position] = distm * (M[leftup] * p[r * 6 + MM] + X[leftup] * p[r * 6 + GapM] + Y[leftup] * p[r * 6 + GapM]);

                X[position] = M[up] * p[r * 6 + MX] + X[up] * p[r * 6 + XX];

                Y[position] = M[left] * p[r * 6 + MY] + Y[left] * p[r * 6 + YY];
            }
            __syncthreads();
        }
    }

    if (shift == 0) {
        NUMBER result = NUMBER(0.0);
        for (int c = 0; c < COLS; c++) {
            int pos = find_index(ROWS, COLS, ROWS - 1, c);
            result += M[pos] + X[pos];
        }

        d_result[index] = result;
    }
}

template<class NUMBER>
__global__ void initialize(char * d_tc_i_all,
                           char * d_tc_d_all,
                           char * d_tc_c_all,
                           NUMBER * d_M,
                           NUMBER * d_X,
                           NUMBER * d_Y,
                           NUMBER * d_p,
                           NUMBER * d_ph2pr,
                           NUMBER * d_matchToMatchProb,
                           NUMBER * d_jacobianLogTable,
                           int * d_ROWS, int * d_COLS,
                           int * d_row_offset, int * d_col_offset, int * d_MXY_offset, int * d_p_offset,
                           NUMBER INITIAL_CONSTANT)
{
    int index = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    int ROWS = d_ROWS[index];
    int COLS = d_COLS[index];
    char * d_tc_i = d_tc_i_all + d_row_offset[index];
    char * d_tc_d = d_tc_d_all + d_row_offset[index];
    char * d_tc_c = d_tc_c_all + d_row_offset[index];

    NUMBER * M = d_M + d_MXY_offset[index];
    NUMBER * X = d_X + d_MXY_offset[index];
    NUMBER * Y = d_Y + d_MXY_offset[index];
    NUMBER * p = d_p + d_p_offset[index];
    p[MM] = NUMBER(0.0);
    p[GapM] = NUMBER(0.0);
    p[MX] = NUMBER(0.0);
    p[XX] = NUMBER(0.0);
    p[MY] = NUMBER(0.0);
    p[YY] = NUMBER(0.0);

    for (int r = 1; r < ROWS; r++)
    {
        int _i = d_tc_i[r-1] & 127;
        int _d = d_tc_d[r-1] & 127;
        int _c = d_tc_c[r-1] & 127;
        p[r * 6 + MM] = set_match_to_match_prob(_i, _d, d_jacobianLogTable, d_matchToMatchProb);
        p[r * 6 + GapM] = NUMBER(1.0) - d_ph2pr[_c];
        p[r * 6 + MX] = d_ph2pr[_i];
        p[r * 6 + XX] = d_ph2pr[_c];
        p[r * 6 + MY] = d_ph2pr[_d];
        p[r * 6 + YY] = d_ph2pr[_c];
    }

    for (int c = 0; c < COLS; c++)
    {
        int index = find_index(ROWS, COLS, 0, c);
        M[index] = NUMBER(0.0);
        X[index] = NUMBER(0.0);
        Y[index] = INITIAL_CONSTANT / (COLS - 1);
    }

    for (int r = 1; r < ROWS; r++)
    {
        int index = find_index(ROWS, COLS, r, 0);
        int last_index = find_index(ROWS, COLS, r - 1, 0);
        M[index] = NUMBER(0.0);
        X[index] = X[last_index] * p[r * 6 + XX];
        Y[index] = NUMBER(0.0);
    }
}

template<class NUMBER>
__global__ void cuda_compute_batch_prob(char * d_tc_rs_all,
                                        char * d_tc_hap_all,
                                        char * d_tc_q_all,
                                        NUMBER * d_M,
                                        NUMBER * d_X,
                                        NUMBER * d_Y,
                                        NUMBER * d_p,
                                        NUMBER * d_ph2pr,
                                        int * d_ROWS, int * d_COLS,
                                        int * d_row_offset, int * d_col_offset, int * d_MXY_offset, int * d_p_offset)
{
    int index = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    int ROWS = d_ROWS[index];
    int COLS = d_COLS[index];
    char * d_tc_rs = d_tc_rs_all + d_row_offset[index];
    char * d_tc_hap = d_tc_hap_all + d_col_offset[index];
    char * d_tc_q = d_tc_q_all + d_row_offset[index];


    NUMBER * M = d_M + d_MXY_offset[index];
    NUMBER * X = d_X + d_MXY_offset[index];
    NUMBER * Y = d_Y + d_MXY_offset[index];
    NUMBER * p = d_p + d_p_offset[index];

    for (int r = 1; r < ROWS; r++)
        for (int c = 1; c < COLS; c++)
        {
            char _rs = d_tc_rs[r-1];
            char _hap = d_tc_hap[c-1];
            int _q = d_tc_q[r-1] & 127;
            NUMBER distm = d_ph2pr[_q];
            if (_rs == _hap || _rs == 'N' || _hap == 'N')
                distm = NUMBER(1.0) - distm;
            else
                distm = distm/3;

            int position = find_index(ROWS, COLS, r, c);
            int left = find_index(ROWS, COLS, r , c - 1);
            int up = find_index(ROWS, COLS, r - 1, c);
            int leftup = find_index(ROWS, COLS, r - 1, c - 1);

            M[position] = distm * (M[leftup] * p[r * 6 + MM] + X[leftup] * p[r * 6 + GapM] + Y[leftup] * p[r * 6 + GapM]);

            X[position] = M[up] * p[r * 6 + MX] + X[up] * p[r * 6 + XX];

            Y[position] = M[left] * p[r * 6 + MY] + Y[left] * p[r * 6 + YY];

//            M[r * COLS + c] = distm * (M[(r-1) * COLS + (c-1)] * p[r * 6 + MM] + X[(r-1) * COLS + (c-1)] * p[r * 6 + GapM] + Y[(r-1) * COLS + (c-1)] * p[r * 6 + GapM]);
//
//            X[r * COLS + c] = M[(r-1) * COLS + c] * p[r * 6 + MX] + X[(r-1) * COLS + c] * p[r * 6 + XX];
//
//            Y[r * COLS + c] = M[r * COLS + (c-1)] * p[r * 6 + MY] + Y[r * COLS + (c-1)] * p[r * 6 + YY];
        }
}

template<class NUMBER>
__global__ void diagonal_compute(char * d_tc_rs_all,
                                 char * d_tc_hap_all,
                                 char * d_tc_q_all,
                                 NUMBER * d_M,
                                 NUMBER * d_X,
                                 NUMBER * d_Y,
                                 NUMBER * d_p,
                                 NUMBER * d_ph2pr,
                                 int * d_ROWS, int * d_COLS,
                                 int * d_row_offset, int * d_col_offset,
                                 int * d_MXY_offset, int * d_p_offset)
{
    int index = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    int shift = threadIdx.x + threadIdx.y * blockDim.x;
    int ROWS = d_ROWS[index];
    int COLS = d_COLS[index];
    char * d_tc_rs = d_tc_rs_all + d_row_offset[index];
    char * d_tc_hap = d_tc_hap_all + d_col_offset[index];
    char * d_tc_q = d_tc_q_all + d_row_offset[index];

    NUMBER * M = d_M + d_MXY_offset[index];
    NUMBER * X = d_X + d_MXY_offset[index];
    NUMBER * Y = d_Y + d_MXY_offset[index];
    NUMBER * p = d_p + d_p_offset[index];

    extern __shared__ __align__(sizeof(NUMBER)) unsigned char temp[];
    NUMBER * shared_p = reinterpret_cast<NUMBER *>(temp);
    for (int z = threadIdx.x; z < ROWS * 6; z += 32){
        shared_p[z] = p[z];
    }
    __syncthreads();

    if (shift < ROWS - 1) {
        int pos_x = 1;
        int pos_y = 1;
        int num_thread = 0;
        for (int i = 0; i < ROWS + COLS - 3; i++) {
            if (i < ROWS - 1) // first triangle
            {
                num_thread ++;
                if (i != 0) {
                    pos_x ++;
                }
            } else if (i < COLS - 1) // second area
            {
                pos_y ++;
            } else // third triangle
            {
                num_thread --;
                pos_y ++;
            }
            if (shift < num_thread) {
                int r = pos_x - shift;
                int c = pos_y + shift;
                // computing
                char _rs = d_tc_rs[r - 1];
                char _hap = d_tc_hap[c - 1];
                int _q = d_tc_q[r - 1] & 127;
                NUMBER distm = d_ph2pr[_q];
                if (_rs == _hap || _rs == 'N' || _hap == 'N')
                    distm = NUMBER(1.0) - distm;
                else
                    distm = distm / 3;

                int position = find_index(ROWS, COLS, r, c);
                int left = find_index(ROWS, COLS, r , c - 1);
                int up = find_index(ROWS, COLS, r - 1, c);
                int leftup = find_index(ROWS, COLS, r - 1, c - 1);

//                NUMBER tempM = distm * (M[leftup] * p[r * 6 + MM] + X[leftup] * p[r * 6 + GapM] + Y[leftup] * p[r * 6 + GapM]);
                M[position] = distm * (M[leftup] * shared_p[r * 6 + MM] + X[leftup] * shared_p[r * 6 + GapM] + Y[leftup] * shared_p[r * 6 + GapM]);

//                NUMBER tempX= M[up] * p[r * 6 + MX] + X[up] * p[r * 6 + XX];
                X[position] = M[up] * shared_p[r * 6 + MX] + X[up] * shared_p[r * 6 + XX];

//                NUMBER tempY = M[left] * p[r * 6 + MY] + Y[left] * p[r * 6 + YY];
                Y[position] = M[left] * shared_p[r * 6 + MY] + Y[left] * shared_p[r * 6 + YY];
            }
            __syncthreads();
        }
    }
}

template<class NUMBER>
__global__ void compute_result(NUMBER * d_M,
                               NUMBER * d_X,
                               int * d_ROWS, int * d_COLS,
                               int * d_MXY_offset,
                               NUMBER * d_result)
{
    int index = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    int ROWS = d_ROWS[index];
    int COLS = d_COLS[index];

    NUMBER * M = d_M + d_MXY_offset[index];
    NUMBER * X = d_X + d_MXY_offset[index];

    NUMBER result = NUMBER(0.0);
    for (int c = 0; c < COLS; c++)
    {
        int pos = find_index(ROWS, COLS, ROWS - 1, c);
        result += M[pos] + X[pos];
    }

    d_result[index] = result;
}

resultSet * cuda_compute_full_prob(dataSet * data) {

//    std::cout << "1" << std::endl;
//    std::cout << "size : " << data->numHaps * data->numReads << std::endl;
    Context<float > ctx;
    resultSet * ret = new resultSet();
    ret->taskID = data->taskID;
    ret->sampleID = data->sampleID;
    ret->numHaps = data->numHaps;
    ret->numReads = data->numReads;
    ret->result = (float *)malloc(ret->numHaps * ret->numReads * sizeof(float));
    // rs q i d c
    int sizeRead = 4096 / ret->numHaps;
    int size = sizeRead * ret->numHaps;
//    std::cout << "size : " << size << std::endl;
//    std::cout << "size read : " << sizeRead << std::endl;

//    std::cout << "2" << std::endl;

    int readsDataLens = 0;
    int hapsDataLens = 0;
    int * readsOffset = (int *)malloc(ret->numReads * sizeof(int));
    int * hapsOffset = (int *)malloc(ret->numHaps * sizeof(int));
    for (int i = 0; i < ret->numReads; i++) {
        readsOffset[i] = readsDataLens;
        readsDataLens += data->readsLens[i];
//        std::cout << "readsLens[] : " << data->readsLens[i] << std::endl;
    }

    for (int i = 0; i < ret->numHaps; i++) {
        hapsOffset[i] = hapsDataLens;
        hapsDataLens += data->hapsLens[i];
//        std::cout << "hapsLens[] : " << data->hapsLens[i] << std::endl;
    }

    int result_offset = 0;
//    std::cout << "3" << std::endl;
    for (int batch_index = 0 ; batch_index < ret->numReads; batch_index += sizeRead) {
        int Vec_size;
        if (batch_index + sizeRead < ret->numReads) {
//            std::cout << "enter 1" << std::endl;
            Vec_size = size;
        } else {
//            std::cout << "enter 2" << std::endl;
            Vec_size = (ret->numReads - batch_index) * ret->numHaps;
        }
        int readSize = Vec_size / ret->numHaps;
//        std::cout << "readSize : " << readSize << std::endl;

        int *d_ROWS;
        int *d_COLS;
        int *ROWS_all = (int *) malloc(Vec_size * sizeof(int));
        int *COLS_all = (int *) malloc(Vec_size * sizeof(int));
        int *d_row_offset;
        int *d_col_offset;
        int *row_offset = (int *) malloc(Vec_size * sizeof(int));
        int *col_offset = (int *) malloc(Vec_size * sizeof(int));
        int *d_MXY_offset;
        int *d_p_offset;
        int *MXY_offset = (int *) malloc(Vec_size * sizeof(int));
        int *p_offset = (int *) malloc(Vec_size * sizeof(int));
        int max_thread = 0;
        int max_cols = 0;

        int R_off = 0;
        int C_off = 0;
        int mxy_off = 0;
        int p_off = 0;
        for (int read_index = 0; read_index < readSize; read_index++) {
            for (int hap_index = 0; hap_index < ret->numHaps; hap_index++) {
                int ROWS = data->readsLens[read_index + batch_index] + 1;
                int COLS = data->hapsLens[hap_index] + 1;
                row_offset[read_index * ret->numHaps + hap_index] = R_off;
                col_offset[read_index * ret->numHaps + hap_index] = C_off;
                MXY_offset[read_index * ret->numHaps + hap_index] = mxy_off;
                p_offset[read_index * ret->numHaps + hap_index] = p_off;
                R_off += (ROWS - 1);
                C_off += (COLS - 1);
                mxy_off += (ROWS * COLS);
                p_off += (ROWS * 6);
                ROWS_all[read_index * ret->numHaps + hap_index] = ROWS;
                COLS_all[read_index * ret->numHaps + hap_index] = COLS;
                if (max_thread < ROWS - 1) {
                    max_thread = ROWS - 1;
                }
                if (max_cols < COLS) {
                    max_cols = COLS;
                }
            }
        }
//        std::cout << "4" << std::endl;
        char *d_tc_rs_all;
        char *d_tc_hap_all;
        char *d_tc_q_all;
        char *d_tc_i_all;
        char *d_tc_d_all;
        char *d_tc_c_all;
        float * d_M;
        float * d_X;
        float * d_Y;
        float * d_p;
        char *tc_rs_all = (char *)malloc(R_off * sizeof(char));
        char *tc_hap_all = (char *)malloc(C_off * sizeof(char));
        char *tc_q_all = (char *)malloc(R_off * sizeof(char));
        char *tc_i_all = (char *)malloc(R_off * sizeof(char));
        char *tc_d_all = (char *)malloc(R_off * sizeof(char));
        char *tc_c_all = (char *)malloc(R_off * sizeof(char));
//        std::cout << "R_OFF : " << R_off << std::endl;
//        std::cout << "C_OFF : " << C_off << std::endl;
//        std::cout << "row_offset[]" << row_offset[(readSize - 1) * ret->numHaps + ret->numHaps -1] << std::endl;
//        std::cout << "col_offset[]" << col_offset[(readSize - 1) * ret->numHaps + ret->numHaps -1] << std::endl;
//        std::cout << "enter memcpy" << std::endl;
        for (int read_index = 0; read_index < readSize; read_index++) {
            for (int hap_index = 0; hap_index < ret->numHaps; hap_index++) {
//                std::cout << "ROWS : " << ROWS_all[read_index * ret->numHaps + hap_index] - 1 << std::endl;
                memcpy(tc_rs_all + row_offset[read_index * ret->numHaps + hap_index], data->readsData + readsOffset[read_index + batch_index], (ROWS_all[read_index * ret->numHaps + hap_index] - 1) * sizeof(char));
//                std::cout << "row_offset : " << row_offset[read_index * ret->numHaps + hap_index] << std::endl;
                memcpy(tc_hap_all + col_offset[read_index * ret->numHaps + hap_index], data->hapsData + hapsOffset[hap_index], (COLS_all[read_index * ret->numHaps + hap_index] - 1) * sizeof(char));
//                std::cout << "read_offset : " << readsOffset[read_index + batch_index] << std::endl;
                memcpy(tc_q_all + row_offset[read_index * ret->numHaps + hap_index], data->readsData + readsOffset[read_index + batch_index] + readsDataLens, (ROWS_all[read_index * ret->numHaps + hap_index] - 1) * sizeof(char));
//                std::cout << "readsDataLen : " << readsDataLens << std::endl;
                memcpy(tc_i_all + row_offset[read_index * ret->numHaps + hap_index], data->readsData + readsOffset[read_index + batch_index] + 2 * readsDataLens, (ROWS_all[read_index * ret->numHaps + hap_index] - 1) * sizeof(char));
                memcpy(tc_d_all + row_offset[read_index * ret->numHaps + hap_index], data->readsData + readsOffset[read_index + batch_index] + 3 * readsDataLens, (ROWS_all[read_index * ret->numHaps + hap_index] - 1) * sizeof(char));
                memcpy(tc_c_all + row_offset[read_index * ret->numHaps + hap_index], data->readsData + readsOffset[read_index + batch_index] + 4 * readsDataLens, (ROWS_all[read_index * ret->numHaps + hap_index] - 1) * sizeof(char));
            }
        }
//        std::cout << "5" << std::endl;
        cudaMalloc((void **) &d_M, mxy_off * sizeof(float));
        cudaMalloc((void **) &d_X, mxy_off * sizeof(float));
        cudaMalloc((void **) &d_Y, mxy_off * sizeof(float));
        cudaMalloc((void **) &d_p, p_off * sizeof(float));
        cudaMalloc((void **) &d_tc_rs_all, R_off * sizeof(char));
        cudaMalloc((void **) &d_tc_hap_all, C_off * sizeof(char));
        cudaMalloc((void **) &d_tc_q_all, R_off * sizeof(char));
        cudaMalloc((void **) &d_tc_i_all, R_off * sizeof(char));
        cudaMalloc((void **) &d_tc_d_all, R_off * sizeof(char));
        cudaMalloc((void **) &d_tc_c_all, R_off * sizeof(char));
        cudaMemcpy(d_tc_rs_all, tc_rs_all, R_off * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_tc_hap_all, tc_hap_all, C_off * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_tc_q_all, tc_q_all, R_off * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_tc_i_all, tc_i_all, R_off * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_tc_d_all, tc_d_all, R_off * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_tc_c_all, tc_c_all, R_off * sizeof(char), cudaMemcpyHostToDevice);

        cudaMalloc((void **) &d_ROWS, Vec_size * sizeof(int));
        cudaMalloc((void **) &d_COLS, Vec_size * sizeof(int));
        cudaMemcpy(d_ROWS, ROWS_all, Vec_size * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_COLS, COLS_all, Vec_size * sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc((void **) &d_row_offset, Vec_size * sizeof(int));
        cudaMalloc((void **) &d_col_offset, Vec_size * sizeof(int));
        cudaMemcpy(d_row_offset, row_offset, Vec_size * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_col_offset, col_offset, Vec_size * sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc((void **) &d_MXY_offset, Vec_size * sizeof(int));
        cudaMalloc((void **) &d_p_offset, Vec_size * sizeof(int));
        cudaMemcpy(d_MXY_offset, MXY_offset, Vec_size * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_p_offset, p_offset, Vec_size * sizeof(int), cudaMemcpyHostToDevice);

        float * d_ph2pr;
        float * d_matchToMatchProb;

        float * d_jacobianLogTable;
        float * d_result;
        float * result = (float *)malloc(Vec_size * sizeof(float));

        cudaMalloc((void **) &d_ph2pr, 128 * sizeof(float));
        cudaMalloc((void **) &d_matchToMatchProb, (((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1) * sizeof(float));
        cudaMalloc((void **) &d_jacobianLogTable, JACOBIAN_LOG_TABLE_SIZE * sizeof(float));
        cudaMalloc((void **) &d_result, Vec_size * sizeof(float));

        cudaMemcpy(d_ph2pr, ctx.ph2pr, 128 * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_matchToMatchProb, ctx.matchToMatchProb, (((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1) * sizeof(float),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_jacobianLogTable, ctx.jacobianLogTable, JACOBIAN_LOG_TABLE_SIZE * sizeof(float),
                   cudaMemcpyHostToDevice);
//        std::cout << "6" << std::endl;
        initialize << < Vec_size, 1 >> >
                                  (d_tc_i_all, d_tc_d_all, d_tc_c_all, d_M, d_X, d_Y, d_p, d_ph2pr, d_matchToMatchProb, d_jacobianLogTable, d_ROWS, d_COLS, d_row_offset, d_col_offset, d_MXY_offset, d_p_offset, ctx.INITIAL_CONSTANT);
        cudaThreadSynchronize();

        int warps = max_thread / 32;
        if (warps * 32 < max_thread) {
            warps++;
        }
        dim3 block(32, warps);
        diagonal_compute << < Vec_size, block, (max_thread + 1) * 6 * sizeof(float) >> >
                                               (d_tc_rs_all, d_tc_hap_all, d_tc_q_all, d_M, d_X, d_Y, d_p, d_ph2pr, d_ROWS, d_COLS, d_row_offset, d_col_offset, d_MXY_offset, d_p_offset);
        cudaThreadSynchronize();

        compute_result << < Vec_size, 1 >> > (d_M, d_X, d_ROWS, d_COLS, d_MXY_offset, d_result);
        cudaThreadSynchronize();

        cudaMemcpy(result, d_result, Vec_size * sizeof(float), cudaMemcpyDeviceToHost);
//        std::cout << "7" << std::endl;
        cudaFree(d_tc_rs_all);
        cudaFree(d_tc_hap_all);
        cudaFree(d_tc_q_all);
        cudaFree(d_tc_i_all);
        cudaFree(d_tc_d_all);
        cudaFree(d_tc_c_all);
//        std::cout << "finish device qidc" << std::endl;
        cudaFree(d_M);
        cudaFree(d_X);
        cudaFree(d_Y);
        cudaFree(d_p);
//        std::cout << "finish MXYp" << std::endl;
        free(tc_rs_all);
        free(tc_hap_all);
        free(tc_q_all);
        free(tc_i_all);
        free(tc_d_all);
        free(tc_c_all);
//        std::cout << "finish host qidc" << std::endl;
        cudaFree(d_ph2pr);
        cudaFree(d_matchToMatchProb);
        cudaFree(d_jacobianLogTable);
        cudaFree(d_result);
        cudaFree(d_ROWS);
        cudaFree(d_COLS);
//        std::cout << "finish device ph2COLS" << std::endl;
        free(ROWS_all);
        free(COLS_all);
//        std::cout << "finish host ROW COL" << std::endl;
        cudaFree(d_row_offset);
        cudaFree(d_col_offset);
//        std::cout << "finish device offset" << std::endl;
        free(row_offset);
        free(col_offset);
//        std::cout << "finish host offset" << std::endl;
        cudaFree(d_MXY_offset);
        cudaFree(d_p_offset);
//        std::cout << "finish deivce MXYp offset" << std::endl;
        free(MXY_offset);
        free(p_offset);
//        std::cout << "finish host MXYp offset" << std::endl;
//        std::cout << "result_offset : " << result_offset << std::endl;
        memcpy(ret->result + result_offset, result, Vec_size * sizeof(float));
        free(result);
        result_offset += Vec_size;
//        std::cout << "8" << std::endl;
    }
//    std::cout << "9" << std::endl;
    free(readsOffset);
    free(hapsOffset);
    return ret;
}

//template double * cuda_compute_full_prob<double>(vector<testcase>& tcVec);
//template float * cuda_compute_full_prob<float>(vector<testcase>& tcVec);

