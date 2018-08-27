//
// Created by shawn on 17-3-7.
//

#include <iostream>
#include <cstdlib>
#include "BAQ.h"
#include <stdio.h>
using namespace std;

BAQ::BAQ() {
    initializeCachedData();
}

void BAQ::Chmm_glocal(const char *ref, int l_ref, const char *query, int query_length, int qstart, int l_query, const char *_iqual, int *state, char *q) {
    if(ref == NULL){
        cout << "BUG: ref sequence is null" << endl;
        exit(1);
    }
    if(query == NULL){
        cout << "BUG: query sequence is null" << endl;
        exit(1);
    }
    if(_iqual == NULL){
        cout << "BUG: query quality vector is null" << endl;
        exit(1);
    }

//    if ( query.length != _iqual.length ) throw new ReviewedGATKException("BUG: read sequence length != qual length");

    if(l_query < 1){
        cout << "BUG: length of query sequence < 0:" << l_query << endl;
        exit(1);
    }
    if(qstart < 0){
        cout << "BUG: query sequence start < 0: " << qstart << endl;
        exit(1);
    }

    int i, k;

    /***initialization***/
//    const int l_ref = ref_length;

    //set band width
//    cout << "l_ref == " << l_ref << "--- l_query == " << l_query << "--- query_length ==" << query_length <<  endl;
    int bw2, bw = l_ref > l_query ? l_ref : l_query;
    if(cb < abs(l_ref - l_query)){
        bw = abs(l_ref - l_query) + 3;
    }
    if(bw > cb) bw = cb;
    if(bw < abs(l_ref - l_query)){
        bw = abs(l_ref - l_query);
    }
    bw2 = bw * 2 + 1;
//    cout << "before new array----" << endl;
//    cout << "m == " << l_query+1 << "--- n == " << bw2*3 + 6 << endl;

//    double f[l_query+1][bw2*3 + 6] = {0};//3*(bw2+2) 3是MDI三种情况，加2是加上start和end的情况
//    double b[l_query+1][bw2*3 + 6] = {0};
//    double s[l_query+2] = {0};
//    cout << "m == " << l_query+1 << "--- n == " << bw2*3 + 6 << endl;

    double ** f = new double *[l_query+1];
    double ** b = new double *[l_query+1];
    double * s = new double [l_query+2];
    for(int i = 0; i < l_query + 1; i++){
        f[i] = new double[bw2*3 + 6];
        b[i] = new double[bw2*3 + 6];
        for(int j = 0; j < bw2*3 + 6; j++){
            f[i][j] = b[i][j] = 0;
        }
        s[i] = 0;
    }
    s[l_query+1] = 0;
    double sM, sI, bM, bI;//sM和sI是从read的最后一个base到end的转移概率；bM是从开始到第一个状态是Match的转移概率，bI是从开始到第一个状态是Insertion的转移概率
    sM = sI = 1. / (2 * l_query + 2);
    bM = (1 - cd) / l_ref; bI = cd / l_ref; // (bM+bI)*l_ref==1

    double m[9];//m 没问题
    m[0*3+0] = (1 - cd - cd) * (1 - sM); m[0*3+1] = m[0*3+2] = cd * (1 - sM);
    m[1*3+0] = (1 - ce) * (1 - sI); m[1*3+1] = ce * (1 - sI); m[1*3+2] = 0.;
    m[2*3+0] = 1 - ce; m[2*3+1] = 0.; m[2*3+2] = ce;

//    cout << "before forward----" << endl;
    /***forward***/
    f[0][set_u(bw, 0, 0)] = s[0] = 1.;//f[0][3] = s[0] = 1， start状态概率？？？
    { // f[1]
        double * fi = f[1];//对fi数组操作就是对f[1]操作，fi数组和f[1]指向内存同一段区域
        double sum;
        int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1, _beg, _end;//beg表示begin
        for (k = beg, sum = 0.; k <= end; ++k) {
            int u;
            double e = calcEpsilon(ref[k-1], query[qstart], _iqual[qstart]);

            u = set_u(bw, 1, k);//set_u转到第一个Match或Insertion可能发生的地方
            fi[u+0] = e * bM; fi[u+1] = EI * bI;
            sum += fi[u] + fi[u+1];
        }
        // rescale
        s[1] = sum;
        _beg = set_u(bw, 1, beg); _end = set_u(bw, 1, end); _end += 2;//其实_end在这里加1就可以
        for (k = _beg; k <= _end; ++k) {
            fi[k] /= sum;//操作目的是让各种情况概率之和为1
        }
    }

    // f[2..l_query]
    for (i = 2; i <= l_query; ++i) {
        double * fi = f[i],  * fi1 = f[i-1];
        double sum;
        int beg = 1, end = l_ref, x, _beg, _end;
        char qyi = query[qstart+i-1];
        x = i - bw; beg = beg > x? beg : x; // band start
        x = i + bw; end = end < x? end : x; // band end
        for (k = beg, sum = 0.; k <= end; ++k) {
            int u, v11, v01, v10;
            double e = calcEpsilon(ref[k-1], qyi, _iqual[qstart+i-1]);//可以把_iqual[qstart+i-1]的值在for循环外先用一个byte变量存放起来，像qyi一样
            u = set_u(bw, i, k); v11 = set_u(bw, i-1, k-1); v10 = set_u(bw, i-1, k); v01 = set_u(bw, i, k-1);
            fi[u+0] = e * (m[0] * fi1[v11+0] + m[3] * fi1[v11+1] + m[6] * fi1[v11+2]);//---fMk(i)=eki*[a00*fMk−1(i−1)+a10*fIk−1(i−1)+a20*fDk−1(i−1)]
            fi[u+1] = EI * (m[1] * fi1[v10+0] + m[4] * fi1[v10+1]);//---fIk(i)
            fi[u+2] = m[2] * fi[v01+0] + m[8] * fi[v01+2];//---fDk(i)
            sum += fi[u] + fi[u+1] + fi[u+2];
            //System.out.println("("+i+","+k+";"+u+"): "+fi[u]+","+fi[u+1]+","+fi[u+2]);
        }
        // rescale
        s[i] = sum;
        _beg = set_u(bw, i, beg); _end = set_u(bw, i, end); _end += 2;
        for (k = _beg, sum = 1./sum; k <= _end; ++k) fi[k] *= sum;
    }
    { // f[l_query+1]
        double sum;
        for (k = 1, sum = 0.; k <= l_ref; ++k) {
            int u = set_u(bw, l_query, k);
            if (u < 3 || u >= bw2*3+3) continue;
            sum += f[l_query][u+0] * sM + f[l_query][u+1] * sI;
        }
        s[l_query+1] = sum; // the last scaling factor
    }
//    cout << "after forward----" << endl;

    /*** backward ***/
    for (k = 1; k <= l_ref; ++k) {
        int u = set_u(bw, l_query, k);
        double * bi = b[l_query];//---可以放for循环外面吧
        if (u < 3 || u >= bw2*3+3) continue;
        bi[u+0] = sM / s[l_query] / s[l_query+1]; bi[u+1] = sI / s[l_query] / s[l_query+1];
    }
    // b[l_query-1..1]
    for (i = l_query - 1; i >= 1; --i) {
        int beg = 1, end = l_ref, x, _beg, _end;
        double * bi = b[i], * bi1 = b[i+1];
        double y = (i > 1)? 1. : 0.;
        char qyi1 = query[qstart+i];
        x = i - bw; beg = beg > x? beg : x;
        x = i + bw; end = end < x? end : x;
        for (k = end; k >= beg; --k) {
            int u, v11, v01, v10;
            u = set_u(bw, i, k); v11 = set_u(bw, i+1, k+1); v10 = set_u(bw, i+1, k); v01 = set_u(bw, i, k+1);
            const double e = (k >= l_ref? 0 : calcEpsilon(ref[k], qyi1, _iqual[qstart+i])) * bi1[v11];
            bi[u+0] = e * m[0] + EI * m[1] * bi1[v10+1] + m[2] * bi[v01+2]; // bi1[v11] has been folded into e.
            bi[u+1] = e * m[3] + EI * m[4] * bi1[v10+1];
            bi[u+2] = (e * m[6] + m[8] * bi[v01+2]) * y;
        }
        // rescale
        _beg = set_u(bw, i, beg); _end = set_u(bw, i, end); _end += 2;
        for (k = _beg, y = 1./s[i]; k <= _end; ++k) bi[k] *= y;
    }

    double pb;
    { // b[0]
        int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1;
        double sum = 0.;
        for (k = end; k >= beg; --k) {
            int u = set_u(bw, 1, k);
            double e = calcEpsilon(ref[k-1], query[qstart], _iqual[qstart]);
            if (u < 3 || u >= bw2*3+3) continue;
            sum += e * b[1][u+0] * bM + EI * b[1][u+1] * bI;
        }
        pb = b[0][set_u(bw, 0, 0)] = sum / s[0]; // if everything works as is expected, pb == 1.0
    }

//    cout << "after backward----" << endl;



    /*** MAP ***/
    for (i = 1; i <= l_query; ++i) {
        double sum = 0., max = 0.;
        const double * fi = f[i], * bi = b[i];
        int beg = 1, end = l_ref, x, max_k = -1;
        x = i - bw; beg = beg > x? beg : x;
        x = i + bw; end = end < x? end : x;
        for (k = beg; k <= end; ++k) {
            const int u = set_u(bw, i, k);//这里初始化u和z可以放在for循环外面吧
            double z;
            sum += (z = (fi[u+0]) * (bi[u+0])); if (z > max) { max = z; max_k = (k-1)<<2 | 0; }//---posterior probability of a read base ci being matching state k̃ (M- or I-typed) is fk̃(i)*bk̃ (i)/p（y）
            sum += (z = (fi[u+1]) * (bi[u+1])); if (z > max) { max = z; max_k = (k-1)<<2 | 1; }
        }
        max /= sum;
//        sum *= s[i]; // if everything works as is expected, sum == 1.0
        if (state != NULL) state[qstart+i-1] = max_k;
        if (q != NULL) {
            k = (int)(-10 * log10(1. - max) + .499); // = 10*log10(1-max) ---乘4.343是因为log的底数是e不是10，加.499是为了四舍五入
            q[qstart+i-1] = (char)(k > 100? 99 : (k < minBaseQual ? minBaseQual : k));//这里稍微有点问题，本来比100大的会变成99，比100小
        }
    }

//    cout << "after map----" << endl;



//    cout << "state--- " << endl;
//    for(int i = 0; i < query_length; i++)
//        cout << state[i] << ", ";
//    cout << endl;
//
//    cout << "q--- " << endl;
//    for(int i = 0; i < query_length; i++)
//        cout << (int)q[i] << ", ";
//    cout << endl;
//    cout << "Chmm_glocal end ----" << endl;
//
//
    for(int i = 0; i < l_query+1; i++){
        delete [] f[i];
        delete [] b[i];
    }
    delete f;
    delete b;
    delete s;
}

int BAQ::set_u(const int b, const int i, const int k) {
    int x = i - b;
    x = x > 0 ? x : 0;
    return (k + 1 - x) * 3;
}

double BAQ::calcEpsilon(char ref, char read, char qualB){
    return EPSILONS[ref][read][(int)qualB];
}

void BAQ::initializeCachedData(){
    for(int i = 0; i < 256; i++){
        qual2prob[i] = pow(10.0, -i/10.);
    }

    for(int i = 0; i < 256; i++)
        for(int j = 0; j < 256; j++)
            for(int q = 0; q <= 93; q++)
                EPSILONS[i][j][q] = 1.0;

    const char * A = "ACGTacgt";
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            for(int q = 0; q <= 93; q++){
                double qual = qual2prob[q < minBaseQual ? minBaseQual : q];
                double e = tolower(A[i]) == tolower(A[j]) ? 1. - qual : qual * EM;
                EPSILONS[(int)A[i]][(int)A[j]][q] = e;
            }
        }
    }
//    for(char b1 : "ACGTacgt"){
//        for(char b2 : "ACGTacgt"){
//            for(int q = 0; q <= 93; q++){
//                double qual = qual2prob[q < minBaseQual ? minBaseQual : q];
//                double e = tolower(b1) == tolower(b2) ? 1. - qual : qual * EM;
//                EPSILONS[(int)b1][(int)b2][q] = e;
//            }
//        }
//    }
}

double BAQ::qual2prob[256] = {0};

double BAQ::cd = pow(10, (-40/10.));
double BAQ::ce = 0.1;
int BAQ::cb = 7;
//    static bool includeClippedBases = false;
char BAQ::minBaseQual = 4;
//    static double DEFAULT_GOP = 40;
double BAQ::EM = 0.33333333333;//emission probability from Match or MisMatch
double BAQ::EI = 0.25;//emission probability from Insertion

double BAQ::EPSILONS[256][256][94] = {0};

void BAQ::findError(double a) {
    if(isnan(a))
        cout << "it's NAN!!!" << endl;
    if(isinf(a))
        cout << "it's INF!!!" << endl;
}
