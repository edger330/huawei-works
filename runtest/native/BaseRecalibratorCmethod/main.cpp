#include <iostream>

#include "BAQ.h"
using namespace std;
int main() {
    BAQ::initializeCachedData();
    cout << "Hmm_glocal C代码测试" << endl;
    int l_ref = 82;
    char ref[l_ref] = {84, 65, 67, 65, 67, 65, 84, 65, 65, 67, 84, 65, 67, 84, 67, 67, 65, 84, 65, 84, 71, 71, 65, 65, 84, 65, 67, 84, 71, 71, 71, 71, 65, 71, 71, 65, 71, 71, 84, 71, 84, 84, 67, 67, 65, 65, 65, 84, 65, 65, 65, 71, 65, 71, 65, 67, 84, 71, 65, 71, 71, 65, 84, 84, 84, 67, 84, 67, 65, 84, 71, 65, 71, 65, 65, 67, 84, 67, 65, 71, 84, 71};
    int query_length = 76;
    char query[query_length] = {65, 67, 65, 84, 65, 65, 67, 84, 65, 67, 84, 67, 67, 65, 84, 65, 84, 71, 71, 65, 65, 84, 65, 67, 84, 71, 71, 71, 71, 84, 71, 71, 71, 71, 71, 84, 71, 84, 84, 67, 67, 65, 65, 65, 84, 65, 65, 65, 71, 65, 71, 65, 67, 84, 71, 65, 71, 71, 65, 84, 84, 84, 67, 84, 67, 65, 84, 71, 65, 71, 65, 65, 67, 84, 67, 78};

    int qstart = 0;
    int l_query = query_length;
    char _iqual[l_query] = {32, 32, 32, 32, 32, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 2};
    int state[l_query]= {0};
    char q[l_query] = {0};
    BAQ * baq = new BAQ();
    for(int i = 0; i < 1000000000; i++){
        baq->Chmm_glocal(ref, l_ref, query, query_length, qstart, l_query, _iqual, state, q);
    }
//    baq->Chmm_glocal(ref, l_ref, query, query_length, qstart, l_query, _iqual, state, q);

    for(int i = 0; i < query_length; i++)
        cout << state[i] << ", ";

    cout << endl;

    for(int i = 0; i < query_length; i++)
        cout << (int)q[i] << ", ";
    cout << endl;
    return 0;
}