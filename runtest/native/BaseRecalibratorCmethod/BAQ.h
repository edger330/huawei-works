//
// Created by shawn on 17-3-7.
//

#ifndef BASERECALIBRATORCMETHOD_BAQ_H
#define BASERECALIBRATORCMETHOD_BAQ_H

//#include "Math.h"
#include <string>
#include <cmath>

class BAQ {
public:
    void Chmm_glocal(const char *ref, int l_ref, const char *query, int query_length, int qstart, int l_query, const char *_iqual, int *state, char *q);
    static void initializeCachedData();
    BAQ();
private:
    static int set_u(const int b, const int i, const int k);
    static double cd;
    static double ce;
    static int cb;
//    static bool includeClippedBases = false;
    static char minBaseQual;
//    static double DEFAULT_GOP = 40;
    static double EM;//emission probability from Match or MisMatch
    static double EI;//emission probability from Insertion
    static double EPSILONS[256][256][94];
    static double calcEpsilon(char ref, char read, char qualB);
    static double qual2prob[256];

    static void findError(double a);
};


#endif //BASERECALIBRATORCMETHOD_BAQ_H
