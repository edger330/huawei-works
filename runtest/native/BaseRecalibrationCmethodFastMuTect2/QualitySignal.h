//
// Created by shawn on 18-3-20.
//

#ifndef BASERECALIBRATIONCMETHOD_QUALITYSIGNAL_H
#define BASERECALIBRATIONCMETHOD_QUALITYSIGNAL_H

#include <vector>

using namespace std;

class QualitySignal {
public:
    static int getQualIndex(int qual, vector<int> line);
};


#endif //BASERECALIBRATIONCMETHOD_QUALITYSIGNAL_H
