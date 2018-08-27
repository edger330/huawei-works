//
// Created by shawn on 18-3-20.
//

#ifndef BASERECALIBRATIONCMETHOD_BASERECALIBRATION_H
#define BASERECALIBRATIONCMETHOD_BASERECALIBRATION_H

#include <vector>

#include "QualitySignal.h"
#include "ContextSignal.h"
#include "CycleSignal.h"
#include "MathUtils.h"

using namespace std;

class BaseRecalibration {
public:
    BaseRecalibration(vector<vector<char> > _tables, int _contextMSize,
                      int _cycleMSize, int _max2Cycle, bool _tableIsNull, int _qualTableSize, vector<int> _qual_line);
    BaseRecalibration();
    void recalibrateCmethod(int readLength, bool NegativeStrandFlag, int readOrderFactor,
                            char bases[], char quals[], int rgSignal, char *Iquals, char *Dquals);
    vector<int> getQualLine();

private:
    int preserveQLessThan;
    vector<vector<char> > tables;
    char IempQS;
    char DempQS;
    bool tableIsNull;
    int qualTableSize;
    int contextMSize;
    int cycleMSize;
    int max2Cycle;
    vector<int> qual_line;
};


#endif //BASERECALIBRATIONCMETHOD_BASERECALIBRATION_H
