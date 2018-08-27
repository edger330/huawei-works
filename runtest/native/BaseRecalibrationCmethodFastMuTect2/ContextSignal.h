//
// Created by shawn on 18-3-20.
//

#ifndef BASERECALIBRATIONCMETHOD_CONTEXTSIGNAL_H
#define BASERECALIBRATIONCMETHOD_CONTEXTSIGNAL_H


class ContextSignal {
public:
    static int getMisMatchContextSignal(char bases[], int bases_length, const int offset, const bool NegativeStrandFlag);
    static int getIndelContextSignal(char bases[], int bases_length, const int offset, const bool NegativeStrandFlag, int contextSignal);
private:
    static int reverseBaseIndexMap(char base);
    static int baseIndexMap(char base);

    const static int mismatchesContextSize = 2;
    const static int indelsContextSize = 3;
};


#endif //BASERECALIBRATIONCMETHOD_CONTEXTSIGNAL_H
