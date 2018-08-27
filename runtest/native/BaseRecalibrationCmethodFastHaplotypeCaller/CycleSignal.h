//
// Created by shawn on 18-3-20.
//

#ifndef BASERECALIBRATIONCMETHOD_CYCLESIGNAL_H
#define BASERECALIBRATIONCMETHOD_CYCLESIGNAL_H

#define CUSHION_FOR_INDELS 4

class CycleSignal {
public:
    static int getCycleSignal(int offset, bool NegativeStrandFlag, int readOrderFactor, int readLength);
    static int getIndelCycleSignal(int offset, int readLength, int cycleKey);
};


#endif //BASERECALIBRATIONCMETHOD_CYCLESIGNAL_H
