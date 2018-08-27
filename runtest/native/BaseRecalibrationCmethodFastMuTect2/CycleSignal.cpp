//
// Created by shawn on 18-3-20.
//
#include <iostream>

#include "CycleSignal.h"

using namespace std;

int CycleSignal::getCycleSignal(int offset, bool NegativeStrandFlag, int readOrderFactor, int readLength) {
    return NegativeStrandFlag ? (((readLength - offset - 1) << 1) + readOrderFactor) : ((offset << 1) + readOrderFactor);
}

int CycleSignal::getIndelCycleSignal(int offset, int readLength, int cycleKey) {
    if(offset < CUSHION_FOR_INDELS || offset >= readLength - CUSHION_FOR_INDELS)
        return -1;
    else
        return (cycleKey << 1) - 16;
}