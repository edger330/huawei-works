//
// Created by shawn on 18-3-20.
//
#include <iostream>
#include "ContextSignal.h"

using namespace std;

int ContextSignal::getMisMatchContextSignal(char *bases, int bases_length, const int offset, const bool NegativeStrandFlag) {
    if (NegativeStrandFlag) {
        if (offset > bases_length - mismatchesContextSize)
            return -1;

        int index = reverseBaseIndexMap(bases[offset + 1]) + (reverseBaseIndexMap(bases[offset]) << 2);//高位是后一个，低位是自己，反转后还是保证高位对应前一个碱基，低位对应后一个碱基
        if(index >= 0)//只要有N碱基存在，index就为负
            return index;
    } else {
        if (offset < mismatchesContextSize - 1)
            return -1;

        int index = baseIndexMap(bases[offset - 1]) + (baseIndexMap(bases[offset]) << 2);//高位是前一个碱基，低位是自己，保证高对应前一个碱基，低位对应后一个碱基
        if(index >= 0)
            return index;
    }
    return -1;
}

int ContextSignal::getIndelContextSignal(char *bases, int bases_length, const int offset, const bool NegativeStrandFlag, int contextSignal) {
    if (NegativeStrandFlag) {
        if (offset > bases_length - indelsContextSize)
            return -1;
        int index = reverseBaseIndexMap(bases[offset + 2]) + (contextSignal << 2);
        if(index >= 0)
            return index;
    } else {
        if (offset < indelsContextSize - 1)
            return -1;
        int index = baseIndexMap(bases[offset - 2]) + (contextSignal << 2);
        if(index >= 0)
            return index;
    }
    return -1;
}

int ContextSignal::reverseBaseIndexMap(char base) {
    switch (base) {
        case 'A':
            return 3;
        case 'C':
            return 2;
        case 'G':
            return 1;
        case 'T':
            return 0;
        default:
            return -100;//N
    }
}

int ContextSignal::baseIndexMap(char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -100;//N
    }
}
