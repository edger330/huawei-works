//
// Created by shawn on 18-3-20.
//

#include "QualitySignal.h"

int QualitySignal::getQualIndex(int qual, vector<int> line) {
    for(int i = (int)(line.size() - 1); i >= 0; i--){//之所以倒着遍历是因为qual为高值的碱基比低值多，倒序遍历比正序遍历的平均次数应该会少
        if(qual == line[i])
            return i;
    }
}