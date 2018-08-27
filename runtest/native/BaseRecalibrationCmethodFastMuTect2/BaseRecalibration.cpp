//
// Created by shawn on 18-3-20.
//

#include <cstring>
#include <iostream>

#include "BaseRecalibration.h"
//#include "QualitySignal.h"
//#include "ContextSignal.h"
//#include "CycleSignal.h"
//#include "MathUtils.h"

#define MAX_RECALIBRATED_Q_SCORE 93

BaseRecalibration::BaseRecalibration(vector<vector<char> > _tables, int _contextMSize,
                                     int _cycleMSize, int _max2Cycle, bool _tableIsNull, int _qualTableSize, vector<int> _qual_line){
    contextMSize = _contextMSize;
    cycleMSize = _cycleMSize;
    max2Cycle = _max2Cycle;
    tableIsNull = _tableIsNull;
    qualTableSize = _qualTableSize;
    qual_line = _qual_line;
    if(!tableIsNull){
        preserveQLessThan = 6;
//    disableIndelQuals = false;
        tables = _tables;
        IempQS = tables[1][qualTableSize - 2];
        DempQS = tables[1][qualTableSize - 1];
    }

//    cout << "contextMSize = " << contextMSize << endl;
//    cout << "cycleMSize = " << cycleMSize << endl;
//    cout << "max2Cycle = " << max2Cycle << endl;
//    cout << "tableIsNull = " << tableIsNull << endl;
//    cout << "qualTableSize = " << qualTableSize << endl;
//    cout << "qual_line.size()  = " << qual_line.size() << endl;
//    for(int i = 0; i < qual_line.size(); i++){
//        cout << qual_line[i] << endl;
//    }
//    cout << "///////////" << endl;

}

void BaseRecalibration::recalibrateCmethod(int readLength, bool NegativeStrandFlag, int readOrderFactor,
                                           char *bases, char *quals, int rgSignal, char *Iquals, char *Dquals) {
//    cout << "In recalibrateCmethod tableIsNull == " << tableIsNull << endl;
    if(tableIsNull){
        memset(Iquals, 45, readLength);
        memset(Dquals, 45, readLength);
        return;
    }
    char origQual;
    int qualKey, cycleKey, contextKey;
    int IcycleKey, IcontextKey;
    char empQS, deltaQCovariates, IdeltaQCovariates, DdeltaQCovariates;

    for(int offset = 0; offset < readLength; offset++){
//        if (offset == 19)
//            int k = 0;
//        cout << offset << "///////////////////////////////////" << endl;
        origQual = quals[offset];
//        cout << "origQual = " << (int)origQual << endl;
//        if(origQual >= preserverQLessThan){//这里去掉origQual >= preserveQLessThan的等号，是为了过滤掉碱基为Ｎ的情况
        qualKey = QualitySignal::getQualIndex(origQual, qual_line);//查询该qual在qualitytable中第几行
//        cout << "qualKey = " << qualKey << endl;
        empQS = tables[1][qualKey];
//        cout << "empQS = " << (int)empQS << endl;
        deltaQCovariates = 0;
        IdeltaQCovariates = 0;
        DdeltaQCovariates = 0;
        contextKey = ContextSignal::getMisMatchContextSignal(bases, readLength, offset, NegativeStrandFlag);
        if(contextKey >= 0){
            IcontextKey = ContextSignal::getIndelContextSignal(bases, readLength, offset, NegativeStrandFlag, contextKey);
            if(IcontextKey >= 0){
                IcontextKey <<= 1;//因为I D 情况　所以乘以２
                IcontextKey += contextMSize;//contextMsize是M情况下Context的行数
//                    cout << "IcontextKey = " << IcontextKey << endl;
                IdeltaQCovariates += tables[2][IcontextKey];
                DdeltaQCovariates += tables[2][IcontextKey + 1];
//                    cout << "IdeltaQCovariates1 = " << (int)IdeltaQCovariates << endl;
//                    cout << "DdeltaQCovariates1 = " << (int)DdeltaQCovariates << endl;
            }
            contextKey += qualKey << 4;
//                cout << "contextKey = " << contextKey << endl;
            deltaQCovariates += tables[2][contextKey];
//                cout << "deltaQCovariates1 = " << (int)deltaQCovariates << endl;
        }

        cycleKey = CycleSignal::getCycleSignal(offset, NegativeStrandFlag, readOrderFactor, readLength);
        if(cycleKey >= 0) {
            IcycleKey = CycleSignal::getIndelCycleSignal(offset, readLength, cycleKey);
            if(IcycleKey >= 0){
                IcycleKey += cycleMSize;//cycleMSize是M情况下Cycle的行数
                IdeltaQCovariates += tables[3][IcycleKey];
                DdeltaQCovariates += tables[3][IcycleKey + 1];
//                    cout << "IdeltaQCovariates2 = " << (int)IdeltaQCovariates << endl;
//                    cout << "DdeltaQCovariates2 = " << (int)DdeltaQCovariates << endl;
            }
            cycleKey += qualKey * max2Cycle;
//                cout << "cycleKey = " << cycleKey << endl;
            deltaQCovariates += tables[3][cycleKey];
//                cout << "deltaQCovariates2 = " << (int)deltaQCovariates << endl;
        }

        quals[offset] = MathUtils::max(MathUtils::min(empQS + deltaQCovariates, MAX_RECALIBRATED_Q_SCORE), 1);
        Iquals[offset] = MathUtils::max(MathUtils::min(IempQS + IdeltaQCovariates, MAX_RECALIBRATED_Q_SCORE), 1);
        Dquals[offset] = MathUtils::max(MathUtils::min(DempQS + DdeltaQCovariates, MAX_RECALIBRATED_Q_SCORE), 1);
//            cout << "empQS + deltaQCovariates = " << (int)(empQS + deltaQCovariates) << endl;
//            cout << "MAX_RECALIBRATED_Q_SCORE = " << (int) MAX_RECALIBRATED_Q_SCORE << endl;
//            quals[offset] = MathUtils::min((empQS + deltaQCovariates), MAX_RECALIBRATED_Q_SCORE);
//            Iquals[offset] = IempQS + IdeltaQCovariates;
//            Dquals[offset] = DempQS + DdeltaQCovariates;
//            if(origQual == preserverQLessThan)
//                quals[offset] = preserverQLessThan;
        if(origQual < preserveQLessThan)
            quals[offset] = origQual;
//        }
    }

//    for(int i = 0; i < readLength; i++){
//        cout << "quals[" << i << "] = " << (int)quals[i] << ", ";
//        cout << "Iquals[" << i << "] = " << (int)Iquals[i] << ", ";
//        cout << "Dquals[" << i << "] = " << (int)Dquals[i] << ", ";
//    }
//    for(int i = 0; i < readLength; i++){
//        cout << (int)quals[i] << ", ";
//    }
//    cout << endl;
//    for(int i = 0; i < readLength; i++){
//        cout << (int)Iquals[i] << ", ";
//    }
//    cout << endl;
//    for(int i = 0; i < readLength; i++){
//        cout << (int)Dquals[i] << ", ";
//    }
//    cout << endl;
//    cout << "/////////////////////////////"<< endl;
    return;
}

BaseRecalibration::BaseRecalibration() {

}


vector<int> BaseRecalibration::getQualLine(){
    return qual_line;
}