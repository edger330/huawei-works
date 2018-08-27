//
// Created by shawn on 18-3-20.
//
#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "org_broadinstitute_gatk_engine_recalibration_BaseRecalibration.h"
#include "BaseRecalibration.h"
#include "ParseGroupFile.h"
#include <sys/time.h>

#define errorModel 0
using namespace std;

BaseRecalibration * case_bqsr;
BaseRecalibration * normal_bqsr;

//JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_engine_recalibration_BaseRecalibration_initializeCfunction
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_tools_walkers_cancer_m2_FastMuTect2_initializeCfunction
(JNIEnv * env, jobject obj, jstring case_tables, jstring normal_tables)
{
    vector<vector<char> > recalTables;
    //load case table
    //todo 实现传递table path的接口
    const char * tables_str = env->GetStringUTFChars(case_tables, 0);
//    const char * tables_str = "";

    // printf("table == %s", tables_str);
    ParseGroupFile p1;

    int _contextMSize, _cycleMSize, _max2Cycle, _qualTableSize;
    bool _tableIsNull;
    vector<int> _qual_line;
    p1.loadTables(recalTables, tables_str, _contextMSize,
        _cycleMSize,  _max2Cycle, _tableIsNull, _qualTableSize, _qual_line); //Load case recalibrateTables here.
    env->ReleaseStringUTFChars(case_tables, tables_str);

    case_bqsr = new BaseRecalibration(recalTables, _contextMSize, _cycleMSize, _max2Cycle, _tableIsNull, _qualTableSize, _qual_line);//todo BaseRecalibration constructor
    //load case table end

    //    recalTables.clear();


    //load normal table
    vector<vector<char> > recalTables2;
//    const char * tables_str2 = "";
    const char * tables_str2 = env->GetStringUTFChars(normal_tables, 0);
    _qual_line.clear();
    p1.loadTables(recalTables2, tables_str2, _contextMSize,
        _cycleMSize,  _max2Cycle, _tableIsNull, _qualTableSize, _qual_line); //Load normal recalibrateTables here.
    env->ReleaseStringUTFChars(normal_tables, tables_str2);


    normal_bqsr = new BaseRecalibration(recalTables2, _contextMSize, _cycleMSize, _max2Cycle, _tableIsNull, _qualTableSize, _qual_line);
    //load normal table end
}


JNIEXPORT jobjectArray JNICALL Java_org_broadinstitute_gatk_tools_walkers_cancer_m2_FastMuTect2_recalibrateCmethod
        (JNIEnv * env, jobject obj, jint readLength, jboolean NegativeStrandFlag, jint readOrderFactor, jbyteArray bases, jbyteArray quals, jint rgSignal)
{
    jbyte c_bases[readLength];
    jbyte c_quals[readLength];
    jbyte Iquals[readLength];
    jbyte Dquals[readLength];
    env->GetByteArrayRegion(bases, 0, readLength, c_bases);
    env->GetByteArrayRegion(quals, 0, readLength, c_quals);

    // time_t start, end;
    // start = time
    // struct timespec prev_time;
    // clock_gettime(CLOCK_REALTIME, &prev_time);

    //todo 这里测试需要判断rgSignal是否等于0, 如果测试无误，可以将else if中的判断去掉
    if(rgSignal == 1){
//        vector<int> test = case_bqsr->getQualLine();
//        cout << "#####case" << endl;
//        for(int i = 0; i < test.size(); i++){
//            cout << test[i] << endl;
//        }
        case_bqsr->recalibrateCmethod(readLength, NegativeStrandFlag, readOrderFactor, (char *)c_bases, (char *)c_quals, rgSignal, (char *)Iquals, (char *)Dquals);
    } else if(rgSignal == 0){
//        vector<int> test = normal_bqsr->getQualLine();
//        cout << "#####normal" << endl;
//        for(int i = 0; i < test.size(); i++){
//            cout << test[i] << endl;
//        }
        normal_bqsr->recalibrateCmethod(readLength, NegativeStrandFlag, readOrderFactor, (char *)c_bases, (char *)c_quals, rgSignal, (char *)Iquals, (char *)Dquals);
    } else{
        cout << "There is an error in recalibrateCmethod." << endl;
    }
    // struct timespec curr_time;
    // clock_gettime(CLOCK_REALTIME, &curr_time);
    // cerr << (uint64_t)((curr_time.tv_sec-prev_time.tv_sec)*1000000000+(curr_time.tv_nsec-prev_time.tv_nsec)) << endl;
    // for(int i = 0; i < readLength; i++)
    //     cerr << (int)c_quals[i] << ", ";
    jbyteArray j_quals = env->NewByteArray(readLength);
    jbyteArray jI_quals = env->NewByteArray(readLength);
    jbyteArray jD_quals = env->NewByteArray(readLength);
    jobjectArray result;
    jclass intArrCls = env->FindClass("[B");
    result = env->NewObjectArray(3, intArrCls, NULL);

    env->SetByteArrayRegion(j_quals, 0, readLength, c_quals);
    env->SetObjectArrayElement(result, 0, j_quals);

    env->SetByteArrayRegion(jI_quals, 0, readLength, Iquals);
    env->SetObjectArrayElement(result, 1, jI_quals);

    env->SetByteArrayRegion(jD_quals, 0, readLength, Dquals);
    env->SetObjectArrayElement(result, 2, jD_quals);

    env->DeleteLocalRef(j_quals);
    env->DeleteLocalRef(jI_quals);
    env->DeleteLocalRef(jD_quals);
    return result;
}