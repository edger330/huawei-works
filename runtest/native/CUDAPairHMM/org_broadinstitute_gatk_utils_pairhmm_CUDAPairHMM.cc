//
// Created by cjw on 17-5-24.
//

#include <vector>
#include <string.h>
#include "org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM.h"
#include "utils.h"
#include "LoadTimeInitializer.h"
//#include "SyncQueue.h"

using namespace std;

//阻塞队列传输数据到gpu pairHMM
//SyncQueue<dataSet> dataQueue(10);
//阻塞队列存储gpu pairHMM返回的数据
//SyncQueue<resultSet> resultQueue(10);

/*
void cudaThread(){
    try {
        while (true) {
            dataSet data;
            dataQueue.Take(data);
            //todo;

        }
    }catch (exception e){
        throw e;
    }
}
static std::thread t1(cudaThread);
*/



JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM_initCUDA
(JNIEnv * env, jobject obj, jclass readDataHolderClass, jclass haplotypeDataHolderClass){
    assert(readDataHolderClass);
    assert(haplotypeDataHolderClass);
    jfieldID fid;
    fid = env->GetFieldID(readDataHolderClass, "readBases", "[B");
    assert(fid && "JNI pairHMM: Could not get FID for readBases");
    g_load_time_initializer.m_readBasesFID = fid;
    fid = env->GetFieldID(readDataHolderClass, "readQuals", "[B");
    assert(fid && "JNI pairHMM: Could not get FID for readQuals");
    g_load_time_initializer.m_readQualsFID = fid;
    fid = env->GetFieldID(readDataHolderClass, "insertionGOP", "[B");
    assert(fid && "JNI pairHMM: Could not get FID for insertionGOP");
    g_load_time_initializer.m_insertionGOPFID = fid;
    fid = env->GetFieldID(readDataHolderClass, "deletionGOP", "[B");
    assert(fid && "JNI pairHMM: Could not get FID for deletionGOP");
    g_load_time_initializer.m_deletionGOPFID = fid;
    fid = env->GetFieldID(readDataHolderClass, "overallGCP", "[B");
    assert(fid && "JNI pairHMM: Could not get FID for overallGCP");
    g_load_time_initializer.m_overallGCPFID = fid;

    fid = env->GetFieldID(haplotypeDataHolderClass, "haplotypeBases", "[B");
    assert(fid && "JNI pairHMM: Could not get FID for haplotypeBases");
    g_load_time_initializer.m_haplotypeBasesFID = fid;
}
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM_doneCUDA
(JNIEnv * env, jobject obj){

//todo;
}

JNIEXPORT jfloatArray JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM_sendDataAndResiveFromCUDA
(JNIEnv * env, jobject obj, jint taskID, jint sampleID, jint numReads, jint numHaplotypes,
    jobjectArray readDataArray, jobjectArray haplotypeDataArray){

    int readLens[numReads];
    int hapLens[numHaplotypes];
    int readLenSum=0,hapLenSum=0;
//    string rs,qs,is,ds,cs,haps;
    dataSet batchData;
//    char * readData;
//    char * hapData;
    jboolean is_copy= JNI_FALSE;
    for(int i=0;i<numReads;++i){
        jobject readObject = env->GetObjectArrayElement(readDataArray, i);
        jbyteArray readBases = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_readBasesFID);
        jsize readLength = env->GetArrayLength(readBases);
        readLens[i]=readLength;
        readLenSum+=readLength;
    }
    for(int i=0;i<numHaplotypes;++i){
        jobject hapObject = env->GetObjectArrayElement(haplotypeDataArray, i);
        jbyteArray hapBases = (jbyteArray)env->GetObjectField(hapObject, g_load_time_initializer.m_haplotypeBasesFID);
        jsize hapLength = env->GetArrayLength(hapBases);
        hapLens[i]=hapLength;
        hapLenSum+=hapLength;
    }

    char * reads = new char[readLenSum*5];
    char * haps = new char[hapLenSum];
    int offset = 0;
    for(int i=0;i<numReads;++i){
        jobject readObject = env->GetObjectArrayElement(readDataArray, i);
        jbyteArray readBases = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_readBasesFID);
        jbyteArray insertionGOP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_insertionGOPFID);
        jbyteArray deletionGOP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_deletionGOPFID);
        jbyteArray overallGCP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_overallGCPFID);
        jbyteArray readQuals = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_readQualsFID);

        jsize readLength = env->GetArrayLength(readBases);
        jbyte* readBasesArray = (jbyte*)env->GetByteArrayElements(readBases, &is_copy);	//order of GET-RELEASE is important
        jbyte* readQualsArray = (jbyte*)env->GetByteArrayElements(readQuals, &is_copy);
        jbyte* insertionGOPArray = (jbyte*)env->GetByteArrayElements(insertionGOP, &is_copy);
        jbyte* deletionGOPArray = (jbyte*)env->GetByteArrayElements(deletionGOP, &is_copy);
        jbyte* overallGCPArray = (jbyte*)env->GetByteArrayElements(overallGCP, &is_copy);

        memcpy(reads+offset, readBasesArray, readLength);
        memcpy(reads+readLenSum+offset, readQualsArray, readLength);
        memcpy(reads+readLenSum*2+offset, insertionGOPArray, readLength);
        memcpy(reads+readLenSum*3+offset, deletionGOPArray, readLength);
        memcpy(reads+readLenSum*4+offset, overallGCPArray, readLength);
        offset+=readLength;
    }

    offset = 0;
    for(int i=0;i<numHaplotypes;++i){
        jobject hapObject = env->GetObjectArrayElement(haplotypeDataArray, i);
        jbyteArray hapBases = (jbyteArray)env->GetObjectField(hapObject, g_load_time_initializer.m_haplotypeBasesFID);
        jsize hapLength = env->GetArrayLength(hapBases);

        jbyte* hapBasesArray = (jbyte*)env->GetByteArrayElements(hapBases, &is_copy);	//order of GET-RELEASE is important

        hapLens[i]=hapLength;
        memcpy(haps+offset,hapBasesArray,hapLength);
        offset += hapLength;
    }
//    hapData= const_cast<char*>(haps.c_str());
    batchData.taskID=taskID;
    batchData.sampleID=sampleID;
    batchData.numReads=numReads;
    batchData.numHaps=numHaplotypes;
    batchData.readsLens=readLens;
    batchData.hapsLens=hapLens;
    batchData.readsData=reads;
    batchData.hapsData=haps;
   //    dataQueue.Put(batchData);
   resultSet * batchResult = cuda_compute_full_prob(&batchData);
    delete[] reads;
    delete[] haps;

    int RtaskID = batchResult->taskID;
    int RsampleID = batchResult->sampleID;
    int RnumReads=batchResult->numReads;
    int RnumHaps=batchResult->numHaps;
    assert(RtaskID==taskID);
    assert(RsampleID==sampleID);
    assert(RnumReads=numReads);
    assert(RnumHaps==numHaplotypes);

    float * result=batchResult->result;

    jfloatArray likehoods;
    likehoods=env->NewFloatArray(RnumReads*RnumHaps);
    env->SetFloatArrayRegion(likehoods,0,RnumReads*RnumHaps,result);
    free(batchResult ->result);
    delete batchResult;
    return likehoods;

}

/*
 * Class:     org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM
 * Method:    resiveFromCUDA
 * Signature: (II)[D
 */
/*JNIEXPORT jfloatArray JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_CUDAPairHMM_resiveFromCUDA
        (JNIEnv * env, jobject obj , jint & taskID, jint & sampleID){

        dataSet data;
        dataQueue.Take(data);
        resultSet * batchResult = cuda_compute_full_prob(&data);

        taskID = batchResult->taskID;
        sampleID = batchResult->sampleID;
        jint numReads=batchResult->numReads;
        jint numHaps=batchResult->numHaps;
        float * result=batchResult->result;

        jfloatArray likehoods;
        likehoods=env->NewFloatArray(numReads*numHaps);
        env->SetFloatArrayRegion(likehoods,0,numReads*numHaps,result);
		free(batchResult ->result);
		free(batchResult);
        return likehoods;
}*/
