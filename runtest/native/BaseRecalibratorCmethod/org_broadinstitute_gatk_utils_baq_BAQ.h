/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class org_broadinstitute_gatk_utils_baq_BAQ */

#ifndef _Included_org_broadinstitute_gatk_utils_baq_BAQ
#define _Included_org_broadinstitute_gatk_utils_baq_BAQ
#ifdef __cplusplus
extern "C" {
#endif
#undef org_broadinstitute_gatk_utils_baq_BAQ_DEBUG
#define org_broadinstitute_gatk_utils_baq_BAQ_DEBUG 0L
#undef org_broadinstitute_gatk_utils_baq_BAQ_DEFAULT_GOP
#define org_broadinstitute_gatk_utils_baq_BAQ_DEFAULT_GOP 40.0
#undef org_broadinstitute_gatk_utils_baq_BAQ_EM
#define org_broadinstitute_gatk_utils_baq_BAQ_EM 0.33333333333
#undef org_broadinstitute_gatk_utils_baq_BAQ_EI
#define org_broadinstitute_gatk_utils_baq_BAQ_EI 0.25
/*
 * Class:     org_broadinstitute_gatk_utils_baq_BAQ
 * Method:    CinitializeCachedData
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_baq_BAQ_CinitializeCachedData
  (JNIEnv *, jobject);

/*
 * Class:     org_broadinstitute_gatk_utils_baq_BAQ
 * Method:    Chmm_glocal
 * Signature: ([BI[BIII[B[I[B)V
 */
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_baq_BAQ_Chmm_1glocal
  (JNIEnv *, jobject, jbyteArray, jint, jbyteArray, jint, jint, jint, jbyteArray, jintArray, jbyteArray);

#ifdef __cplusplus
}
#endif
#endif
