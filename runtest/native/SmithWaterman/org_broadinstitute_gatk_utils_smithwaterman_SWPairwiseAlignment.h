/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class org_broadinstitute_gatk_utils_smithwaterman_SWPairwiseAlignment */

#ifndef _Included_org_broadinstitute_gatk_utils_smithwaterman_SWPairwiseAlignment
#define _Included_org_broadinstitute_gatk_utils_smithwaterman_SWPairwiseAlignment
#ifdef __cplusplus
extern "C" {
#endif

/*
 * Class:     org_broadinstitute_gatk_utils_smithwaterman_SWPairwiseAlignment
 * Method:    calculateMatrix
 * Signature: ([B[BLjava/util/ArrayList;Ljava/util/ArrayList;)I
 */
JNIEXPORT jint JNICALL Java_org_broadinstitute_gatk_utils_smithwaterman_SWPairwiseAlignment_calculateMatrix
  (JNIEnv * env, jobject obj, jbyteArray reference, jbyteArray alternate,
        jint w_match, jint w_mismatch, jint w_open, jint w_extend, jstring strategyName,
        jobject cigars_op, jobject cigars_len);

#ifdef __cplusplus
}
#endif
#endif
