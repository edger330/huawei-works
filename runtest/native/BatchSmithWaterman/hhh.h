/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class org_broadinstitute_gatk_utils_smithwaterman_BatchSWPairwiseAlignment */

#ifndef _Included_org_broadinstitute_gatk_utils_smithwaterman_BatchSWPairwiseAlignment
#define _Included_org_broadinstitute_gatk_utils_smithwaterman_BatchSWPairwiseAlignment
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     org_broadinstitute_gatk_utils_smithwaterman_BatchSWPairwiseAlignment
 * Method:    jniInit
 * Signature: (Lorg/broadinstitute/gatk/utils/smithwaterman/BatchSWPairwiseAlignment/OVERHANG_STRATEGY;Lorg/broadinstitute/gatk/utils/smithwaterman/Parameters;Z)V
 */
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_smithwaterman_BatchSWPairwiseAlignment_jniInit
  (JNIEnv * env, jobject obj,jobject strategy, jobject parameter, jboolean cutoff,jobject cigar, jobject cigarelement,jobject op);

/*
 * Class:     org_broadinstitute_gatk_utils_smithwaterman_BatchSWPairwiseAlignment
 * Method:    jniAlign
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;)V
 */
JNIEXPORT jstring JNICALL Java_org_broadinstitute_gatk_utils_smithwaterman_BatchSWPairwiseAlignment_jniAlign
  (JNIEnv *env, jobject obj, jobject haplotypes, jobject reads, jobject cigars, jobject alignment_offsets);

#ifdef __cplusplus
}
#endif
#endif
