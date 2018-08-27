#include "org_broadinstitute_gatk_utils_baq_BAQ.h"
#include "BAQ.h"

#include<iostream>
BAQ * baq;
using namespace std;

JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_baq_BAQ_CinitializeCachedData
  (JNIEnv * env, jobject obj)
{
	
	baq = new BAQ();
}


JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_baq_BAQ_Chmm_1glocal
  (JNIEnv * env, jobject obj, jbyteArray ref, jint l_ref, jbyteArray query, jint query_length, jint qstart, jint l_query, jbyteArray iqual, jintArray state, jbyteArray q)
{
	jbyte c_ref[l_ref];
	jbyte c_query[query_length];
	jbyte c_iqual[query_length];
	jint c_state[query_length];
	jbyte c_q[query_length];
	env->GetByteArrayRegion(ref, 0, l_ref, c_ref);
	env->GetByteArrayRegion(query, 0, query_length, c_query);
	env->GetByteArrayRegion(iqual, 0, query_length, c_iqual);

	// for(int i = 0; i < l_ref; i++)
	// 	cout << (int)c_ref[i] << ", ";

	// cout << endl;
	// for(int i = 0; i < query_length; i++)
	// 	cout << (int)c_query[i] << ", ";

	// cout << endl;
	// for(int i = 0; i < query_length; i++)
	// 	cout << (int)c_iqual[i] << ", ";
	// cout << endl;
	
	baq->Chmm_glocal((char *)c_ref, l_ref, (char *)c_query, query_length, qstart, l_query, (char *)c_iqual, (int *)c_state, (char *)c_q);

//	 cout << "Chmm_glocal----"<< endl;
	// for(int i = 0; i < query_length; i++)
	// 	cout << c_state[i] << ", ";

	// cout << endl;
	// for(int i = 0; i < query_length; i++)
	// 	cout << (int)c_q[i] << ", ";

	// cout << endl;
//	 cout << "Chmm_glocal end ----" << endl;
	env->SetIntArrayRegion(state, 0, query_length, c_state);
	env->SetByteArrayRegion(q, 0, query_length, c_q);
//    env->DeleteLocalRef((int*)c_state);
//    env->DeleteLocalRef((char *)c_q);
}

