//
// Created by cjw on 17-5-24.
//
#include "utils.h"

//Static members in ContextBase
template<>
bool ContextBase<double>::staticMembersInitializedFlag = false;
template<>
double ContextBase<double>::jacobianLogTable[JACOBIAN_LOG_TABLE_SIZE] = { };
template<>
double ContextBase<double>::matchToMatchProb[((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] = { };
template<>
bool ContextBase<float>::staticMembersInitializedFlag = false;
template<>
float ContextBase<float>::jacobianLogTable[JACOBIAN_LOG_TABLE_SIZE] = { };
template<>
float ContextBase<float>::matchToMatchProb[((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] = { };
