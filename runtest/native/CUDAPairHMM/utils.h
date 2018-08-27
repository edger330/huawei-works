//
// Created by cjw on 17-5-24.
//

#ifndef CUDAPAIRHMM_UTILS_H
#define CUDAPAIRHMM_UTILS_H


#include "common_data_structure.h"
#include <iostream>
#include <vector>

struct dataSet{
    int taskID,sampleID;
    int numReads,numHaps;
    int * readsLens;
    int * hapsLens;
    const char * readsData;
    const char * hapsData;
};
struct resultSet{
    int taskID,sampleID;
    int numReads,numHaps;
    float * result;
};
resultSet * cuda_compute_full_prob(dataSet * data);

#endif //CUDAPAIRHMM_UTILS_H
