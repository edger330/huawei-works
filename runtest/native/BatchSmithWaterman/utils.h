#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>
#include <string>
using namespace std;

typedef struct {
	int w_open;
	int w_extend;
	int w_match;
	int w_mismatch;
}Parameter;

enum State{
        MATCH,
        INSERTION,
        DELETION,
        CLIP
};

typedef struct testcase{
    int refLen;
    int altLen;
    char * reference;
    char * alternate;
};

inline const char* ToString(State v)
{
    switch (v)
    {
        case MATCH:   return "MATCH";
        case INSERTION:   return "INSERTION";
        case DELETION: return "DELETION";
        case CLIP:      return "CLIP";
        default: return "NONE";
    }
}

inline string int_to_string(int a){
    stringstream ss;
    ss<<a;
    string str = ss.str();
    return str;
}

#endif