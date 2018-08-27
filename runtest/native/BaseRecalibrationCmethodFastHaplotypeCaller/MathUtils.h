//
// Created by shawn on 18-3-20.
//

#ifndef BASERECALIBRATIONCMETHOD_MATHUTILS_H
#define BASERECALIBRATIONCMETHOD_MATHUTILS_H


class MathUtils {
public:
    static char max(char a, char b){
        if(a > b)
            return a;
        else
            return b;
    };

    static char min(char a, char b){
        if(a < b)
            return a;
        else
            return b;
    };

};


#endif //BASERECALIBRATIONCMETHOD_MATHUTILS_H
