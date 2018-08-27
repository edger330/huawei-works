//
// Created by shawn on 18-3-20.
//

#ifndef BASERECALIBRATIONCMETHOD_PARSEGROUPFILE_H
#define BASERECALIBRATIONCMETHOD_PARSEGROUPFILE_H

#include <vector>

using namespace std;

class ParseGroupFile {
public:
    void loadTables(vector<vector<char> > & _tables, const char * tables_str, int & _contextMSize,
                    int & _cycleMSize, int & _max2Cycle, bool & _tableIsNull, int & _qualTableSize, vector<int> & qual_line);

};


#endif //BASERECALIBRATIONCMETHOD_PARSEGROUPFILE_H
