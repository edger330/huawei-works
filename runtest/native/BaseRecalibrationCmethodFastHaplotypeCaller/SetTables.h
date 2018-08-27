//
// Created by shawn on 18-3-20.
//

#ifndef BASERECALIBRATIONCMETHOD_SETTABLES_H
#define BASERECALIBRATIONCMETHOD_SETTABLES_H

#include <vector>
#include <iostream>
#include <cstdlib>

using namespace std;
class SetTables {
public:

    static void setRGRecalTable(vector<string> node, vector<char> &recalTable);
    static void setQSRecalTable(vector<string> node, vector<char> &recalTable);
    static void setCovRecalTable(vector<string> node, vector<char> &recalTable, vector<char> qualRecalTable, vector<int> line);
    static void setCycleRecalTable(vector<string> node, vector<char> &recalTable, vector<char> qualRecalTable, int maxCycle, vector<int> line);

};


#endif //BASERECALIBRATIONCMETHOD_SETTABLES_H
