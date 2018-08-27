//
// Created by shawn on 18-3-20.
//

#include "SetTables.h"
#include "QualitySignal.h"

void SetTables::setRGRecalTable(vector<string> node, vector<char> &recalTable) {
    recalTable.push_back((char)atoi(node[2].c_str()));
}

void SetTables::setQSRecalTable(vector<string> node, vector<char> &recalTable) {
    recalTable.push_back((char)atoi(node[3].c_str()));
}

void SetTables::setCovRecalTable(vector<string> node, vector<char> &recalTable, vector<char> qualRecalTable, vector<int> line) {
    if(node[4] == "M")
        recalTable.push_back((char)(atoi(node[5].c_str()) - qualRecalTable[QualitySignal::getQualIndex(atoi(node[1].c_str()), line)]));
    else if(node[4] == "I")
        recalTable.push_back((char)(atoi(node[5].c_str()) - qualRecalTable[(qualRecalTable.size()-2)]));
    else if(node[4] == "D")
        recalTable.push_back((char)(atoi(node[5].c_str()) - qualRecalTable[(qualRecalTable.size()-1)]));
    else
        cout << "something wrong happens in setTables::setCovRecalTable function." << endl;
}

void SetTables::setCycleRecalTable(vector<string> node, vector<char> &recalTable, vector<char> qualRecalTable, int maxCycle, vector<int> line){
    int key;
    int sign = 0;
    int cycle = atoi( node[2].c_str() );
//    cout << cycle << endl;
    if (cycle < 0){
        sign = 1;
    }

    if(node[4] == "M"){
        key = QualitySignal::getQualIndex( atoi( node[1].c_str() ), line ) * 2 * maxCycle;
        key += ( abs(cycle) - 1 ) * 2 + sign;
        recalTable[key] = (char)( atoi( node[5].c_str() ) - qualRecalTable[ QualitySignal::getQualIndex( atoi(node[1].c_str()), line ) ]);
    } else if(node[4] == "I"){
        key = (maxCycle * 2 * ( (int) qualRecalTable.size() - 2));
        key += ( abs(cycle) - 5 ) * 4 + sign * 2;
        recalTable[key] = (char)(atoi(node[5].c_str()) - qualRecalTable[(qualRecalTable.size()-2)]);
    } else if(node[4] == "D"){
        key = (maxCycle * 2 * ( (int) qualRecalTable.size() - 2));
        key += ( abs(cycle) - 5 ) * 4 + sign * 2 + 1;
        recalTable[key] = (char)(atoi(node[5].c_str()) - qualRecalTable[(qualRecalTable.size()-1)]);
    } else {
        cout << "something wrong happens in setTables::setCycleRecalTable function." << endl;
    }

}