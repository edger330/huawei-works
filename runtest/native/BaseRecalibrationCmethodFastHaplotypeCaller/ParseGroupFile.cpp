//
// Created by shawn on 18-3-20.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include "ParseGroupFile.h"
#include "SetTables.h"

void ParseGroupFile::loadTables(vector<vector<char> > & _tables, const char * tables_str, int & _contextMSize,
                                int & _cycleMSize, int & _max2Cycle, bool & _tableIsNull, int & _qualTableSize, vector<int> & qual_line){
    //recalibrateCmethod没有用到 maxCycle与rgTableSize所以
    int maxCycle = 0;
    int rgTableSize = 0;
    _contextMSize = 0;
    _cycleMSize = 0;
    _max2Cycle = 0;
    _tableIsNull = false;
    _qualTableSize = 0;

    string line;
    ifstream in(tables_str);
    if(!in.is_open()) {
        cout << "can't load " << tables_str << "!";
        exit(1);
    }
    int counter = 0;
    string word;
    vector<string> node;
    vector<char> recalTable;

    int signal = 0;

    int cycleTable_index = 0;//标记cycleTable从表中哪一行开始

    while(getline(in, line)) {
        cycleTable_index++;
        if(counter == 4 && signal == (_qualTableSize - 2)*16 + 2*64) {//将context table加入到tables中 todo signal随quality值变化
            _tables.push_back(recalTable);//C++ vector中的push_back方法是值拷贝值而不是地址拷贝
            recalTable.clear();//将上一个表的信息传入_tables结构后清空表，为读取下一个表做准备

            counter++;//开始读取cycle表
        }
        if(line.length() <= 1) {//遇到空行，说明接下来要开始读一张新表了
            counter++;
            if(counter > 2 && counter < 6) {
                _tables.push_back(recalTable);
                recalTable.clear();
            }
            if(counter < 5 && counter >= 2) {//忽略每张表的header信息
                for(int i = 0; i <= 3; i++) {
                    getline(in, line);
                    cycleTable_index++;
                }
            }
        }
        if(counter < 2)
            continue;//忽略GATK argument table 和　quantization table

        istringstream record(line);
        while(record >> word) {//把record中的一行字段逐个赋值给word，传入node中
            node.push_back(word);
        }
        switch (counter) {
            case 2:
                if(node.size() == 0){
//                    cout << "_tableIsNull == " << _tableIsNull << endl;
                    _tableIsNull = true;
//                    cout << "now _tableIsNull == " << _tableIsNull << endl;
                    return;
                }
//                cout << "node.size()== " << node.size() << endl;
                SetTables::setRGRecalTable(node, recalTable);
                rgTableSize++;
                break;
            case 3:
                SetTables::setQSRecalTable(node, recalTable);
                qual_line.push_back(atoi(node[1].c_str()));
                _qualTableSize++;
                break;
            case 4://目前已经采取一定方法，可以对任意多的quality种类进行识别读取，但Insertion quality和deletion　quality只能各有一行　todo 将context表和cycle表的划分方法规范的重写
                SetTables::setCovRecalTable(node, recalTable, _tables[1], qual_line);
                signal++;//标记Context表和cycle表的边界
                break;
            case 5:
                cycleTable_index--;
                if( abs( atoi(node[2].c_str()) ) > maxCycle ){
                    maxCycle = abs( atoi(node[2].c_str()) );
                }


//                SetTables::setCycleRecalTable(node, recalTable, _tables[1]);
                break;
            default:
//                cout << "table doesn't exit!" << endl;
                cout << "load table end!" << endl;
        }
        node.clear();
    }

    ifstream in2(tables_str);
    if(!in2.is_open()) {
        cout << "can't load table " << tables_str << " again !";
        exit(1);
    }
    int index = 1;

    int size = maxCycle * 2 * (_qualTableSize - 2) + (maxCycle - 8 ) * 2 * 2;//预先计算好recalCycleTable的大小
//    cout << "size = " << size << endl;
    vector<char> recalCycleTable(size, 0);

    while(getline(in2, line)) {
        if (index < cycleTable_index){
            index++;
        } else if(line.length() <= 1){//当遇到了表末尾的空行
            break;
        } else if(index >= cycleTable_index){//当读取到Cycle表的第一行
            istringstream record(line);
            while(record >> word) {//把record中的一行字段逐个赋值给word，传入node中
                node.push_back(word);
            }
            SetTables::setCycleRecalTable(node, recalCycleTable, _tables[1], maxCycle, qual_line);
            node.clear();
        }
    }

    _tables.push_back(recalCycleTable);
    recalCycleTable.clear();


    rgTableSize /= 3;//每个readGroup都有MDI三种情况
    _contextMSize = 16 * (_qualTableSize - 2);//M的节点数,之所以减二是qualTable中最后的两行且只有最后两行是I和D的情况 todo 增强扩展性用mismatchContextSize来代替常数16


    //todo maxCycle计算方式
//    maxCycle = ((int)_tables[3].size() + 32)/(2 * _qualTableSize);//计算得到的公式，　MaxCycle就是Cycle的最大可能值 这里认为　I D 情况总是比 M 情况首尾少4个
    _cycleMSize = maxCycle * 2 * (_qualTableSize - 2);//M的节点数,之所以减二是qualTable中最后的两行且只有最后两行是I和D的情况 todo 增强扩展性用最大cycle值来代替常数76
    _max2Cycle = 2 * maxCycle;
    //test
//    for(int i = 0; i < 4; i++){
//        for(int j = 0; j < _tables[i].size(); j++)
//            cout << i << "--" << j << " == " << (int)_tables[i][j] << endl;
//    }
//    cout << "-----------end------------- " << endl;

}

