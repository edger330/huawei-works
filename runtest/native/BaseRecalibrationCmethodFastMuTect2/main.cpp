#include <iostream>
#include <vector>
#include "BaseRecalibration.h"
#include "ParseGroupFile.h"

using namespace std;


int main() {
    vector<vector<char>> recalTables;

    ParseGroupFile p1;
    int _contextMSize;
    int _cycleMSize;
    int _max2Cycle;
    bool _tableIsNull;
    int _qualTableSize;
    vector<int> qual_line(0);
    p1.loadTables(recalTables, "/home/shawn/bqsr_io/tiny_output/1_case.grp", _contextMSize, _cycleMSize, _max2Cycle, _tableIsNull, _qualTableSize, qual_line);
    cout << "_contextMSize = " << _contextMSize << endl;
    cout << "_cycleMSize = " << _cycleMSize << endl;
    cout << "_max2Cycle = " <<_max2Cycle << endl;
    cout << "_tableIsNull = " << _tableIsNull << endl;
    cout << "_qualTableSize = " << _qualTableSize << endl;
    cout << "qual_line = " << endl;
    for(int i = 0; i < qual_line.size(); i++){
        cout<< qual_line[i] << endl;
    }
//    cout << _contextMSize << _cycleMSize << _max2Cycle << _tableIsNull << _qualTableSize << endl;

    BaseRecalibration bqsr(recalTables, _contextMSize, _cycleMSize, _max2Cycle, _tableIsNull, _qualTableSize, qual_line);

    BaseRecalibration a;
    int readNum = 1;
    int readLength = 50;
    bool NegativeStrandFlag = false;
    int readOrderFactor = 0;
    int rgSignal = 0;
    char bases[] = {71, 71, 67, 65, 71, 65, 71, 71, 84, 84, 71, 67, 78, 71, 84, 71, 65, 71, 67, 67, 65, 65, 71, 65, 84, 67, 65, 84, 71, 67, 67, 65, 67, 84, 71, 67, 65, 67, 84, 67, 67, 65, 71, 67, 84, 71, 71, 71, 71, 67};
    char quals[] = {37, 37, 37, 37, 37, 36, 31, 37, 37, 37, 36, 25, 6, 37, 37, 38, 25, 38, 37, 37, 36, 36, 37, 29, 30, 37, 36, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37};
//    char bases[] = {71, 71, 67, 67, 65, 71, 65, 67, 65, 71, 65, 84, 84, 67, 65, 67, 65, 71, 67, 65, 71, 65, 65, 84, 84, 67, 84, 65, 67, 67, 65, 71, 65, 67, 65, 84, 84, 67, 65, 65, 65, 71, 65, 65, 84, 71, 84, 67, 84, 84, 67, 84, 84, 84, 67, 65, 84, 84, 67, 65, 65, 71, 71, 65, 65, 71, 65, 65, 65, 84, 71, 65, 84, 65, 67, 67};
//    char quals[] = {32, 32, 32, 32, 32, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 32, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 32, 36, 36, 36, 36, 36, 36, 32, 36, 32, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36};
    char Iquals[50] = {0};
    char Dquals[50] = {0};
    for(int i = 0; i < readNum; i++){
        bqsr.recalibrateCmethod(readLength, NegativeStrandFlag, readOrderFactor, bases, quals, rgSignal, Iquals, Dquals);
//        string quals_string(quals);
    }


    return 0;
}
