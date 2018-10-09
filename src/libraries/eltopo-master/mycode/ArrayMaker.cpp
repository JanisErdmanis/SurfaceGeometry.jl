#include "ArrayMaker.h"
#include <iostream>

using namespace std;

ArrayMaker::ArrayMaker(int iNum, float fNum) {
    cout << "Got arguments: " << iNum << ", and " << fNum << endl;
    iNumber = iNum;
    fNumber = fNum;
    fArr = new float[iNumber];
}

float* ArrayMaker::fillArr() {
    cout << "Filling the array" << endl;
    for (int i=0; i < iNumber; i++) {
        fArr[i] = fNumber;
        fNumber *= 2;
    } 
    return fArr;
}
