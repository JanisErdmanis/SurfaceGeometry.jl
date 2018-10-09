#ifndef ARRAYMAKER_H
#define ARRAYMAKER_H

class ArrayMaker
{
    private:
        int iNumber;
        float fNumber;
        float* fArr;
    public:
        ArrayMaker(int, float);
        float* fillArr();
};

#endif
