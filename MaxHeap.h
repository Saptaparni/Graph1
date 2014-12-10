//
//  MaxHeap1.h
//  Graph1
//
//  Created by Saptaparni Kumar on 11/27/14.
//  Copyright (c) 2014 Saptaparni Kumar. All rights reserved.
//

#ifndef Graph1_Heap1_h
#define Graph1_Heap1_h
#include "vector"
using namespace std;

class MaxHeap
{
private:
    vector<pair<double,int> > m_heap;
    void BubbleDown(int index);
    void BubbleUp(int index);
    void Heapify();
    void Heapify(int index);
    
public:
    MaxHeap(int* array, int length);
    MaxHeap(const vector<pair<double,int> >& p_heap);
    MaxHeap();
    int getSize();
    
    void Insert(double weight, int newValue);
    pair<double, int> Maximum();
    void DeleteMax();
    void heapsort();
    void MaxHeapify(int fromIndex, int toIndex);
    void PrintHeap();
    bool deleteElem(int u, double dist);
    
};

#endif
