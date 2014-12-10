//
//  MaxEdgeHeap.h
//  Graph1
//
//  Created by Saptaparni Kumar on 11/30/14.
//  Copyright (c) 2014 Saptaparni Kumar. All rights reserved.
//

#ifndef Graph1_MaxEdgeHeap_h
#define Graph1_MaxEdgeHeap_h
#include <set>
#include <iostream>
using namespace std;

class MaxEdgeSortedHeap
{
    
public:
    set<pair<double ,pair<int, int> > > SortedEdges;
//    MaxHeap(int* array, int length);
//    MaxHeap(const vector<pair<double,int> >& p_heap);
//    MaxHeap();
//    int getSize();
//    
//    
    void Insert(pair<double ,pair<int, int > >  newValue)
    {
        SortedEdges.insert(newValue);
    }
    int Size()
    {
        return SortedEdges.size();
    }
    
    void getMax(int& u, int& v, double& dist)
    {
        if (SortedEdges.size() > 0) {
            dist = SortedEdges.rbegin()->first;
            u = SortedEdges.rbegin()->second.first;
            v = SortedEdges.rbegin()->second.second;
            SortedEdges.erase(make_pair(dist, make_pair(u, v)));
        }
    }
    
//    pair<double, int> Maximum();
//    void DeleteMax();
//    void heapsort();
//    void PrintHeap();
    
};



#endif
