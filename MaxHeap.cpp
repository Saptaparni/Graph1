//
//  Graph.cpp
//  Graph1
//
//  Created by Saptaparni Kumar on 11/27/14.
//  Copyright (c) 2014 Saptaparni Kumar. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <queue>
#include "MaxHeap.h"
using namespace std;
MaxHeap::MaxHeap(int* array, int length) : m_heap(length)
{
//    for(int i = 0; i < length; ++i)
//    {
//        m_heap[i] = array[i];
//    }
//    
//    Heapify();
}

MaxHeap::MaxHeap(const vector<pair<double,int> >& vector) : m_heap(vector)
{
    Heapify();
}

MaxHeap::MaxHeap()
{
}

void MaxHeap::Heapify()
{
    int length = m_heap.size();
    for(int i=length-1; i>=0; --i)
    {
        BubbleDown(i);
    }
}
void MaxHeap::Heapify( int index)
{
    //int length = m_heap.size();
    for(int i=index; i>=0; --i)
    {
        BubbleDown(i);
    }
}

void MaxHeap::MaxHeapify(int fromIndex, int toIndex)
{
    
    int j;
    pair<double,int> temp = m_heap[fromIndex];
    j = 2*fromIndex;
    while (j <= toIndex)
    {
        if (j < toIndex && m_heap[j+1] > m_heap[j])
            j = j+1;
        if (temp > m_heap[j])
            break;
        else if (temp <= m_heap[j])
        {
            m_heap[j/2] = m_heap[j];
            j = 2*j;
        }
    }
    m_heap[j/2] = temp;
    return;
}

void MaxHeap::BubbleDown(int index)
{
    int length = m_heap.size();
    int leftChildIndex = 2*index + 1;
    int rightChildIndex = 2*index + 2;
    
    if(leftChildIndex >= length)
        return; //index is a leaf
    
    int maxIndex = index;
    
    if(m_heap[index] < m_heap[leftChildIndex])
    {
        maxIndex = leftChildIndex;
    }
    
    if((rightChildIndex < length) && (m_heap[maxIndex] < m_heap[rightChildIndex]))
    {
        maxIndex = rightChildIndex;
    }
    
    if(maxIndex != index)
    {
        //need to swap
        pair<double,int> temp = m_heap[index];
        m_heap[index] = m_heap[maxIndex];
        m_heap[maxIndex] = temp;
        BubbleDown(maxIndex);
    }

}

void MaxHeap::BubbleUp(int index)
{
    if(index == 0)
        return;
    
    int parentIndex = (index-1)/2;
    
    if(m_heap[parentIndex] < m_heap[index])
    {
        pair<double,int> temp = m_heap[parentIndex];
        m_heap[parentIndex] = m_heap[index];
        m_heap[index] = temp;
        BubbleUp(parentIndex);
    }
}

int MaxHeap::getSize()
{
    return m_heap.size();
}

void MaxHeap::Insert(double weight, int newValue)
{
    int length = m_heap.size();
    m_heap.resize(length+1);
    m_heap[length] = make_pair(weight,newValue);
    
    BubbleUp(length);
}

pair<double,int> MaxHeap::Maximum()
{
    return m_heap[0];
}

void MaxHeap::DeleteMax()
{
    int length = m_heap.size();
    
    if(length == 0)
    {
        return;
    }
    
    m_heap[0] = m_heap[length-1];
    m_heap.pop_back();
    
    BubbleDown(0);
}

void MaxHeap::heapsort()
{
    int length = m_heap.size();
    //pair<double,int> temp;
    
    for (int i = length; i >= 1; i--)
    {
        
        Heapify(length);
    }
}

bool MaxHeap::deleteElem(int u, double dist){
//    pair<double ,pair<int, int> > val1 = make_pair(dist, make_pair(u, v));
//    pair<double ,pair<int, int> > val2 = make_pair(dist, make_pair(v, u));
//    pair<double ,pair<int, int> > delete1;
//    
    if (m_heap.size() > 0) {
        for (vector<pair<double,int> >::iterator it1 = m_heap.begin();
             it1 != m_heap.end(); it1++)
        //if (m_heap.(make_pair(dist, u)))!= m_heap.end()) {
            if (*it1 == make_pair(dist, u)) {
                m_heap.erase(it1);
                return true;
            }
        
    }
    return false;
}

void MaxHeap::PrintHeap()
{
    int length = m_heap.size();
    for (int i= 0; i<length; ++i) {
        cout<< "("<<m_heap[i].second <<", "<< m_heap[i].first<<")\n";
    }
}