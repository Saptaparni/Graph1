//
//  main.cpp
//  Graph1
//
//  Created by Saptaparni Kumar on 11/25/14.
//  Copyright (c) 2014 Saptaparni Kumar. All rights reserved.
//

#include <iostream>
#include <vector>
#include <time.h>
#include <fstream>
#include "GraphNetwork.h"
#include "MaxHeap.h"
#include "MaxEdgeSortedHeap.h"

int main()
{
    std::clock_t start;
    double duration;
    
    start = std::clock();
    Graph g1(5000,6);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Sparse Graph Generation duration: "<< duration <<'\n';
    
    //g1.PrintGraph("/Users/saptaparnikumar/Desktop/XCode/Graph1/Graph1/test1.txt");
    g1.PrintGraph("sparse1.txt");
    
    
    start = std::clock();
    Graph g2(5000,1000);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Dense Graph Generation duration: "<< duration <<'\n';
    
    //g2.PrintGraph("/Users/saptaparnikumar/Desktop/XCode/Graph1/Graph1/test2.txt");
    g2.PrintGraph("sparse2.txt");
   
    //Testing for no heap dijkstra
    std::vector<double> max_distance;
    std::vector<int> previous;
    
    //for (int i=0; i<5; i++) {
       
    srand(time(NULL));
    int u = rand() % 5000 ;
    int v = rand() % 5000 ;
    //cout<<"\n\n***************************Algorithm Set "<<i << " ******************************* \n";
    cout<<"\n#############################Sparse Graph##########################\n";
    //Sparse Graph: Testing for no heap dijkstra
    g1.DijkstraComputeMaxCapacityPathsWithoutHeap(u, v, max_distance, previous);
    //Sparse Graph: Testing for heap + dijkstra
    g1.DijkstraComputeMaxCapacityPathsWithHeap(u, v, max_distance, previous);
    //Sparse Graph: Testing for kruskal's using heapsort
    //For the values, check file Kruskal1.txt
    //map<pair<int, int>, double> Tree1 =
    g1.KruskalComputeMaxCapacityPathsWithHeap(u, v, "Kruskal1.txt");
    
    cout<<"\n\n#############################Dense Graph##########################\n";
    //Dense Graph: Testing for no heap dijkstra
    g2.DijkstraComputeMaxCapacityPathsWithoutHeap(u, v, max_distance, previous);
    //Dense Graph: Testing for heap dijkstra
    g2.DijkstraComputeMaxCapacityPathsWithHeap(u, v, max_distance, previous);
    //Dense Graph: Testing for kruskal's using heapsort
    //For the values, check file Kruskal2.txt
    //map<pair<int, int>, double> Tree2 =
    g2.KruskalComputeMaxCapacityPathsWithHeap(u, v, "Kruskal2.txt");
  //  }

    return 0;
}
