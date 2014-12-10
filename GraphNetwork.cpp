//
//  GraphNetwork.cpp
//  Graph1
//
//  Created by Saptaparni Kumar on 11/27/14.
//  Copyright (c) 2014 Saptaparni Kumar. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <algorithm>
#include "MaxHeap.h"
#include "MaxEdgeSortedHeap.h"
#include "GraphNetwork.h"

using namespace std;

typedef vector<std::vector<neighbor> > adjacency_list_t;
enum status {Unseen, Fringe, Intree};
std::clock_t start;
double duration;

struct subset
{
    int parent;
    int rank;
};

int Find(int vertex, vector<int>& dad)
{
    
    if(dad[vertex]!=vertex)
    {
        dad[vertex]=Find(dad[vertex],dad);
    }
    return dad[vertex];
}

void Union(int r1, int r2,vector<int>& dad, vector<int>& rank)
{
    if(rank[r1]>rank[r2])
        dad[r1]=r2;
    else if(rank[r1]<rank[r2])
        dad[r2]=r1;
    else
    {
        dad[r1]=r2;
        rank[r1]=rank[r1]+1;
    }
}


// = false;
void DFSUtil(int v, int u, bool visited[], map<pair<int, int>, double> STree,
             list<int>& pathDFS, bool& found, double& weight)
{
    // Mark the current node as visited and print it
    //cout<<"Calling DFSUtil";
    
    visited[v] = true;
//    pathDFS.push_back(v);
//    if (v == 5) {
//        return;
//    }
    //cout << v << " ";
    int current;
    // Recur for all the vertices adjacent to this vertex
    list<int>::iterator i;
    for(map<pair<int, int>, double>::iterator it1 = STree.begin();
                it1 != STree.end(); ++it1)
            {
                if (it1->first.first == v || it1->first.second == v) {
                    //path.push_back(it1->first);
                    current = v == it1->first.first ? it1->first.second : it1->first.first;
                    if(!visited[current]){
                        DFSUtil(current, u, visited, STree, pathDFS, found, weight);
                        if (found){
                            pathDFS.push_front(v);
                            weight = min(weight, it1->second);
                            return;
                        }
                    }
                    
                }
            }
    if (v == u) {
        pathDFS.push_front(v);
        found =true;
        return;
    }
}

// DFS traversal of the vertices reachable from v. It uses recursive DFSUtil()
void Graph::DFS(int v, int u,const map<pair<int, int>,
                double>& STree, list<int>& pathDFS, double& weight)
{
    
    // Mark all the vertices as unvisited
    bool *visited = new bool[V];
    for(int i = 0; i < V; i++)
        visited[i] = false;
    bool found = false;
    DFSUtil(v, u, visited, STree, pathDFS, found, weight);
}



void Graph::CreateGraph(int degree)
{
    long tot = V*degree;
    set<int> edge ;
    for (int i=0; i<tot; i++) edge.insert(i);
    //cout<<"********************Creating Spanning Tree**********************\n";
    set<int> unvisited;
    for (int i=0; i<V; i++) unvisited.insert(i);
    map<int,int> visited;
    srand(time(NULL));
    int visit = rand() % V ;
    visited.insert(make_pair(visit,0));
    unvisited.erase(visit);
    while (!unvisited.empty()) {
        //choose a vertex from visited
        visit = rand() % V ;
        while (visited.find(visit)== visited.end() || visited.find(visit)->second >degree-1) {
            visit = rand() % V ;
        }
        int unvisit = rand() % V ;
        while( unvisited.find(unvisit) == unvisited.end()
              || Edge.find(make_pair(visit,unvisit)) != Edge.end()
              || Edge.find(make_pair(unvisit,visit)) != Edge.end())
        {
            unvisit = rand() % V ;
        }
        //now we have visit and unvisit
        int iter1= visit*degree;
        while (edge.find(iter1) == edge.end() && iter1 < (visit+1)*degree) {
            iter1++;
        }
        
        edge.erase(edge.find(iter1));
        int iter2= unvisit*degree;
        while (edge.find(iter2) == edge.end() && iter2 < (unvisit+1)*degree) {
            iter2++;
        }
        edge.erase(edge.find(iter2));
        (visited.find(visit)->second)++;
        visited.insert(make_pair(unvisit, 1));
        unvisited.erase(unvisit);
        addEdge(visit,unvisit);
        pair<int,int> edge_12 = make_pair(visit,unvisit);
        Edge.insert(edge_12);
        int weight = rand() % V*2+1 ;
        EdgeWeights.insert(make_pair(edge_12,weight));
    }
    
    //***********Add rest of the edges*********************
    while(!edge.empty())
    {
        int vert_11 = rand() % tot ;
        while( edge.find(vert_11) == edge.end())
            vert_11 = rand() % tot ;
        edge.erase(edge.find(vert_11));
        int vert_12 = vert_11/degree;
        int vert_21 = rand() % tot ;
        int vert_22 = vert_21/degree;
        while( edge.find(vert_21) == edge.end() || vert_22 == vert_12
              || Edge.find(make_pair(vert_12,vert_22)) != Edge.end()
              || Edge.find(make_pair(vert_22,vert_12)) != Edge.end())
        {
            vert_21 = rand() % tot ;
            vert_22 = vert_21/degree;
            
        }
        edge.erase(edge.find(vert_21));
        addEdge(vert_12, vert_22);
        pair<int,int> edge_12 = make_pair(vert_12,vert_22);
        Edge.insert(edge_12);
        int weight = rand() % V*2 +1;
        EdgeWeights.insert(make_pair(edge_12,weight));
    }
}


Graph::Graph(int ver, int degree)
{
    V = ver;
    adj = new list<int>[V];
    CreateGraph(degree);
}
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

void Graph::PrintGraph(string fileName)
{
    
    ofstream textfile;
    textfile.open (fileName);
    for(map<pair<int, int>,int>::iterator it1 = EdgeWeights.begin(); it1!= EdgeWeights.end();++it1)
    {
        textfile<< it1->first.first <<" to " << it1->first.second << ": weight " << it1->second<<",\t";
    }
    textfile.close();
}




void Graph::DijkstraComputeMaxCapacityPathsWithoutHeap(int source, int destination,
                                                       std::vector<double>& distance,
                                                       std::vector<int>& dad)
{
    start = std::clock();
    set<pair<double,int> > fringes;
    vector<int> status;
    //Step 1
    status.resize(V, Unseen);
    distance.resize(V, min_weight);
    dad.resize(V,-2);
    //Step 2
    status[source] = Intree;
    dad[source] = -1;
    distance[source] = 0;
    //Step 3
    for (list<int>::iterator it = adj[source].begin(); it != adj[source].end(); it++)
    {
        int neighb = *it;
        status[neighb] = Fringe;
        dad[neighb] = source;
        map<pair<int, int>,int>::iterator it1 = EdgeWeights.find(make_pair(source, neighb));
        if (it1 == EdgeWeights.end()) {
            it1 = EdgeWeights.find(make_pair(neighb,source));
        }
        if (it1 != EdgeWeights.end())
            distance[neighb] = it1->second;
        else
            cout<< "\n error \n";
        fringes.insert(make_pair(distance[neighb], neighb));
    }
    
    //Step 4
    // Add fringes to the queue here.
    while (!fringes.empty()) {
        double dist = fringes.rbegin()->first;
        int v = fringes.rbegin()->second;// get the largest
        status[v] = Intree;
        fringes.erase(make_pair(dist, v));
        //for each edge [v,w] do:
        for (list<int>::iterator it = adj[v].begin(); it != adj[v].end(); it++)
        {
            //preprocessing:
            int w = *it;
            double weight = -1;
            map<pair<int, int>,int>::iterator it1 = EdgeWeights.find(make_pair(v, w));
            if (it1 == EdgeWeights.end()) {
                it1 = EdgeWeights.find(make_pair(w,v));
            }
            if (it1 != EdgeWeights.end())
                weight = it1->second;
            else
                cout<< "\n error \n";
            if (status[w] == Unseen) {
                status[w] = Fringe;
                dad[w] = v;
                distance[w] = std::min(distance[v], weight);
                fringes.insert(make_pair(distance[w],w));
            }
            else if(status[w] == Fringe &&
                    distance[w] < min(distance[v], weight))
            {
                fringes.erase(make_pair(distance[w], w));
                distance[w] = min(distance[v], weight);
                dad[w] = v;
                fringes.insert(make_pair(distance[w], w));
            }
        }
        
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"\nDijkstra: Duration without heap : "<< duration <<'\n';
    
    std::cout << "Dijkstra: Distance without heap : " << distance[destination] << std::endl;
    std::list<int> path = DijkstraGetShortestPathTo(destination, dad);
    std::cout << "Path : ";
    std::copy(path.begin(), path.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
}



void Graph::DijkstraComputeMaxCapacityPathsWithHeap(int source, int destination,
                                                    std::vector<double>& distance,
                                                    std::vector<int>& dad)
{
    start = std::clock();
    MaxHeap fringes;
    vector<int> status;
    //Step 1
    status.resize(V, Unseen);
    distance.resize(V, min_weight);
    dad.resize(V,-2);
    //Step 2
    status[source] = Intree;
    dad[source] = -1;
    distance[source] = 0;
    //Step 3
    for (list<int>::iterator it = adj[source].begin(); it != adj[source].end(); it++)
    {
        int neighb = *it;
        status[neighb] = Fringe;
        dad[neighb] = source;
        map<pair<int, int>,int>::iterator it1 = EdgeWeights.find(make_pair(source, neighb));
        if (it1 == EdgeWeights.end()) {
            it1 = EdgeWeights.find(make_pair(neighb,source));
        }
        if (it1 != EdgeWeights.end())
            distance[neighb] = it1->second;
        else
            cout<< "\n error \n";
        fringes.Insert(distance[neighb], neighb);
    }
    
    //Step 4
    // Add fringes to the heap here.
    while (fringes.getSize() >0) {
        pair<double, int> max = fringes.Maximum();
        //double dist = max.first;
        int v = max.second;// get the largest
        status[v] = Intree;
        fringes.DeleteMax();
        //for each edge [v,w] do:
        for (list<int>::iterator it = adj[v].begin(); it != adj[v].end(); it++)
        {
            //preprocessing:
            int w = *it;
            double weight = -1;
            map<pair<int, int>,int>::iterator it1 = EdgeWeights.find(make_pair(v, w));
            if (it1 == EdgeWeights.end()) {
                it1 = EdgeWeights.find(make_pair(w,v));
            }
            if (it1 != EdgeWeights.end())
                weight = it1->second;
            else
                cout<< "\n error \n";
            if (status[w] == Unseen) {
                status[w] = Fringe;
                dad[w] = v;
                distance[w] = std::min(distance[v], weight);
                fringes.Insert(distance[w],w);
            }
            else if(status[w] == Fringe &&
                    distance[w] < min(distance[v], weight))
            {
                bool success = fringes.deleteElem( w, distance[w]);
                distance[w] = min(distance[v], weight);
                dad[w] = v;
                if (success) {
                    fringes.Insert(distance[w], w);
                }
                else
                    cout<< "\nError Dijkstra Heap";
                
                
            }
        }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"\nDijkstra: Duration with heap : "<< duration <<'\n';
    
    std::cout << "Dijkstra: Distance with heap : " << distance[destination] << std::endl;
    std::list<int> path = DijkstraGetShortestPathTo(destination, dad);
    std::cout << "Path : ";
    std::copy(path.begin(), path.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
}




std::list<int> Graph::DijkstraGetShortestPathTo(int vertex, const std::vector<int> &previous)
{
    std::list<int> path;
    for ( ; vertex != -1; vertex = previous[vertex])
        path.push_front(vertex);
    return path;
}


void Graph::KruskalComputeMaxCapacityPathsWithHeap (int source, int destination,
                                     string fileName)
{
    ofstream textfile;
    textfile.open (fileName);
    start = std::clock();
    //Empty Spanning Tree
    map<pair<int, int>, double> STree;
    //Makeset for every vertex
    vector<int> dad(V);
    vector<int> rank(V,0);
    for (int i=0; i<V; i++) {
        dad[i] = i;
    }
    //Sort Edges using heapSort
    MaxEdgeSortedHeap SortedEdges;
    for(map<pair<int, int>,int>::iterator it1 = EdgeWeights.begin(); it1!= EdgeWeights.end();++it1)
    {
        SortedEdges.Insert(make_pair(it1->second,make_pair(it1->first.first,
                                                           it1->first.second)));
    }
    int dad1, dad2;
    double dist;
    int u=-1,v=-1;
    while (SortedEdges.Size() >0) {
        
        SortedEdges.getMax(u, v, dist);
        
        dad1 = Find(u, dad);
        dad2 = Find(v, dad);
        if (dad1 != dad2) {
            pair<int,int> edge1 = u<v ? make_pair(u, v): make_pair(v, u);
            pair<pair<int,int>,double> node = make_pair(edge1,dist);
            pair<map<pair<int,int>,double>::iterator,bool> ret ;
            ret = STree.insert(node);
            if (ret.second == false) {
                cout << "error\n";
            }
        }
        Union(dad1, dad2, dad, rank);
    }
    
    list<int> pathDFS;
    double weight = max_weight;
    DFS(source, destination, STree, pathDFS, weight);
    cout<<"\nKruskal: Duration with heapSort : "<< duration <<'\n';
    cout<<"Kruskal: Distance with heapSort :"<< weight << " \nPath : ";
    for (list<int>::iterator it1 = pathDFS.begin(); it1 != pathDFS.end(); it1++) {
        cout<< *it1 << " ";
    }
    
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    textfile<<"\nKruskal using HeapSort\n";
    for(map<pair<int, int>, double>::iterator it1 = STree.begin();
        it1 != STree.end(); ++it1)
        textfile << "("<<it1->first.first<<","<< it1->first.second<<"): ";//<<it1->second<<"   ";
    textfile.close();
}






