#include<iostream>
#include <list>
#include <set>
#include <ctime>
#include <cstdlib>
#include <utility>
#include <map>
#include <limits> // for numeric_limits
#include <algorithm>
#include <iterator>
#include <vector>
#include "MaxHeap.h"
using namespace std;

const double max_weight = std::numeric_limits<double>::infinity();
const double min_weight = 0-std::numeric_limits<double>::infinity();
struct neighbor {
    int target;
    double weight;
    neighbor(int arg_target, double arg_weight)
    : target(arg_target), weight(arg_weight) { }
};

// A class that represents an undirected graph
class Graph
{
    int V;    // No. of vertices
    list<int> *adj;    // A dynamic array of adjacency lists
    set<pair<int, int> > Edge; //Set of Edges
    map<pair<int, int>,int> EdgeWeights;// Their corresponding weights
    void CreateGraph(int degree);
    
public:
    // Constructor and destructor
    Graph(int V)   {this->V = V; adj = new list<int>[V]; }
    Graph(int ver, int degree);
    ~Graph() { delete [] adj; } // To avoid memory leak
    
    // function to add an edge to graph
    void addEdge(int v, int w);
    
    // Method to check if all non-zero degree vertices are connected
    bool isConnected();
    void DFS(int v, int u, const map<pair<int, int>, double>& STree,
             list<int>& pathDFS, double& weight);
    
    // Function to do DFS starting from v. Used in isConnected();
    //void DFSUtil(int v, bool visited[]);
    void PrintGraph(string fileName);
    void DijkstraComputeMaxCapacityPathsWithoutHeap(int source, int destination,
                                                    std::vector<double> &distance,
                                                    std::vector<int> &dad);
    
    void DijkstraComputeMaxCapacityPathsWithHeap(int source, int destination,
                                                 std::vector<double>& distance,
                                                 std::vector<int>& dad);
    std::list<int> DijkstraGetShortestPathTo(int vertex,
                                             const std::vector<int> &previous);
    void KruskalComputeMaxCapacityPathsWithHeap(int source, int destination, string fileName);
};
