#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <algorithm> 

struct Node {
    Node* next = NULL; // Pointer to the next node in the graph. If null, then root node
    std::vector<Node*> neighbors; // A vector of neighbors
    int ClusterSize = 1;  // If root, store the size of the cluster. Default to 0.
    Node() = default; // Default constructor
};

class LinkedGraph{
    private:
    int N; // Number of nodes
    std::vector<Node> nodes;   // A vector of nodes

    public:
    // Void constructor (no edges, every node is a root node)
    LinkedGraph(int N);

    // Sampler from a scale-free distribution
    void RandomPercolationSF(int N, double alpha, int kmin, int kmax);

    std::vector<int> getClusterDistribution() const;

    std::vector<int> getDegreeDistribution() const;

    int getLCC();

    void addEdgeSFPR(std::vector<int>&);

    void addEdgesSFPR(int n, std::vector<int>&);
    
    void addEdgeSF(std::vector<int>&);

    void addEdgesSF(int n, std::vector<int>&);


    double getAverageClusterSize() const;

    int getSecondLargestClusterSize() const;

    Node* findRoot(Node* node);

    void addRandomEdge();

    void addRandomEdges(int numEdges);

    void addRandomProductRule();

    void addRandomEdgesProductRule(int);

    void addRandomSumRule();

    void addRandomEdgesSumRule(int);

    void addEdgeBFRule();

    void addEdgesBFRule(int);
};