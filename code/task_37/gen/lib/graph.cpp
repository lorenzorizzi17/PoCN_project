#include "../include/graph.hpp"
#include<map>
#include<assert.h>
#include<numeric>
#include<algorithm>

LinkedGraph::LinkedGraph(int n) : N(n) {
    nodes.reserve(N); // Reserve space for N nodes
    for (int i = 0; i < N; i++){
        nodes.push_back(Node());
    }
}


std::vector<int> LinkedGraph::getClusterDistribution() const {
    std::vector<int> result;
    for (const auto& node : nodes) {
        if (node.next == nullptr) {
            result.push_back(node.ClusterSize); // If it's a root node, add its cluster size
        }
    }
    // discard the largest one
    auto it = std::max_element(result.begin(), result.end(), [](int a, int b) { return a < b; });
    return result;
}

double LinkedGraph::getAverageClusterSize() const {
    // first of all, get the cluster distribution
    std::vector<int> clusterSizes = getClusterDistribution();
    // build a map to count occurrences of each cluster size (like table() in R)
    std::map<int, int> occurrences;
    for (int cl : clusterSizes) {
        occurrences[cl]++;
    }
    // Avoid the LCC, clearly
    double numerator=0 ; double denominator = 0;
    for(auto it = occurrences.begin(); it != std::prev(occurrences.end()); ++it) {
        numerator += it->first* it->first* it->second; 
        denominator += it->second*it->first; 
    }
    return numerator / denominator;
}

int LinkedGraph::getSecondLargestClusterSize() const {
    std::vector<int> clusterSizes = getClusterDistribution();
    if (clusterSizes.size() < 2) {
        return 0; // Not enough clusters to find the second largest
    }
    std::sort(clusterSizes.begin(), clusterSizes.end(), std::greater<int>());
    return clusterSizes[1]; // Return the second largest cluster size
}

int LinkedGraph::getLCC(){
    std::vector<int> clusterSizes = getClusterDistribution();
    return *max_element(clusterSizes.begin(), clusterSizes.end(), [](int a, int b) { return a < b;});
}

Node* LinkedGraph::findRoot(Node* node) {
    if (node->next == nullptr)
        return node;

    // Path compression
    node->next = findRoot(node->next);
    return node->next;
}
void LinkedGraph::addEdgeSF(std::vector<int>& stubs) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, stubs.size() - 1);
    if (stubs.size() < 2) {
        return;
    }
    int i1, i2;
    do {
        i1 = dis(gen);
        i2 = dis(gen);
    } while (i1 == i2);

    auto it1 = stubs.begin();
    std::advance(it1, i1);

    auto it2 = stubs.begin();
    std::advance(it2, i2);

    if (it2 < it1) {
        std::swap(it1, it2);
    }
    int u = *it1;
    int v = *it2;

    Node* rootU = findRoot(&nodes[u]);  // Find the root of node u
    Node* rootV = findRoot(&nodes[v]);  // Find the root of node v
    int S1 = rootU->ClusterSize; // Size of the cluster of u
    int S2 = rootV->ClusterSize; // Size of the cluster of v

    nodes[u].neighbors.push_back(&nodes[v]); // Add v as a neighbor of u
    nodes[v].neighbors.push_back(&nodes[u]); // Add u as a neighbor of v

    if (rootU != rootV) { // If they are not in the same cluster
        // Merge the two clusters
        if (S1 < S2) {
            rootU->next = rootV; // Point u to v
            rootV->ClusterSize += S1; // Update the size of the cluster of v
        } else {
            rootV->next = rootU; // Point v to u
            rootU->ClusterSize += S2; // Update the size of the cluster of u
        }
    }

    stubs.erase(it2);  
    stubs.erase(it1);  
}

void LinkedGraph::addEdgeSFPR(std::vector<int>& stubs) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, stubs.size() - 1);

    if (stubs.size() < 4) {
        return;
    }
    // choose 4 stubs randomly
    int i1, i2, i3, i4;
    do {
        i1 = dis(gen);
        i2 = dis(gen);
        i3 = dis(gen);
        i4 = dis(gen);
    } while (i1 == i2 || i1 == i3 || i1 == i4 || i2 == i3 || i2 == i4 || i3 == i4);

    // Get the iterators for the stubs
    auto it1 = stubs.begin();
    std::advance(it1, i1);
    auto it2 = stubs.begin();
    std::advance(it2, i2);
    auto it3 = stubs.begin();
    std::advance(it3, i3);
    auto it4 = stubs.begin();
    std::advance(it4, i4);

    // Ensure that the stubs are in increasing order (so that erasing them does not mess things up)
    if (it2 < it1) {
        std::swap(it1, it2);
    }
    if (it4 < it3) {
        std::swap(it3, it4);
    }

    int stub1 = *it1, stub2 = *it2;
    int stub3 = *it3, stub4 = *it4;

    // Find the roots of all 4 candidate nodes ...
    Node* rootStub1 = findRoot(&nodes[stub1]);
    Node* rootStub2 = findRoot(&nodes[stub2]);
    Node* rootStub3 = findRoot(&nodes[stub3]);
    Node* rootStub4 = findRoot(&nodes[stub4]);
    // ... and their cluster size
    int S1 = rootStub1->ClusterSize;
    int S2 = rootStub2->ClusterSize;
    int S3 = rootStub3->ClusterSize;
    int S4 = rootStub4->ClusterSize;
    // product rule
    if (S1 * S2 < S3 * S4 || (S1 * S2 == S3 * S4 && dis(gen) % 2 == 0)) { //choose stub1 and stub2
        nodes[stub1].neighbors.push_back(&nodes[stub2]);
        nodes[stub2].neighbors.push_back(&nodes[stub1]);

        if(rootStub1 != rootStub2) { // If they are not in the same cluster
            // Merge the two clusters
            if (S1 < S2) {
                rootStub1->next = rootStub2; // Point stub1 to stub2
                rootStub2->ClusterSize += S1; // Update the size of the cluster of stub2
            } else {
                rootStub2->next = rootStub1; // Point stub2 to stub1
                rootStub1->ClusterSize += S2; // Update the size of the cluster of stub1
            }
        }
        // Erase the stubs that have been used
        stubs.erase(it2);  
        stubs.erase(it1);  
    } else {
        nodes[stub3].neighbors.push_back(&nodes[stub4]);
        nodes[stub4].neighbors.push_back(&nodes[stub3]);

        if(rootStub3 != rootStub4) {
            // Merge the two clusters
            if (S3 < S4) {
                rootStub3->next = rootStub4; // Point stub3 to stub4
                rootStub4->ClusterSize += S3; // Update the size of the cluster of stub4
            } else {
                rootStub4->next = rootStub3; // Point stub4 to stub3
                rootStub3->ClusterSize += S4; // Update the size of the cluster of stub3
            }
        }

        stubs.erase(it4);  
        stubs.erase(it3);  
    }
}

void LinkedGraph::addRandomEdge() {
    // first of all, select two random nodes
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, N - 1);
    int u = dis(gen);
    int v = dis(gen);

    Node* rootU = findRoot(&nodes[u]);  // Find the root of node u
    Node* rootV = findRoot(&nodes[v]);  // Find the root of node v
    int S1 = rootU->ClusterSize; // Size of the cluster of u
    int S2 = rootV->ClusterSize; // Size of the cluster of v

    nodes[u].neighbors.push_back(&nodes[v]); // Add v as a neighbor of u
    nodes[v].neighbors.push_back(&nodes[u]); // Add u as a neighbor of v

    if (rootU != rootV) { // If they are not in the same cluster
        // Merge the two clusters
        if (S1 < S2) {
            rootU->next = rootV; // Point u to v
            rootV->ClusterSize += S1; // Update the size of the cluster of v
        } else {
            rootV->next = rootU; // Point v to u
            rootU->ClusterSize += S2; // Update the size of the cluster of u
        }
    }
}





void LinkedGraph::addEdgesSF(int n, std::vector<int>& stubs) {
    for (int i = 0; i < n; ++i) {
        this->addEdgeSF(stubs);
    }
}

void LinkedGraph::addEdgesSFPR(int n, std::vector<int>& stubs) {
    for (int i = 0; i < n; ++i) {
        this->addEdgeSFPR(stubs);
    }
}

std::vector<int> LinkedGraph::getDegreeDistribution() const {
    std::vector<int> degreeDistribution;
    for (const auto& node : nodes) {
        int k = node.neighbors.size(); // Get the degree of the node
        degreeDistribution.push_back(k); // Increment the count for this cluster size
    }
    return degreeDistribution;
}

void LinkedGraph::addRandomEdges(int numEdges) {
    for (int i = 0; i < numEdges; ++i) {
        addRandomEdge();
    }
}

void LinkedGraph::addRandomProductRule(){
    // Select two random nodes
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, N - 1);

    // first pair (one edge)
    int u1 = dis(gen);
    int u2 = dis(gen);
    // second pair (one edge)
    int v1 = dis(gen);
    int v2 = dis(gen);

    Node* rootU1 = findRoot(&nodes[u1]);  // Find the root of node u1
    Node* rootU2 = findRoot(&nodes[u2]);  // Find the root of node u2
    int s1 = rootU1->ClusterSize; // Size of the cluster of u1
    int s2 = rootU2->ClusterSize; // Size of the cluster of u2
    int product1 = s1 * s2; // Product of the sizes of the clusters

    Node* rootV1 = findRoot(&nodes[v1]);  // Find the root of node u1
    Node* rootV2 = findRoot(&nodes[v2]);  // Find the root of node u2
    int s1v = rootV1->ClusterSize; // Size of the cluster of u1
    int s2v = rootV2->ClusterSize; // Size of the cluster of u2
    int product2 = s1v * s2v; // Product of the sizes of the clusters

    if (product1 < product2) {
        nodes[u1].neighbors.push_back(&nodes[u2]); // Add v as a neighbor of u
        nodes[u2].neighbors.push_back(&nodes[u1]); // Add u as a neighbor of v
        // Select the first edge, u1-u2 and discard the first edge
        if (rootU1 != rootU2) { // If they are not in the same cluster
            rootU1->next = rootU2; // Point u1 to u2
            rootU2->ClusterSize += s1; // Update the size of the cluster of u2
        }
    } else {
        nodes[v1].neighbors.push_back(&nodes[v2]); // Add v as a neighbor of u
        nodes[v2].neighbors.push_back(&nodes[v1]); // Add u as a neighbor of v
        // Select the second edge, v1-v2 and discard the first edge
        if (rootV1 != rootV2) { // If they are not in the same cluster
            rootV1->next = rootV2; // Point v1 to v2
            rootV2->ClusterSize += s1v; // Update the size of the cluster of v2
        }
    }
}

void LinkedGraph::addRandomEdgesProductRule(int numEdges) {
    for (int i = 0; i < numEdges; ++i) {
        this->addRandomProductRule();
    }
}

void LinkedGraph::addRandomSumRule(){
    // Select two random nodes
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, N - 1);

    // first pair (one edge)
    int u1 = dis(gen);
    int u2 = dis(gen);
    // second pair (one edge)
    int v1 = dis(gen);
    int v2 = dis(gen);

    Node* rootU1 = findRoot(&nodes[u1]);  // Find the root of node u1
    Node* rootU2 = findRoot(&nodes[u2]);  // Find the root of node u2
    int s1 = rootU1->ClusterSize; // Size of the cluster of u1
    int s2 = rootU2->ClusterSize; // Size of the cluster of u2
    int product1 = s1 + s2; // Product of the sizes of the clusters

    Node* rootV1 = findRoot(&nodes[v1]);  // Find the root of node u1
    Node* rootV2 = findRoot(&nodes[v2]);  // Find the root of node u2
    int s1v = rootV1->ClusterSize; // Size of the cluster of u1
    int s2v = rootV2->ClusterSize; // Size of the cluster of u2
    int product2 = s1v + s2v; // Product of the sizes of the clusters

    if (product1 < product2) {
        // Select the first edge, u1-u2 and discard the first edge
        nodes[u1].neighbors.push_back(&nodes[u2]); // Add v as a neighbor of u
        nodes[u2].neighbors.push_back(&nodes[u1]); // Add u as a neighbor of v
        if (rootU1 != rootU2) { // If they are not in the same cluster
            rootU1->next = rootU2; // Point u1 to u2
            rootU2->ClusterSize += s1; // Update the size of the cluster of u2
        }
    } else {
        // Select the second edge, v1-v2 and discard the first edge
        nodes[v1].neighbors.push_back(&nodes[v2]); // Add v as a neighbor of u
        nodes[v2].neighbors.push_back(&nodes[v1]); // Add u as a neighbor of v
        if (rootV1 != rootV2) { // If they are not in the same cluster
            rootV1->next = rootV2; // Point v1 to v2
            rootV2->ClusterSize += s1v; // Update the size of the cluster of v2
        }
    }
}

void LinkedGraph::addRandomEdgesSumRule(int numEdges) {
    for (int i = 0; i < numEdges; ++i) {
        this->addRandomSumRule();
    }
}


void LinkedGraph::addEdgeBFRule(){
    // Select two random nodes
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, N - 1);

    // first pair (one edge)
    int u1 = dis(gen);
    int u2 = dis(gen);

    if (nodes[u1].next == nullptr && nodes[u2].next == nullptr && nodes[u2].ClusterSize == 1 && nodes[u1].ClusterSize == 1 && u1 != u2) { //so both of them are root nodes
        nodes[u1].neighbors.push_back(&nodes[u2]); // Add v as a neighbor of u
        nodes[u2].neighbors.push_back(&nodes[u1]); // Add u as a neighbor of v
        nodes[u1].next = &nodes[u2]; // Point u1 to u2
        nodes[u2].ClusterSize += 1; // Update the size of the cluster of u2
        return;
    } // Otherwise, find a new pair and connect it

    this->addRandomEdge();
}

void LinkedGraph::addEdgesBFRule(int numEdges) {
    for (int i = 0; i < numEdges; ++i) {
        this->addEdgeBFRule();
    }
}

