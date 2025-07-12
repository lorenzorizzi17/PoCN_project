#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <numeric>

#include "include/graph.hpp"


int main() {
    int N = 1e6;
    int rep = 10;

    int mCritical = 5*1e5;
    int mValues[] = {mCritical - N/5, mCritical - N/10, mCritical, mCritical + N/10, mCritical + N/5};

    int type = 0;

    for (int m : mValues) {
        std::map<int, std::vector<int>> cluster_counts; 
        
        for (int r = 0; r < rep; ++r) {
            LinkedGraph graph(N);
            switch (type) {
                case 0:
                    graph.addRandomEdges(m);
                    break;
                case 1:
                    graph.addRandomEdgesProductRule(m);
                    break;
                case 2: 
                    graph.addRandomEdgesSumRule(m);
                    break;
                case 3: 
                    graph.addEdgesBFRule(m);
                    break;
                default:
                    std::cerr << "Unknown type!" << std::endl;
                    return EXIT_FAILURE;
            }
            std::vector<int> dist = graph.getClusterDistribution();

            std::map<int, int> count_per_size;
            for (int size : dist)
                count_per_size[size]++;

            for (const auto& [size, count] : count_per_size)
                cluster_counts[size].push_back(count);

        }


        std::ofstream file("../data/ClusterDistribution/ER_" + std::to_string(m) + ".txt");
        for (const auto& [size, counts] : cluster_counts) {
            double avg = std::accumulate(counts.begin(), counts.end(), 0.0) / counts.size();
            file << size << " " << avg << std::endl;
        }
        file.close();
    }
}
