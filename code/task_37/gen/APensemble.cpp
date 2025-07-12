#include <fstream>

#include "include/graph.hpp"

int main(){
    int N = 1e6;
    int m = N/2;

    for(int i = 0; i < 1; i++){
        LinkedGraph graph(N); // Create a new graph with N nodes
        graph.addRandomEdgesProductRule(m); // Add m random edges
        std::vector<int> dd = graph.getDegreeDistribution();
        std::ofstream file("./data/DegreeDistribution/ER.txt");
        for (int j = 0; j < dd.size(); j++) {
            file << dd[j] << std::endl;
        }
        file.close();
    }
}