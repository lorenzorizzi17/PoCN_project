#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "./include/graph.hpp"

int main() {


    int const rep = 32;
    int step = 100; // Step size for adding edges

    int types[4] = {0, 1, 2, 3}; // 0 for ER, 1 for PR, 2 for SR, 3 for BF

    omp_set_num_threads(8);
    for(int type : types){

        std::ofstream file;
        if (type == 0) {
            file.open("../data/deltaM/ER.txt");
        } else if (type == 1) {
            file.open("../data/deltaM/PR.txt");
        } else if (type == 2) {
            file.open("../data/deltaM/SR.txt");
        } else if (type == 3) {
            file.open("../data/deltaM/BF.txt");
        }

        for(int N = 1e4; N <= 1e6; N *= 2) {
            int sumDeltaM = 0; // Sum of deltaM values
            #pragma omp parallel for reduction(+:sumDeltaM)
            for (int r = 0; r < rep; r++) {
                LinkedGraph graph(N); // Create a graph with 10 nodes
                int deltaM = 0;
                bool found0 = false; // Flag to check if we found the first deltaM where LCC > 0.5*N
                bool found1 = false; // Flag to check if we found the first deltaM where LCC > 0.9*N
                for(int i = 0; i < N; i+= step) {
                    if (type == 0) {
                        graph.addRandomEdges(step); // Add a random edge for ER
                    } else if (type == 1) {
                        graph.addRandomEdgesProductRule(step); // Add edges for Product Rule
                    } else if (type == 2) {
                        graph.addRandomEdgesSumRule(step); // Add edges for Sum Rule
                    } else if (type == 3) {
                        graph.addEdgesBFRule(step); // Add edges for BF Rule
                    }
                    int lcc = graph.getLCC(); // Get the LCC size
                    if (lcc > std::sqrt(N) && !found0) {
                        found0 = true;
                        deltaM = i;
                    }
                    if (found0 && lcc > 0.5 * N && !found1) {
                        deltaM = i-deltaM;
                        found1 = true;
                        break;
                    }
                }
                sumDeltaM += deltaM; // Add the deltaM value to the sum
            }
            file << N << " " << double(sumDeltaM) / double(rep)<< std::endl;
        }
    }
    
}