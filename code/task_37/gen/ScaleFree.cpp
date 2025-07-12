#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "./include/graph.hpp"

std::vector<int> generatePowerLawDegrees(int N, double alpha, int kmin, int kmax) {
    std::vector<double> cdf;
    double norm = 0.0;

    // Build the unnormalized CDF
    for (int k = kmin; k <= kmax; ++k) {
        norm += std::pow(k, -alpha);
        cdf.push_back(norm);
    }

    // Normalize the CDF
    for (double& val : cdf) val /= norm;

    // Setup RNG
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Prepare to generate degrees
    std::vector<int> degrees;
    degrees.reserve(N);

    // Generate degrees based on the CDF using inverse transform sampling
    for (int i = 0; i < N; ++i) {
        double r = dis(gen);
        auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
        int k = std::distance(cdf.begin(), it) + kmin;
        degrees.push_back(k);
    }

    return degrees;
}


std::vector<int> generateStubs(const std::vector<int>& degrees) {
    std::vector<int> stubs;
    for (int i = 0; i < degrees.size(); ++i) {
        stubs.insert(stubs.end(), degrees[i], i); 
    }
    return stubs;
}


double kappa_from_degree_sequence(const std::vector<int>& degrees) {
    if (degrees.empty()) return 0.0;

    double sum_k = 0.0;
    double sum_k2 = 0.0;

    for (int k : degrees) {
        sum_k += k;
        sum_k2 += static_cast<double>(k) * k;
    }

    return sum_k2 / sum_k;
}

int main(){
    int N = 1e6; // Number of nodes
    int rep = 4;
    int steps = 150;

    int types[1] = {0}; // 0 for RG, 1 for PR
    double alphas[1] = {3.8}; // Different alpha values
    int kmins[1] = {1}; // Different kmin values

    for(int type : types){  // 0 for RG, 1 for PR
        for (double alpha : alphas) { // Different alpha values
            for (int kmin : kmins) { // Different kmin values

                int kmax = kmin*static_cast<int>(std::pow(N, 1.0 / (alpha - 1)));
                std::vector<int> degrees = generatePowerLawDegrees(N, alpha, kmin, kmax);
                std::vector<int> stubs = generateStubs(degrees);

                std::cout << "Working with type = " << type << ", alpha = " << alpha << ", kmin = " << kmin << std::endl;
                std::cout << "kappa = " << kappa_from_degree_sequence(degrees) << std::endl;

                double avgDeg = double(stubs.size()) / double(N);
                std::string tp;
                if (type == 0) {
                    tp = "RG";
                } else {
                    tp = "PR";
                }

                std::string filename = "../data/LCC/"+ tp +"SF_gamma" + std::to_string(alpha) + "_kmin" + std::to_string(kmin) + ".txt";
                std::ofstream fileLCC(filename);
                filename = "../data/SLCC/"+ tp +"SF_gamma" + std::to_string(alpha) + "_kmin" + std::to_string(kmin) + ".txt";
                std::ofstream fileSLCC(filename);
                filename = "../data/AvgClusterSize/"+ tp +"SF_gamma" + std::to_string(alpha) + "_kmin" + std::to_string(kmin) + ".txt";
                std::ofstream fileAvgClusterSize(filename);

                int mMax = avgDeg*N/2; // Maximum number of edges to be added (fulling the stubs)
                
                omp_set_num_threads(rep); // Set the number of threads for OpenMP
                std::vector<std::vector<int>> lccGlobal; // Largest connected component size
                std::vector<std::vector<int>> secondLargestClusterSizeGlobal; // Second largest connected component size
                std::vector<std::vector<double>> avgClusterSizeGlobal; // Average cluster size

                #pragma omp parallel for
                for (int thread = 0; thread < rep; ++thread){
                    auto localstub = stubs;
                    LinkedGraph graph(N); // Create an empty graph with N nodes
                    std::vector<int> lcc; // Largest connected component size   ];
                    std::vector<int> secondLargestClusterSize; // Second largest connected component size
                    std::vector<double> avgClusterSize; // Average cluster size
                    for(int i = 0; i < mMax; i += steps) {
                        if (type == 0) {
                            graph.addEdgesSF(steps,localstub);
                        } else {
                            graph.addEdgesSFPR(steps,localstub);
                        }
                        lcc.push_back(graph.getLCC());
                        secondLargestClusterSize.push_back(graph.getSecondLargestClusterSize());
                        avgClusterSize.push_back(graph.getAverageClusterSize());
                    }
                    # pragma omp critical
                    {
                       lccGlobal.push_back(lcc);
                       secondLargestClusterSizeGlobal.push_back(secondLargestClusterSize);
                       avgClusterSizeGlobal.push_back(avgClusterSize);
                    }
                }
                // Calculate the average LCC, second largest cluster size, and average cluster size across all threads
                for (int i = 0; i < mMax; i += steps) {
                    double lcc = 0.0;
                    double secondLargestClusterSize = 0.0;
                    double avgClusterSize = 0.0;

                    for (int thread = 0; thread < rep; ++thread) {
                        lcc += lccGlobal[thread][i/steps];
                        secondLargestClusterSize += secondLargestClusterSizeGlobal[thread][i/steps];
                        avgClusterSize += avgClusterSizeGlobal[thread][i/steps];
                    }

                    lcc /= rep;
                    secondLargestClusterSize /= rep;
                    avgClusterSize /= rep;

                    fileLCC << double(i)/double(mMax) << " " << lcc/double(N) << std::endl;
                    fileSLCC << double(i)/double(mMax) << " " << secondLargestClusterSize/double(N) << std::endl;
                    fileAvgClusterSize << double(i)/double(mMax) << " " << avgClusterSize << std::endl;
                }   
            }
        }
    }
}
