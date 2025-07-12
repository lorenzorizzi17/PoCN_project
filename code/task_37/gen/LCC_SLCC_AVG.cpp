#include <fstream>
#include "./include/graph.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// This main functions is used to compute the LCC, the second largest cluster size and the average cluster size
/// for different graph topologies:
/// 1. Erdos-Renyi
/// 2. Product Rule
/// 3. Sum Rule
/// 4. BF Rule
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){
    int N = 1e6;      // Number of nodes in the graph
    int rep = 5;       // Number of independent realizations for each graph

    int step = 150; 
    int tmax = 1*N;
    LinkedGraph graph(N);
    
    std::ofstream ERfile("../data/LCC/ER.txt");
    std::ofstream ProductRuleFile("../data/LCC/PR.txt");
    std::ofstream SumRuleFile("../data/LCC/SR.txt");
    std::ofstream BFfile("../data/LCC/BF.txt");

    std::ofstream ERsecondFile("../data/SLCC/ERsecond.txt");
    std::ofstream ProductRulesecondFile("../data/SLCC/PRsecond.txt");
    std::ofstream SumRulesecondFile("../data/SLCC/SRsecond.txt");
    std::ofstream BFsecondFile("../data/SLCC/BFsecond.txt");

    std::ofstream ERaverageFile("../data/AvgClusterSize/ERaverage.txt");
    std::ofstream ProductRuleAverageFile("../data/AvgClusterSize/PRaverage.txt");
    std::ofstream SumRuleAverageFile("../data/AvgClusterSize/SRaverage.txt");
    std::ofstream BFaverageFile("../data/AvgClusterSize/BFaverage.txt");
    

    double ERvector [tmax/step] = {0}; 
    double ERsecond [tmax/step] = {0};
    double ERaverage [tmax/step] = {0}; 
    for (int r = 0; r < rep; r++) {
        LinkedGraph graph(N); // Create a graph with 10 nodes
        for(int i = 0; i < tmax; i+= step) {
            graph.addRandomEdges(step); // Add a random edge
            ERvector[int(i/step)] += ((graph.getLCC()/double(N))) / double(rep); // Get the LCC size
            ERsecond[int(i/step)] += ((graph.getSecondLargestClusterSize()/double(N))) / double(rep); // Get the second cluster size
            ERaverage[int(i/step)] += ((graph.getAverageClusterSize())) / double(rep); // Get the average cluster size
        }
        std::cout << "Repetition " << r+1 << " done." << std::endl;
    }
    std::cout << "ER done." << std::endl;
   

    // start again with a new graph
    double PRvector [tmax/step] = {0}; 
    double PRsecond[tmax/step] = {0};
    double PRaverage[tmax/step] = {0};
    for (int r = 0; r < rep; r++) {
        graph = LinkedGraph(N);
        for(int i = 0; i < tmax; i+= step) {
            graph.addRandomEdgesProductRule(step); // Add a random edge
            PRvector[int(i/step)] += ((graph.getLCC()/double(N))) / double(rep); // Get the LCC size
            PRsecond[int(i/step)] += ((graph.getSecondLargestClusterSize()/double(N))) / double(rep); // Get the second cluster size
            PRaverage[int(i/step)] += ((graph.getAverageClusterSize())) / double(rep); // Get the average cluster size
        }
        std::cout << "Repetition " << r+1 << " done." << std::endl;
    }
    std::cout << "Product Rule done." << std::endl;

    // start again with a new graph
    double SRvector [tmax/step] = {0}; 
    double SRsecond[tmax/step] = {0};
    double SRaverage[tmax/step] = {0};
    for (int r = 0; r < rep; r++) {
        graph = LinkedGraph(N);
        for(int i = 0; i < tmax; i+= step) {
            graph.addRandomEdgesSumRule(step); // Add a random edge
            SRvector[int(i/step)] += ((graph.getLCC()/double(N))) / double(rep); // Get the LCC size
            SRsecond[int(i/step)] += ((graph.getSecondLargestClusterSize()/double(N))) / double(rep); // Get the second cluster size
            SRaverage[int(i/step)] += ((graph.getAverageClusterSize())) / double(rep); // Get the average cluster size
        }
        std::cout << "Repetition " << r+1 << " done." << std::endl;
    }
    std::cout << "Sum Rule done." << std::endl;

    // start again with a new graph
    double BFvector [tmax/step] = {0}; 
    double BFsecond[tmax/step] = {0};
    double BFaverage[tmax/step] = {0};
    for (int r = 0; r < rep; r++) {
        graph = LinkedGraph(N);
        for(int i = 0; i < tmax; i+= step) {
            graph.addEdgesBFRule(step); // Add a random edge
            BFvector[int(i/step)] += ((graph.getLCC()/double(N))) / double(rep); // Get the LCC size
            BFsecond[int(i/step)] += ((graph.getSecondLargestClusterSize()/double(N))) / double(rep); // Get the second cluster size
            BFaverage[int(i/step)] += ((graph.getAverageClusterSize())) / double(rep); // Get the average cluster size
        }
        std::cout << "Repetition " << r+1 << " done." << std::endl;
    }
    std::cout << "BF Rule done." << std::endl;

    // Write the results to the files
    for (int i = 0; i < tmax/step; i++) {
        ERfile << step*i << " " << ERvector[i] << std::endl;
        ProductRuleFile << i*step << " " << PRvector[i] << std::endl;
        SumRuleFile << step*i << " " << SRvector[i] << std::endl;
        BFfile << i*step << " " << BFvector[i] << std::endl;
        // Write the second cluster sizes
        ERsecondFile << step*i << " " << ERsecond[int(i)] << std::endl;
        ProductRulesecondFile << step*i << " " << PRsecond[int(i)] << std::endl;
        SumRulesecondFile << step*i << " " << SRsecond[i] << std::endl;
        BFsecondFile << i*step << " " << BFsecond[(i)] << std::endl;
        // Write the average cluster sizes
        ERaverageFile << step*i << " " << ERaverage[i] << std::endl;
        ProductRuleAverageFile << step*i << " " << PRaverage[i] << std::endl;
        SumRuleAverageFile << step*i << " " << SRaverage[i] << std::endl;
        BFaverageFile << i*step << " " << BFaverage[i] << std::endl;
    }
}