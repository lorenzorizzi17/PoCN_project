#include <iostream>
#include <fstream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////
//// Displays the Large Connected Component for different rules (ER, PR, SR, BF)
//// when the number of edges is varied
////////////////////////////////////////////////////////////////////////////////////////
void plotLCC(){
    gStyle->SetOptStat(0);
    int N = 1e6;

    TGraph* graphs[4];
    int colors[4] = {kBlue+1, kRed+1, kGreen+2, kOrange+7};

    graphs[0] = new TGraph("../data/LCC/ER.txt");
    graphs[1] = new TGraph("../data/LCC/PR.txt");
    graphs[2] = new TGraph("../data/LCC/SR.txt");
    graphs[3] = new TGraph("../data/LCC/BF.txt");

    for(int i = 0; i < 4; i++) {
        graphs[i]->SetTitle("Order parameter S; m / N; S");
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerSize(0.7);
        graphs[i]->SetMinimum(0);
        graphs[i]->SetMaximum(1);
        graphs[i]->GetXaxis()->SetLimits(0, 1);
    }
    // Divide by N the x vectors
    for (int i = 0; i < 4; i++) {
        int n = graphs[i]->GetN();
        double *x = graphs[i]->GetX();
        for (int j = 0; j < n; j++) {
            x[j] /= double(N);
        }
    }

    graphs[0]->Draw("ALP");
    graphs[1]->Draw("LP SAME");
    graphs[2]->Draw("LP SAME");
    graphs[3]->Draw("LP SAME");
    

    TLegend *legend = new TLegend(0.15, 0.65, 0.4, 0.85);
    legend->SetFillColorAlpha(kWhite, 0.8);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    for (int i = 0; i < 4; i++){
        legend->AddEntry(graphs[i], (i == 0 ? "Erdos-Renyi" : (i == 1 ? "Product Rule" : (i == 2 ? "Sum Rule" : "BF Rule"))), "lp");
    }
    legend->Draw();
}


////////////////////////////////////////////////////////////////////////////////////////
//// Displays the Average Cluster Size for different rules (ER, PR, SR, BF)
//// when the number of edges is varied
////////////////////////////////////////////////////////////////////////////////////////
void plotAvgClusterSize(){
    gStyle->SetOptStat(0);
    int N = 1e6;

    TCanvas* c = new TCanvas();
    c->SetLogy();
    TGraph* graphs[4];
    int colors[4] = {kBlue+1, kRed+1, kGreen+2, kOrange+7};

    graphs[0] = new TGraph("../data/AvgClusterSize/ERaverage.txt");
    graphs[1] = new TGraph("../data/AvgClusterSize/PRaverage.txt");
    graphs[2] = new TGraph("../data/AvgClusterSize/SRaverage.txt");
    graphs[3] = new TGraph("../data/AvgClusterSize/BFaverage.txt");

    for(int i = 0; i < 4; i++) {
        graphs[i]->SetTitle("Average cluster size; m / N; #chi");
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerSize(0.7);
        graphs[i]->SetMinimum(0.99);
        graphs[i]->SetMaximum(10000);
        graphs[i]->GetXaxis()->SetLimits(0.3, 1);
    }
    // Divide by N the x vectors
    for (int i = 0; i < 4; i++) {
        int n = graphs[i]->GetN();
        double *x = graphs[i]->GetX();
        for (int j = 0; j < n; j++) {
            x[j] /= double(N);
        }
    }

    graphs[0]->Draw("ALP");
    graphs[1]->Draw("LP SAME");
    graphs[2]->Draw("LP SAME");
    graphs[3]->Draw("LP SAME");

    TLegend *legend = new TLegend(0.15, 0.65, 0.4, 0.85);
    legend->SetFillColorAlpha(kWhite, 0.8);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    for (int i = 0; i < 4; i++){
        legend->AddEntry(graphs[i], (i == 0 ? "Erdos-Renyi" : (i == 1 ? "Product Rule" : (i == 2 ? "Sum Rule" : "BF Rule"))), "lp");
    }
    legend->Draw();
}


////////////////////////////////////////////////////////////////////////////////////////
//// Displays the finite cluster size distribution for different values of m
//// for a given type of graph (ER, PR, SR, BF)
////////////////////////////////////////////////////////////////////////////////////////

void HistoClusterDistribution(std::string type, std::string m1, std::string m2, std::string m3, std::string m4, std::string m5){
    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas();
    c->SetLogx();
    c->SetLogy();

    const int numGraphs = 5;
    std::string mValues[numGraphs] = {m1, m2, m3, m4, m5};
    std::vector<double> x_vals[numGraphs];
    std::vector<double> y_vals[numGraphs];

    // gathering data
    for (int i = 0; i < numGraphs; ++i) {
        std::ifstream file("../data/ClusterDistribution/"+type+"_" + mValues[i] + ".txt");
        double size, count;
        while (file >> size >> count) {
            if (size > 0 && count > 0) { 
                x_vals[i].push_back(size);
                y_vals[i].push_back(count);
            }
        }
        file.close();
    }

    TGraph* graphs[numGraphs];
    int colors[numGraphs] = {kBlue+1, kRed+1, kGreen+2, kOrange+7, kMagenta+2};

    for (int i = 0; i < numGraphs; ++i) {
        graphs[i] = new TGraph(x_vals[i].size(), x_vals[i].data(), y_vals[i].data());
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetMarkerSize(.8);
        graphs[i]->SetLineColor(colors[i]);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetTitle("Finite cluster size distribution (PR);Cluster size s;Occurrences");
        graphs[i]->GetXaxis()->SetLimits(1, 5*1e2);
    }

    graphs[0]->Draw("AP");
    for (int i = 1; i < numGraphs; ++i) {
        graphs[i]->Draw("LP SAME");
    }

    TLegend *legend = new TLegend(0.15, 0.65, 0.4, 0.85);
    for (int i = 0; i < numGraphs; ++i) {
        legend->AddEntry(graphs[i], ("m/N = " + std::to_string(std::stod(mValues[i]) / 1e6)).c_str(), "lp");
    }
    legend->SetFillColorAlpha(kWhite, 0.8);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->Draw();
}


void plotSF(std::string fam, std::string type, double alpha1, double alpha2, double alpha3, double alpha4, int kmin){
    gStyle->SetOptStat(0); // Niente box statistiche

    TGraph *gr = new TGraph(("../data/"+fam+"/"+type+"SF_gamma"+ std::to_string(alpha1)+"_kmin"+std::to_string(kmin)+".txt").c_str());
    gr->SetTitle("Order parameter S (scale free nets); m / M; S");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue+1);
    gr->SetLineColor(kBlue+1);
    gr->SetLineWidth(2);
    gr->SetMarkerSize(0.7);
    gr->SetMinimum(0);
    gr->SetMaximum(1);
    gr->GetXaxis()->SetLimits(0, 1);
    gr->Draw("ALP");  // A=axis, L=line, P=points

    TGraph *gr2 = new TGraph(("../data/"+fam+"/"+type+"SF_gamma"+ std::to_string(alpha2)+"_kmin"+std::to_string(kmin)+".txt").c_str());
    gr2->SetMarkerStyle(23);
    gr2->SetMarkerColor(kRed+1);
    gr2->SetLineColor(kRed+1);
    gr2->SetLineWidth(2);
    gr2->SetMarkerSize(0.7);
    gr2->SetMinimum(0);
    gr2->SetMaximum(1);
    gr2->GetXaxis()->SetLimits(0, 1);
    gr2->Draw("LP SAME");

    TGraph *gr3 = new TGraph(("../data/"+fam+"/"+type+"SF_gamma"+ std::to_string(alpha3)+"_kmin"+std::to_string(kmin)+".txt").c_str());
    gr3->SetMarkerStyle(23);
    gr3->SetMarkerColor(kOrange+7);
    gr3->SetLineColor(kOrange+7);
    gr3->SetLineWidth(2);
    gr3->SetMarkerSize(0.7);
    gr3->SetMinimum(0);
    gr3->SetMaximum(1);
    gr3->GetXaxis()->SetLimits(0, 1);
    gr3->Draw("LP SAME");

    TGraph *gr4 = new TGraph(("../data/"+fam+"/"+type+"SF_gamma"+ std::to_string(alpha4)+"_kmin"+std::to_string(kmin)+".txt").c_str());
    gr4->SetMarkerStyle(23);
    gr4->SetMarkerColor(kGreen+7);
    gr4->SetLineColor(kGreen+7);
    gr4->SetLineWidth(2);
    gr4->SetMarkerSize(0.7);
    gr4->SetMinimum(0);
    gr4->SetMaximum(1);
    gr4->GetXaxis()->SetLimits(0, 1);
    gr4->Draw("LP SAME");

    //draw a legend
    TLegend *legend = new TLegend(0.15, 0.65, 0.4, 0.85);
    legend->AddEntry(gr, ("gamma = " + std::to_string(alpha1)+ ", kmin =" +std::to_string(kmin)).c_str() , "lp");
    legend->AddEntry(gr2, ("gamma = " + std::to_string(alpha2)+ ", kmin ="+std::to_string(kmin)).c_str(), "lp");
    legend->AddEntry(gr3, ("gamma = " + std::to_string(alpha3)+ ", kmin ="+ std::to_string(kmin)).c_str(), "lp");
    legend->AddEntry(gr4, ("gamma = " + std::to_string(alpha4)+ ", kmin ="+ std::to_string(kmin)).c_str(), "lp");
    legend->SetFillColorAlpha(kWhite, 0.8);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->Draw();


    // Draw two vertical lines at x = 0.388 and x = 0.6
    TLine *line1 = new TLine(1./(2.80898-1), 0, 1./(2.80898-1), 0.5);
    line1->SetLineColor(kOrange+2);
    line1->SetLineStyle(2); // Dashed line
    line1->SetLineWidth(2);
    line1->Draw("same");
}





void HistoDD() {
    std::ifstream file("./data/DegreeDistribution/ER.txt");

    TH1F* hist = new TH1F("DegDist", "Degree Distribution", 21, 0, 21); 

    int value;
    while (file >> value) {
        hist->Fill(value);
    }
    file.close();

    TCanvas* c = new TCanvas("c", "Degree Distribution", 900, 700);

    hist->SetTitle("Log-Log plot of Degree Distribution;Degree;Occurrences");
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.0);
    hist->SetLineColor(kBlue + 2);
    hist->SetLineWidth(2);
    hist->SetStats(false);

    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleOffset(1.2);

    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleOffset(1.4);

    hist->GetXaxis()->SetRangeUser(1, 100);
    c->SetLogx();
    c->SetLogy();

    hist->Draw();
// fit a power law to the histogram
    TF1* fit = new TF1("fit", "[0]*pow(x, [1])", 10, 20);
    fit->SetParameters(4.54541e7, -3.2);  // [0] = A, [1] = alpha (esponente)
    hist->Fit(fit, "R"); // Fit the histogram with the power law function       


    c->SetGridx();
    c->SetGridy();

    c->Update();
}


///////////////////////////////////////////////
// Displays the scaling law for Delta M
// for different rules (ER, PR, SR, BF)
// when the number of edges is varied
///////////////////////////////////////////////

void plotDeltaM(double exponent){
    gStyle->SetOptStat(0); 

    TGraph *graphs[4];
    int colors[4] = {kBlue+1, kRed+1, kOrange+1, kGreen+1};
    graphs[0] = new TGraph("../data/deltaM/ER.txt");
    graphs[1] = new TGraph("../data/deltaM/PR.txt");
    graphs[2] = new TGraph("../data/deltaM/SR.txt");
    graphs[3] = new TGraph("../data/deltaM/BF.txt");

    for (int i = 0; i < 4; i++) {
        graphs[i]->SetTitle("Scaling law for #Delta;N;#Delta / N");
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerSize(0.7);
        graphs[i]->SetMinimum(-0.05);
        graphs[i]->SetMaximum(0.28);
        int n = graphs[i]->GetN();       
        double *x = graphs[i]->GetX();         
        double *y = graphs[i]->GetY();          
        for (int i = 0; i < n; ++i) {
            if (x[i] != 0)               
                y[i] = y[i] / TMath::Power(x[i],exponent);      
            else
                y[i] = 0;                
        }
    }

    graphs[0]->Draw("ALP");  
    graphs[1]->Draw("LP SAME");  
    graphs[2]->Draw("LP SAME");  
    graphs[3]->Draw("LP SAME");  

    TLegend *legend = new TLegend(0.15, 0.65, 0.4, 0.85);
    legend->AddEntry(graphs[0], "ER", "lp");
    legend->AddEntry(graphs[1], "BF", "lp");
    legend->AddEntry(graphs[2], "PR", "lp");
    legend->AddEntry(graphs[3], "SR", "lp");

    legend->SetFillColorAlpha(kWhite, 0.8);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->Draw();

    gPad->Modified();
    gPad->Update();
}


