#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <string>

// Definição da função de massa invariante
double calcular_massa_invariante(const std::vector<float>& pt, const std::vector<float>& eta, const std::vector<float>& phi) {
    if (pt.size() >= 2) {
        return sqrt(2 * pt[0] * pt[1] * (TMath::CosH(eta[0] - eta[1]) - TMath::Cos(phi[0] - phi[1])));
    }
    return -1.0; 
}

void Zmasscutfit() {
    // Diretórios para análise
    std::vector<std::string> diretorios = {
        "/opendata/eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/2430000",
        "/opendata/eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/2520000",
        "/opendata/eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/80000"
    };

    std::vector<double> e_massas_invariantes;

    for (const auto& diretorio : diretorios) {
        for (const auto& entry : std::filesystem::directory_iterator(diretorio)) {
            std::string file_path = entry.path();
            TFile file(file_path.c_str(), "READ");
            if (!file.IsOpen()) continue;

            TTreeReader reader("Events", &file);
            TTreeReaderArray<float> Electron_pt(reader, "Electron_pt");
            TTreeReaderArray<float> Electron_eta(reader, "Electron_eta");
            TTreeReaderArray<float> Electron_phi(reader, "Electron_phi");
            TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
            TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
            TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
            TTreeReaderArray<float> Tau_pt(reader, "Tau_pt");
            TTreeReaderArray<float> Tau_eta(reader, "Tau_eta");
            TTreeReaderArray<float> Tau_phi(reader, "Tau_phi");

            while (reader.Next()) {
                std::vector<std::pair<float, int>> leptons;
                for (size_t i = 0; i < Electron_pt.GetSize(); ++i) {
                    if (Electron_pt[i] > 20 && fabs(Electron_eta[i]) < 2.5) {
                        leptons.emplace_back(Electron_pt[i], i); 
                    }
                }
                for (size_t i = 0; i < Muon_pt.GetSize(); ++i) {
                    if (Muon_pt[i] > 20 && fabs(Muon_eta[i]) < 2.5) {
                        leptons.emplace_back(Muon_pt[i], i + Electron_pt.GetSize()); 
                    }
                }
                for (size_t i = 0; i < Tau_pt.GetSize(); ++i) {
                    if (Tau_pt[i] > 20 && fabs(Tau_eta[i]) < 2.5) {
                        leptons.emplace_back(Tau_pt[i], i + Electron_pt.GetSize() + Muon_pt.GetSize()); 
                    }
                }

                std::sort(leptons.rbegin(), leptons.rend());

                if (leptons.size() >= 2) {
                    int idx1 = leptons[0].second;
                    int idx2 = leptons[1].second;

                    float pt1, eta1, phi1, pt2, eta2, phi2;

                    if (idx1 < Electron_pt.GetSize()) {
                        pt1 = Electron_pt[idx1];
                        eta1 = Electron_eta[idx1];
                        phi1 = Electron_phi[idx1];
                    } else if (idx1 < Electron_pt.GetSize() + Muon_pt.GetSize()) {
                        pt1 = Muon_pt[idx1 - Electron_pt.GetSize()];
                        eta1 = Muon_eta[idx1 - Electron_pt.GetSize()];
                        phi1 = Muon_phi[idx1 - Electron_pt.GetSize()];
                    } else {
                        pt1 = Tau_pt[idx1 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                        eta1 = Tau_eta[idx1 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                        phi1 = Tau_phi[idx1 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                    }

                    if (idx2 < Electron_pt.GetSize()) {
                        pt2 = Electron_pt[idx2];
                        eta2 = Electron_eta[idx2];
                        phi2 = Electron_phi[idx2];
                    } else if (idx2 < Electron_pt.GetSize() + Muon_pt.GetSize()) {
                        pt2 = Muon_pt[idx2 - Electron_pt.GetSize()];
                        eta2 = Muon_eta[idx2 - Electron_pt.GetSize()];
                        phi2 = Muon_phi[idx2 - Electron_pt.GetSize()];
                    } else {
                        pt2 = Tau_pt[idx2 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                        eta2 = Tau_eta[idx2 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                        phi2 = Tau_phi[idx2 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                    }

                    std::vector<float> pt_values = {pt1, pt2};
                    std::vector<float> eta_values = {eta1, eta2};
                    std::vector<float> phi_values = {phi1, phi2};
                    double massa_invariante = calcular_massa_invariante(pt_values, eta_values, phi_values);

                    if (massa_invariante >= 0) {
                        e_massas_invariantes.push_back(massa_invariante);
                    }
                }
            }
        }
    }

    // 1. Criamos o histograma com os dados de massa invariante
    TH1F* hMassaInvariante = new TH1F("hMassaInvariante", "Distribuição de Massa Invariante", 150, 60, 120);
    for (const auto& massa : e_massas_invariantes) {
        if (massa >= 0) hMassaInvariante->Fill(massa);
    }

    // 2. Definimos a variável para o RooFit
    RooRealVar x("x", "Massa Invariante (GeV/c^{2})", 60, 120);
    
    // 3. Criamos o RooDataHist a partir do histograma
    RooDataHist data("data", "Dados de Massa Invariante", RooArgList(x), hMassaInvariante);

    // 4. Definimos as funções para o sinal (Crystal Ball) e fundo (Exponencial)
    RooRealVar mean("mean", "Média", 91.2, 90, 93); // Média do Z boson
    RooRealVar sigma("sigma", "Desvio padrão", 2.3, 1.5, 3.0); // Largura do pico
    RooRealVar alpha("alpha", "Parâmetro alpha", 1.3, 0.8, 2.5); // Parâmetro alpha (cauda)
    RooRealVar n("n", "Parâmetro n", 7, 3, 15); // Parâmetro n (cauda)
    RooCBShape signal("signal", "Função Crystal Ball", x, mean, sigma, alpha, n);

    // Função de fundo exponencial
    RooRealVar tau("tau", "Parâmetro Tau (Fundo Exponencial)", -0.2, -1.0, 0.0);
    RooExponential background("background", "Fundo Exponencial", x, tau);

    // 5. Criamos o modelo total como uma soma (signal + background)
    RooRealVar fracSignal("fracSignal", "Fração de Sinal", 0.6, 0.2, 0.4); // Fração do sinal
    RooAddPdf model("model", "Sinal + Fundo", RooArgList(signal, background), RooArgList(fracSignal));

    // 6. Ajuste do modelo aos dados
    model.fitTo(data, RooFit::PrintLevel(-1), RooFit::Range(60, 120));

    // 7. Criamos o gráfico para o ajuste
    RooPlot* frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame, RooFit::LineColor(kBlue)); // Ajuste total
    model.plotOn(frame, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed)); // Sinal (Crystal Ball)
    model.plotOn(frame, RooFit::Components("background"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed)); // Fundo (Exponencial)

    // 8. Criação do canvas e exibição do gráfico
    TCanvas* canvas = new TCanvas("canvasFit", "Ajuste do Sinal + Fundo", 800, 600);
    frame->SetTitle("Ajuste de Massa Invariante com Sinal + Fundo");
    frame->GetXaxis()->SetTitle("Massa Invariante (GeV/c^{2})");
    frame->GetYaxis()->SetTitle("Eventos");
    frame->Draw();

    // 9. Adiciona a legenda
    TLegend* legend = new TLegend(0.7, 0.6, 0.9, 0.9);
    legend->AddEntry(frame->getObject(0), "Dados", "P");
    legend->AddEntry(frame->getObject(1), "Ajuste Total", "L");
    legend->AddEntry(frame->getObject(2), "Sinal (Crystal Ball)", "L");
    legend->AddEntry(frame->getObject(3), "Fundo (Exponencial)", "L");
    legend->Draw();

    // 10. Salva a imagem
    canvas->SaveAs("mass_Z_fit.png");
    
    // Libera a memória
    delete canvas;
    delete frame;
    delete legend;
    delete hMassaInvariante;
}
