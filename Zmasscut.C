#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <string>    

// Função para calcular a massa invariante com dois leptons
double calcular_massa_invariante(const std::vector<float>& pt, const std::vector<float>& eta, const std::vector<float>& phi) {
    if (pt.size() >= 2) {
        return sqrt(2 * pt[0] * pt[1] * (TMath::CosH(eta[0] - eta[1]) - TMath::Cos(phi[0] - phi[1])));
    }
    return -1.0;  // Retorna -1 caso não tenha pelo menos dois leptons
}

void Zmasscut() {
    // Diretórios de entrada
    std::vector<std::string> diretorios = {
        "/opendata/eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/2430000",
        "/opendata/eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/2520000"
    };

    // Vetor para armazenar as massas invariantes
    std::vector<double> e_massas_invariantes;

    // Inicializando o histograma de massa invariante
    TH1F* hMassaInvariante = new TH1F("hMassaInvariante", "Invariant Mass Distribution", 200, 70, 110);

    // Processar os arquivos
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
            TTreeReaderArray<float> Jet_pt(reader, "Jet_pt");
            TTreeReaderArray<float> Jet_eta(reader, "Jet_eta");
            TTreeReaderArray<float> Jet_phi(reader, "Jet_phi");
            TTreeReaderArray<float> Tau_pt(reader, "Tau_pt");
            TTreeReaderArray<float> Tau_eta(reader, "Tau_eta");
            TTreeReaderArray<float> Tau_phi(reader, "Tau_phi");

            while (reader.Next()) {
                // Vetor de leptons com pT e eta
                std::vector<std::pair<float, int>> leptons; // (pT, índice)
                for (size_t i = 0; i < Electron_pt.GetSize(); ++i) {
                    if (Electron_pt[i] > 20 && fabs(Electron_eta[i]) < 2.5) {
                        leptons.emplace_back(Electron_pt[i], i); // Elétrons
                    }
                }
                for (size_t i = 0; i < Muon_pt.GetSize(); ++i) {
                    if (Muon_pt[i] > 20 && fabs(Muon_eta[i]) < 2.5) {
                        leptons.emplace_back(Muon_pt[i], i + Electron_pt.GetSize()); // Múons
                    }
                }
                for (size_t i = 0; i < Tau_pt.GetSize(); ++i) {
                    if (Tau_pt[i] > 20 && fabs(Tau_eta[i]) < 2.5) {
                        leptons.emplace_back(Tau_pt[i], i + Electron_pt.GetSize() + Muon_pt.GetSize()); // Taus
                    }
                }

                // Ordenar leptons por pT em ordem decrescente
                std::sort(leptons.rbegin(), leptons.rend());

                // Selecionar os dois leptons com maior pT
                if (leptons.size() >= 2) {
                    int idx1 = leptons[0].second;
                    int idx2 = leptons[1].second;

                    float pt1, eta1, phi1, pt2, eta2, phi2;

                    // Atribuir valores de pt, eta, phi para cada lepton
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

                    // Calcular a massa invariante
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

    // Plot da Massa Invariante
    TCanvas* canvas = new TCanvas("canvasInvariantMass", "Invariant Mass Distribution", 800, 600);
    hMassaInvariante->SetLineColor(kBlack);
    canvas->SetLogy();
    for (const auto& massa : e_massas_invariantes) {
        if (massa >= 0) hMassaInvariante->Fill(massa);
    }
    hMassaInvariante->Draw();
    hMassaInvariante->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    hMassaInvariante->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("invariant_mass_Z_v2.png");

    // Limpeza
    delete canvas;
    delete hMassaInvariante;
}
