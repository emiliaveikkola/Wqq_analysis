#include <TH3D.h>
#include <TH1D.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <iostream>

// Function to process and save histograms from different files (scaled MC, unscaled MC, and data)
void preprocessing_pt() {
    
    // Define file names
    std::string scaledMCFileName = "../output/output_MCSCALEDRun2_tagprobe_test5.root";
    std::string unscaledMCFileName = "../output/output_MCRun2_tagprobe_test4.root";
    std::string dataFileName = "../output/output_DATARun2_tagprobe_test4.root";
    std::string outputFileName = "../processed_histograms/processed_histograms_pt.root";
    
    // Open the scaled MC file
    TFile* scaledMCFile = new TFile(scaledMCFileName.c_str(), "READ");
    if (!scaledMCFile || scaledMCFile->IsZombie()) {
        std::cerr << "Error: Could not open scaled MC file " << scaledMCFileName << std::endl;
        return;
    }

    // Open the unscaled MC file
    TFile* unscaledMCFile = new TFile(unscaledMCFileName.c_str(), "READ");
    if (!unscaledMCFile || unscaledMCFile->IsZombie()) {
        std::cerr << "Error: Could not open unscaled MC file " << unscaledMCFileName << std::endl;
        scaledMCFile->Close();
        delete scaledMCFile;
        return;
    }

    // Open the data file
    TFile* dataFile = new TFile(dataFileName.c_str(), "READ");
    if (!dataFile || dataFile->IsZombie()) {
        std::cerr << "Error: Could not open data file " << dataFileName << std::endl;
        scaledMCFile->Close();
        unscaledMCFile->Close();
        delete scaledMCFile;
        delete unscaledMCFile;
        return;
    }

    // Open the output file to save processed histograms
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not open output file " << outputFileName << std::endl;
        scaledMCFile->Close();
        unscaledMCFile->Close();
        dataFile->Close();
        delete scaledMCFile;
        delete unscaledMCFile;
        delete dataFile;
        return;
    }

    // List of histogram names and unique suffixes for scaled MC
    std::vector<std::pair<std::string, std::string>> histograms = {
        //Pt
        {"h3PtFlavorPairs_DATAMC_RESP", "RESP"},
        {"h3PtFlavorPairs_DATAMC_FSR_0.996", "FSR0996"},
        {"h3PtFlavorPairs_DATAMC_FSR_0.995", "FSR0995"}, // -0.5%
        {"h3PtFlavorPairs_DATAMC_FSR_1.005", "FSR1005"}, // 0.5%
        {"h3PtFlavorPairs_DATAMC_FSR_1.010", "FSR1010"}, // 1%
        {"h3PtFlavorPairs_DATAMC_FSR_1.020", "FSR1020"}, // 2%
        {"h3PtFlavorPairs_DATAMC_FSR_0.990", "FSR0990"}, // -1%
        {"h3PtFlavorPairs_DATAMC_FSR_0.980", "FSR0980"}, // -2%
        {"h3PtFlavorPairs_DATAMC_JER80", "JER080"},
        {"h3PtFlavorPairs_DATAMC_JER90", "JER090"},
        {"h3PtFlavorPairs_DATAMC_JER105", "JER105"},
        {"h3PtFlavorPairs_DATAMC_JER110", "JER110"},
        {"h3PtFlavorPairs_DATAMC_JER120", "JER120"},
        {"h3PtFlavorPairs_DATAMC_JER_pt80", "JER080pt"},
        {"h3PtFlavorPairs_DATAMC_JER_pt90", "JER090pt"},
        {"h3PtFlavorPairs_DATAMC_JER_pt105", "JER105pt"},
        {"h3PtFlavorPairs_DATAMC_JER_pt110", "JER110pt"},
        {"h3PtFlavorPairs_DATAMC_JER_pt120", "JER120pt"},
        {"h3PtFlavorPairs_DATAMC_ISR40", "ISR040"},
        {"h3PtFlavorPairs_DATAMC_ISR62", "ISR062"},
        {"h3PtFlavorPairs_DATAMC_ISR70", "ISR070"},
        {"h3PtFlavorPairs_DATAMC_ISR77", "ISR077"},
        {"h3PtFlavorPairs_DATAMC_ISR120", "ISR120"},
        {"h3PtFlavorPairs_DATAMC_ISR130", "ISR130"},
        {"h3PtFlavorPairs_DATAMC_ISR160", "ISR160"}
    };

    // Process scaled MC histograms
    for (const auto& [histName, suffix] : histograms) {
        TH3D* h3 = dynamic_cast<TH3D*>(scaledMCFile->Get(histName.c_str()));
        if (!h3) {
            std::cerr << "Error: Histogram " << histName << " not found in scaled MC file." << std::endl;
            continue;
        }

        // Projections for all required categories with specified ranges
        std::vector<std::pair<std::string, TH1D*>> projections = {
            {"h3all_" + suffix, h3->ProjectionZ(("h3all_" + suffix).c_str(), 1, 3, 1, 3)},

            {"h3gencsall_" + suffix, h3->ProjectionZ(("h3gencsall_" + suffix).c_str(), 1, 1, 1, 3)},
            {"h3genudall_" + suffix, h3->ProjectionZ(("h3genudall_" + suffix).c_str(), 2, 2, 1, 3)},
            {"h3genxall_" + suffix, h3->ProjectionZ(("h3genxall_" + suffix).c_str(), 3, 3, 1, 3)},

            {"h3tagcsall_" + suffix, h3->ProjectionZ(("h3tagcsall_" + suffix).c_str(), 1, 3, 1, 1)},
            {"h3tagudall_" + suffix, h3->ProjectionZ(("h3tagudall_" + suffix).c_str(), 1, 3, 2, 2)},
            {"h3tagxall_" + suffix, h3->ProjectionZ(("h3tagxall_" + suffix).c_str(), 1, 3, 3, 3)},

            {"h3gencstagcs_" + suffix, h3->ProjectionZ(("h3gencstagcs_" + suffix).c_str(), 1, 1, 1, 1)},
            {"h3gencstagud_" + suffix, h3->ProjectionZ(("h3gencstagud_" + suffix).c_str(), 1, 1, 2, 2)},
            {"h3gencstagx_" + suffix, h3->ProjectionZ(("h3gencstagx_" + suffix).c_str(), 1, 1, 3, 3)},

            {"h3genudtagcs_" + suffix, h3->ProjectionZ(("h3genudtagcs_" + suffix).c_str(), 2, 2, 1, 1)},
            {"h3genudtagud_" + suffix, h3->ProjectionZ(("h3genudtagud_" + suffix).c_str(), 2, 2, 2, 2)},
            {"h3genudtagx_" + suffix, h3->ProjectionZ(("h3genudtagx_" + suffix).c_str(), 2, 2, 3, 3)},

            {"h3genxtagcs_" + suffix, h3->ProjectionZ(("h3genxtagcs_" + suffix).c_str(), 3, 3, 1, 1)},
            {"h3genxtagud_" + suffix, h3->ProjectionZ(("h3genxtagud_" + suffix).c_str(), 3, 3, 2, 2)},
            {"h3genxtagx_" + suffix, h3->ProjectionZ(("h3genxtagx_" + suffix).c_str(), 3, 3, 3, 3)}
        };

        // Normalize each projection by the integral of "h3all" from scaled MC
        double allIntegral = projections[0].second->Integral();
        if (allIntegral != 0) {
            for (auto& [name, proj] : projections) {
                proj->Scale(1. / allIntegral);
            }
        }

        // Write each projection to the output file
        outputFile->cd();
        for (auto& [name, proj] : projections) {
            proj->Write();
        }
    }

    // Process unscaled MC histograms (normalized by h3all_unscaled)
    TH3D* h3PtFlavorPairs_DATAMC_MC = (TH3D*)unscaledMCFile->Get("h3PtFlavorPairs_DATAMC");
    if (h3PtFlavorPairs_DATAMC_MC) {
        TH1D* h3all_unscaled = h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3all_unscaled", 1, 3, 1, 3);
        double unscaledAllIntegral = h3all_unscaled->Integral();

        std::vector<std::pair<std::string, TH1D*>> unscaledProjections = {
            {"h3all_unscaled", h3all_unscaled},
            {"h3tagcsall_unscaled", h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3tagcsall_unscaled", 1, 3, 1, 1)},
            {"h3tagudall_unscaled", h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3tagudall_unscaled", 1, 3, 2, 2)},
            {"h3tagxall_unscaled", h3PtFlavorPairs_DATAMC_MC->ProjectionZ("h3tagxall_unscaled", 1, 3, 3, 3)}
        };

        // Normalize each unscaled projection by the integral of h3all_unscaled
        if (unscaledAllIntegral != 0) {
            for (auto& [name, proj] : unscaledProjections) {
                proj->Scale(1.0 / unscaledAllIntegral);
            }
        }

        // Write unscaled projections to the output file
        outputFile->cd();
        for (auto& [name, proj] : unscaledProjections) {
            proj->Write();
        }
    }

    // Process data histograms (normalized by h3all_data)
    TH3D* h3PtFlavorPairs_DATAMC_DATA = (TH3D*)dataFile->Get("h3PtFlavorPairs_DATAMC");
    if (h3PtFlavorPairs_DATAMC_DATA) {
        TH1D* h3all_data = h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);
        double dataAllIntegral = h3all_data->Integral();

        std::vector<std::pair<std::string, TH1D*>> dataProjections = {
            {"h3all_data", h3all_data},
            {"h3tagcsall_data", h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1)},
            {"h3tagudall_data", h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2)},
            {"h3tagxall_data", h3PtFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3)}
        };

        // Normalize each data projection by the integral of h3all_data
        if (dataAllIntegral != 0) {
            for (auto& [name, proj] : dataProjections) {
                proj->Scale(1.0 / dataAllIntegral);
            }
        }

        // Write data projections to the output file
        outputFile->cd();
        for (auto& [name, proj] : dataProjections) {
            proj->Write();
        }
    }

    // Clean up
    outputFile->Close();
    scaledMCFile->Close();
    unscaledMCFile->Close();
    dataFile->Close();
    delete outputFile;
    delete scaledMCFile;
    delete unscaledMCFile;
    delete dataFile;
}