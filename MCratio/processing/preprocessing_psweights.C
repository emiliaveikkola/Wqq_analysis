#include <TH3D.h>
#include <TH1D.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <iostream>

void preprocessing_psweights() {
    // Define file names
    std::string inputFileName = "../output/MCSCALEDRun2_scales.root";
    std::string outputFileName = "../processed_histograms/psweights.root";

    // Open the input file
    TFile* inputFile = new TFile(inputFileName.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return;
    }

    // Open the output file
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not open output file " << outputFileName << std::endl;
        inputFile->Close();
        delete inputFile;
        return;
    }

    // List of PS weight scale suffixes
    std::vector<std::string> suffixes = {
        "fsr:murfac=0.707", "fsr:murfac=1.414", "fsr:murfac=0.5",
        "fsr:murfac=2.0", "fsr:murfac=0.25", "fsr:murfac=4.0", "fsr:g2gg:murfac=0.5",
        "fsr:g2gg:murfac=2.0", "fsr:g2qq:murfac=0.5", "fsr:g2qq:murfac=2.0",
        "fsr:q2qg:murfac=0.5", "fsr:q2qg:murfac=2.0", "fsr:x2xg:murfac=0.5",
        "fsr:x2xg:murfac=2.0", "fsr:g2gg:cns=-2.0", "fsr:g2gg:cns=2.0",
        "fsr:g2qq:cns=-2.0", "fsr:g2qq:cns=2.0", "fsr:q2qg:cns=-2.0",
        "fsr:q2qg:cns=2.0", "fsr:x2xg:cns=-2.0", "fsr:x2xg:cns=2.0",
        "isr:murfac=0.707", "isr:murfac=1.414", "isr:murfac=0.5", "isr:murfac=2.0",
        "isr:murfac=0.25", "isr:murfac=4.0", "isr:g2gg:murfac=0.5", "isr:g2gg:murfac=2.0",
        "isr:g2qq:murfac=0.5", "isr:g2qq:murfac=2.0", "isr:q2qg:murfac=0.5",
        "isr:q2qg:murfac=2.0", "isr:x2xg:murfac=0.5", "isr:x2xg:murfac=2.0",
        "isr:g2gg:cns=-2.0", "isr:g2gg:cns=2.0", "isr:g2qq:cns=-2.0", "isr:g2qq:cns=2.0",
        "isr:q2qg:cns=-2.0", "isr:q2qg:cns=2.0", "isr:x2xg:cns=-2.0", "isr:x2xg:cns=2.0"
    };

    // Process both Mass and Pt histograms
    std::vector<std::string> histTypes = {"h3MassFlavorPairs_DATAMC"};

    for (const auto& histType : histTypes) {
        for (const auto& suffix : suffixes) {
            // Construct histogram name
            std::string histName = histType + "_" + suffix;
            TH3D* h3 = dynamic_cast<TH3D*>(inputFile->Get(histName.c_str()));
            if (!h3) {
                std::cerr << "Error: Histogram " << histName << " not found in input file." << std::endl;
                continue;
            }

            // Define projections
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

            // Normalize each projection by its own integral
            for (auto& [name, proj] : projections) {
                double integral = proj->Integral();
                if (integral != 0) {
                    proj->Scale(1.0 / integral);
                } else {
                    std::cerr << "Warning: Integral of " << name << " is zero; skipping normalization." << std::endl;
                }
            }

            // Write projections to output file
            outputFile->cd();
            for (auto& [name, proj] : projections) {
                proj->Write();
            }
        }
    }

    // Clean up
    outputFile->Close();
    inputFile->Close();
    delete outputFile;
    delete inputFile;
}