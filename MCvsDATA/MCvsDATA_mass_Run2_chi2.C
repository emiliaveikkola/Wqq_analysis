#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include <algorithm>
#include "../minitools/tdrstyle_mod22.C"

// Initialize vectors to store chi-squared values for each category
std::vector<double> chi2_cs_values;
std::vector<double> chi2_ud_values;
std::vector<double> chi2_x_values;
std::vector<double> chi2_all_values;

void MCvsDATA_mass_Run2_chi2() {
    TFile *file = new TFile("../output/output_MCRun2.root", "READ");
    // Open your data and MC files
    TH3D* h3MassFlavorPairs_DATAMC_MC = (TH3D*)file->Get("h3MassFlavorPairs_DATAMC");
    // Project Z for different indices and draw
    TH1D* h3all = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3all", 1, 3, 1, 3);

    TH1D* h3gencsall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencsall", 1, 1, 1, 3);
    TH1D* h3genudall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudall", 2, 2, 1, 3);
    TH1D* h3genxall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxall", 3, 3, 1, 3);

    TH1D* h3tagcsall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3tagcsall", 1, 3, 1, 1);
    TH1D* h3tagudall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3tagudall", 1, 3, 2, 2);
    TH1D* h3tagxall = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3tagxall", 1, 3, 3, 3);

    TH1D* h3gencstagcs = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagcs", 1, 1, 1, 1);
    TH1D* h3gencstagud = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagud", 1, 1, 2, 2);
    TH1D* h3gencstagx = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3gencstagx", 1, 1, 3, 3);
    
    TH1D* h3genudtagcs = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagcs", 2, 2, 1, 1);
    TH1D* h3genudtagud = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagud", 2, 2, 2, 2);
    TH1D* h3genudtagx = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genudtagx", 2, 2, 3, 3);
    
    TH1D* h3genxtagcs = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagcs", 3, 3, 1, 1);
    TH1D* h3genxtagud = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagud", 3, 3, 2, 2);
    TH1D* h3genxtagx = h3MassFlavorPairs_DATAMC_MC->ProjectionZ("h3genxtagx", 3, 3, 3, 3);

    TFile *file2 = new TFile("../output/output_DATARun2.root", "READ");
    TH3D* h3MassFlavorPairs_DATAMC_DATA = (TH3D*)file2->Get("h3MassFlavorPairs_DATAMC");

    TH1D* h3all_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3all_data", 1, 3, 1, 3);

    TH1D* h3tagcsall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagcsall_data", 1, 3, 1, 1);
    TH1D* h3tagudall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagudall_data", 1, 3, 2, 2);
    TH1D* h3tagxall_data = h3MassFlavorPairs_DATAMC_DATA->ProjectionZ("h3tagxall_data", 1, 3, 3, 3);


    h3tagcsall_data->Scale(1./h3all_data->Integral());
    h3tagudall_data->Scale(1./h3all_data->Integral());
    h3tagxall_data->Scale(1./h3all_data->Integral());

    h3all_data->Scale(1./h3all_data->Integral());

    // Define different sets of scale factors (replace this with your actual sets)
std::vector<std::vector<double>> scale_factors = {
    {0.86433, 1.16679, 0.2, 0.2, 1.02945, 0.963425, 0.776569, 1.06566, 0.858425},
    {0.86433, 1.14027, 0.4, 0.2, 1.04424, 0.918284, 0.776569, 1.06873, 0.846995},
    {0.86433, 1.11375, 0.6, 0.2, 1.05903, 0.873137, 0.776569, 1.07179, 0.835574},
    {0.86433, 1.08723, 0.8, 0.2, 1.07382, 0.827995, 0.776569, 1.07486, 0.824145},
    {0.86433, 1.06071, 1.0, 0.2, 1.0886, 0.782854, 0.776569, 1.07792, 0.812715},
    {0.86433, 1.03418, 1.2, 0.2, 1.1034, 0.737707, 0.776569, 1.08098, 0.801295},
    {0.86433, 1.00766, 1.4, 0.2, 1.11818, 0.692565, 0.776569, 1.08405, 0.789866},
    {0.86433, 0.981139, 1.6, 0.2, 1.13297, 0.647425, 0.776569, 1.08711, 0.778436},
    {0.86433, 0.954617, 1.8, 0.2, 1.14776, 0.602278, 0.776569, 1.09018, 0.767014},
    {0.845747, 1.17511, 0.2, 0.4, 1.02508, 0.963425, 0.786861, 1.06439, 0.858425},
    {0.845747, 1.14859, 0.4, 0.4, 1.03987, 0.918284, 0.786861, 1.06745, 0.846995},
    {0.845747, 1.12207, 0.6, 0.4, 1.05466, 0.873137, 0.786861, 1.07052, 0.835574},
    {0.845747, 1.09554, 0.8, 0.4, 1.06945, 0.827995, 0.786861, 1.07358, 0.824145},
    {0.845747, 1.06902, 1.0, 0.4, 1.08424, 0.782854, 0.786861, 1.07665, 0.812715},
    {0.845747, 1.0425, 1.2, 0.4, 1.09903, 0.737707, 0.786861, 1.07971, 0.801295},
    {0.845747, 1.01598, 1.4, 0.4, 1.11382, 0.692565, 0.786861, 1.08277, 0.789866},
    {0.845747, 0.989455, 1.6, 0.4, 1.12861, 0.647425, 0.786861, 1.08584, 0.778436},
    {0.845747, 0.962932, 1.8, 0.4, 1.1434, 0.602278, 0.786861, 1.0889, 0.767014},
    {0.827164, 1.18342, 0.2, 0.6, 1.02072, 0.963425, 0.797154, 1.06311, 0.858425},
    {0.827164, 1.1569, 0.4, 0.6, 1.0355, 0.918284, 0.797154, 1.06618, 0.846995},
    {0.827164, 1.13038, 0.6, 0.6, 1.05029, 0.873137, 0.797154, 1.06924, 0.835574},
    {0.827164, 1.10386, 0.8, 0.6, 1.06508, 0.827995, 0.797154, 1.0723, 0.824145},
    {0.827164, 1.07734, 1.0, 0.6, 1.07987, 0.782854, 0.797154, 1.07537, 0.812715},
    {0.827164, 1.05081, 1.2, 0.6, 1.09466, 0.737707, 0.797154, 1.07843, 0.801295},
    {0.827164, 1.02429, 1.4, 0.6, 1.10945, 0.692565, 0.797154, 1.0815, 0.789866},
    {0.827164, 0.99777, 1.6, 0.6, 1.12424, 0.647425, 0.797154, 1.08456, 0.778436},
    {0.827164, 0.971248, 1.8, 0.6, 1.13903, 0.602278, 0.797154, 1.08763, 0.767014},
    {0.80858, 1.19174, 0.2, 0.8, 1.01635, 0.963425, 0.807447, 1.06184, 0.858425},
    {0.80858, 1.16522, 0.4, 0.8, 1.03114, 0.918284, 0.807447, 1.0649, 0.846995},
    {0.80858, 1.1387, 0.6, 0.8, 1.04593, 0.873137, 0.807447, 1.06796, 0.835574},
    {0.80858, 1.11217, 0.8, 0.8, 1.06072, 0.827995, 0.807447, 1.07103, 0.824145},
    {0.80858, 1.08565, 1.0, 0.8, 1.07551, 0.782854, 0.807447, 1.07409, 0.812715},
    {0.80858, 1.05913, 1.2, 0.8, 1.0903, 0.737707, 0.807447, 1.07716, 0.801295},
    {0.80858, 1.03261, 1.4, 0.8, 1.10508, 0.692565, 0.807447, 1.08022, 0.789866},
    {0.80858, 1.00608, 1.6, 0.8, 1.11987, 0.647425, 0.807447, 1.08329, 0.778436},
    {0.80858, 0.979563, 1.8, 0.8, 1.13466, 0.602278, 0.807447, 1.08635, 0.767014},
    {0.789997, 1.20005, 0.2, 1.0, 1.01198, 0.963425, 0.81774, 1.06056, 0.858425},
    {0.789997, 1.17353, 0.4, 1.0, 1.02677, 0.918284, 0.81774, 1.06363, 0.846995},
    {0.789997, 1.14701, 0.6, 1.0, 1.04156, 0.873137, 0.81774, 1.06669, 0.835574},
    {0.789997, 1.12049, 0.8, 1.0, 1.05635, 0.827995, 0.81774, 1.06975, 0.824145},
    {0.789997, 1.09397, 1.0, 1.0, 1.07114, 0.782854, 0.81774, 1.07282, 0.812715},
    {0.789997, 1.06744, 1.2, 1.0, 1.08593, 0.737707, 0.81774, 1.07588, 0.801295},
    {0.789997, 1.04092, 1.4, 1.0, 1.10072, 0.692565, 0.81774, 1.07895, 0.789866},
    {0.789997, 1.0144, 1.6, 1.0, 1.11551, 0.647425, 0.81774, 1.08201, 0.778436},
    {0.789997, 0.987878, 1.8, 1.0, 1.1303, 0.602278, 0.81774, 1.08507, 0.767014},
    {0.771414, 1.20837, 0.2, 1.2, 1.00762, 0.963425, 0.828032, 1.05928, 0.858425},
    {0.771414, 1.18185, 0.4, 1.2, 1.0224, 0.918284, 0.828032, 1.06235, 0.846995},
    {0.771414, 1.15533, 0.6, 1.2, 1.03719, 0.873137, 0.828032, 1.06541, 0.835574},
    {0.771414, 1.1288, 0.8, 1.2, 1.05198, 0.827995, 0.828032, 1.06848, 0.824145},
    {0.771414, 1.10228, 1.0, 1.2, 1.06677, 0.782854, 0.828032, 1.07154, 0.812715},
    {0.771414, 1.07576, 1.2, 1.2, 1.08156, 0.737707, 0.828032, 1.0746, 0.801295},
    {0.771414, 1.04924, 1.4, 1.2, 1.09635, 0.692565, 0.828032, 1.07767, 0.789866},
    {0.771414, 1.02272, 1.6, 1.2, 1.11114, 0.647425, 0.828032, 1.08073, 0.778436},
    {0.771414, 0.996193, 1.8, 1.2, 1.12593, 0.602278, 0.828032, 1.0838, 0.767014},
    {0.75283, 1.21669, 0.2, 1.4, 1.00325, 0.963425, 0.838325, 1.05801, 0.858425},
    {0.75283, 1.19016, 0.4, 1.4, 1.01804, 0.918284, 0.838325, 1.06107, 0.846995},
    {0.75283, 1.16364, 0.6, 1.4, 1.03283, 0.873137, 0.838325, 1.06414, 0.835574},
    {0.75283, 1.13712, 0.8, 1.4, 1.04762, 0.827995, 0.838325, 1.0672, 0.824145},
    {0.75283, 1.1106, 1.0, 1.4, 1.06241, 0.782854, 0.838325, 1.07027, 0.812715},
    {0.75283, 1.08407, 1.2, 1.4, 1.0772, 0.737707, 0.838325, 1.07333, 0.801295},
    {0.75283, 1.05755, 1.4, 1.4, 1.09198, 0.692565, 0.838325, 1.07639, 0.789866},
    {0.75283, 1.03103, 1.6, 1.4, 1.10677, 0.647425, 0.838325, 1.07946, 0.778436},
    {0.75283, 1.00451, 1.8, 1.4, 1.12156, 0.602278, 0.838325, 1.08252, 0.767014},
    {0.734247, 1.225, 0.2, 1.6, 0.998882, 0.963425, 0.848618, 1.05673, 0.858425},
    {0.734247, 1.19848, 0.4, 1.6, 1.01367, 0.918284, 0.848618, 1.0598, 0.846995},
    {0.734247, 1.17196, 0.6, 1.6, 1.02846, 0.873137, 0.848618, 1.06286, 0.835574},
    {0.734247, 1.14543, 0.8, 1.6, 1.04325, 0.827995, 0.848618, 1.06592, 0.824145},
    {0.734247, 1.11891, 1.0, 1.6, 1.05804, 0.782854, 0.848618, 1.06899, 0.812715},
    {0.734247, 1.09239, 1.2, 1.6, 1.07283, 0.737707, 0.848618, 1.07205, 0.801295},
    {0.734247, 1.06587, 1.4, 1.6, 1.08762, 0.692565, 0.848618, 1.07512, 0.789866},
    {0.734247, 1.03935, 1.6, 1.6, 1.10241, 0.647425, 0.848618, 1.07818, 0.778436},
    {0.734247, 1.01282, 1.8, 1.6, 1.1172, 0.602278, 0.848618, 1.08125, 0.767014},
    {0.715664, 1.23332, 0.2, 1.8, 0.994516, 0.963425, 0.858911, 1.05546, 0.858425},
    {0.715664, 1.20679, 0.4, 1.8, 1.0093, 0.918284, 0.858911, 1.05852, 0.846995},
    {0.715664, 1.18027, 0.6, 1.8, 1.02409, 0.873137, 0.858911, 1.06158, 0.835574},
    {0.715664, 1.15375, 0.8, 1.8, 1.03888, 0.827995, 0.858911, 1.06465, 0.824145},
    {0.715664, 1.12723, 1.0, 1.8, 1.05367, 0.782854, 0.858911, 1.06771, 0.812715},
    {0.715664, 1.1007, 1.2, 1.8, 1.06846, 0.737707, 0.858911, 1.07078, 0.801295},
    {0.715664, 1.07418, 1.4, 1.8, 1.08325, 0.692565, 0.858911, 1.07384, 0.789866},
    {0.715664, 1.04766, 1.6, 1.8, 1.09804, 0.647425, 0.858911, 1.07691, 0.778436},
    {0.715664, 1.02114, 1.8, 1.8, 1.11283, 0.602278, 0.858911, 1.07997, 0.767014}
};


    // Loop over each set of scale factors and calculate chi-squared for each
    for (size_t i = 0; i < scale_factors.size(); ++i) {
        const auto& factors = scale_factors[i];

        TH1D* h3all_scaled = (TH1D*)h3all->Clone("h3all_scaled");
        TH1D* h3all_scaled2 = (TH1D*)h3all->Clone("h3all_scaled2");

        TH1D* h3gencsall_scaled = (TH1D*)h3gencsall->Clone("h3gencsall_scaled");
        TH1D* h3genudall_scaled = (TH1D*)h3genudall->Clone("h3genudall_scaled");
        TH1D* h3genxall_scaled = (TH1D*)h3genxall->Clone("h3genxall_scaled");

        TH1D* h3tagcsall_scaled = (TH1D*)h3tagcsall->Clone("h3tagcsall_scaled");
        TH1D* h3tagudall_scaled = (TH1D*)h3tagudall->Clone("h3tagudall_scaled");
        TH1D* h3tagxall_scaled = (TH1D*)h3tagxall->Clone("h3tagxall_scaled");

        TH1D* h3gencstagcs_scaled = (TH1D*)h3gencstagcs->Clone("h3gencstagcs_scaled");
        TH1D* h3gencstagud_scaled = (TH1D*)h3gencstagud->Clone("h3gencstagud_scaled");
        TH1D* h3gencstagx_scaled = (TH1D*)h3gencstagx->Clone("h3gencstagx_scaled");

        TH1D* h3genudtagcs_scaled = (TH1D*)h3genudtagcs->Clone("h3genudtagcs_scaled");
        TH1D* h3genudtagud_scaled = (TH1D*)h3genudtagud->Clone("h3genudtagud_scaled");
        TH1D* h3genudtagx_scaled = (TH1D*)h3genudtagx->Clone("h3genudtagx_scaled");

        TH1D* h3genxtagcs_scaled = (TH1D*)h3genxtagcs->Clone("h3genxtagcs_scaled");
        TH1D* h3genxtagud_scaled = (TH1D*)h3genxtagud->Clone("h3genxtagud_scaled");
        TH1D* h3genxtagx_scaled = (TH1D*)h3genxtagx->Clone("h3genxtagx_scaled");

        // Apply scaling factors
        h3gencstagcs_scaled->Scale(factors[0]);
        h3gencstagx_scaled->Scale(factors[1]);
        h3gencstagud_scaled->Scale(factors[2]);
        h3genudtagcs_scaled->Scale(factors[3]);
        h3genudtagx_scaled->Scale(factors[4]);
        h3genudtagud_scaled->Scale(factors[5]);
        h3genxtagcs_scaled->Scale(factors[6]);
        h3genxtagx_scaled->Scale(factors[7]);
        h3genxtagud_scaled->Scale(factors[8]);

        // Reset h3all_scaled to zero and sum the scaled histograms
        h3all_scaled->Reset();
        h3all_scaled->Add(h3gencstagcs_scaled);
        h3all_scaled->Add(h3gencstagud_scaled);
        h3all_scaled->Add(h3gencstagx_scaled);
        h3all_scaled->Add(h3genudtagcs_scaled);
        h3all_scaled->Add(h3genudtagud_scaled);
        h3all_scaled->Add(h3genudtagx_scaled);
        h3all_scaled->Add(h3genxtagcs_scaled);
        h3all_scaled->Add(h3genxtagud_scaled);
        h3all_scaled->Add(h3genxtagx_scaled);

        // Reset and sum the histograms for each category
        h3tagcsall_scaled->Reset();
        h3tagcsall_scaled->Add(h3gencstagcs_scaled);
        h3tagcsall_scaled->Add(h3genudtagcs_scaled);
        h3tagcsall_scaled->Add(h3genxtagcs_scaled);

        h3tagudall_scaled->Reset();
        h3tagudall_scaled->Add(h3gencstagud_scaled);
        h3tagudall_scaled->Add(h3genudtagud_scaled);
        h3tagudall_scaled->Add(h3genxtagud_scaled);

        h3tagxall_scaled->Reset();
        h3tagxall_scaled->Add(h3gencstagx_scaled);
        h3tagxall_scaled->Add(h3genudtagx_scaled);
        h3tagxall_scaled->Add(h3genxtagx_scaled);

        // Final scaling based on the integral of h3all_scaled
        double scale_factor2 = 1.0 / h3all_scaled->Integral();
        h3gencstagx_scaled->Scale(scale_factor2);
        h3genudtagx_scaled->Scale(scale_factor2);
        h3genxtagx_scaled->Scale(scale_factor2);
        h3gencstagud_scaled->Scale(scale_factor2);
        h3genudtagud_scaled->Scale(scale_factor2);
        h3genxtagud_scaled->Scale(scale_factor2);
        h3gencstagcs_scaled->Scale(scale_factor2);
        h3genudtagcs_scaled->Scale(scale_factor2);
        h3genxtagcs_scaled->Scale(scale_factor2);
        h3tagcsall_scaled->Scale(scale_factor2);
        h3tagudall_scaled->Scale(scale_factor2);
        h3tagxall_scaled->Scale(scale_factor2);
        h3all_scaled->Scale(scale_factor2);

        // Calculate chi-squared for each category
        TH1D *h_tagcsvsdata_scaled = (TH1D*)h3tagcsall_scaled->Clone("h_tagcsvsdata_scaled");
        h_tagcsvsdata_scaled->Divide(h3tagcsall_data);
        TF1* f7m = new TF1("f7m", "[0]", 70, 100);
        f7m->FixParameter(0, 1);
        h_tagcsvsdata_scaled->Fit(f7m, "RN");
        double chi27m = f7m->GetChisquare();
        chi2_cs_values.push_back(chi27m);  // Store the chi-squared for 'cs'


        TH1D *h_tagudvsdata_scaled = (TH1D*)h3tagudall_scaled->Clone("h_tagudvsdata_scaled");
        h_tagudvsdata_scaled->Divide(h3tagudall_data);
        TF1* f6m = new TF1("f6m", "[0]", 70, 100);
        f6m->FixParameter(0, 1);
        h_tagudvsdata_scaled->Fit(f6m, "RN");
        double chi26m = f6m->GetChisquare();
        chi2_ud_values.push_back(chi26m);  // Store the chi-squared for 'ud'

        TH1D *h_tagxvsdata_scaled = (TH1D*)h3tagxall_scaled->Clone("h_tagxvsdata_scaled");
        h_tagxvsdata_scaled->Divide(h3tagxall_data);
        TF1* f5m = new TF1("f5m", "[0]", 70, 100);
        f5m->FixParameter(0, 1);
        h_tagxvsdata_scaled->Fit(f5m, "RN");
        double chi25m = f5m->GetChisquare();
        chi2_x_values.push_back(chi25m);  // Store the chi-squared for 'x'

        TH1D *h_allvsdata_scaled = (TH1D*)h3all_scaled->Clone("h_allvsdata_scaled");
        h_allvsdata_scaled->Divide(h3all_data);
        TF1* f8m = new TF1("f8m", "[0]", 70, 100);
        f8m->FixParameter(0, 1);
        h_allvsdata_scaled->Fit(f8m, "RN");
        double chi28m = f8m->GetChisquare();
        chi2_all_values.push_back(chi28m);  // Store the chi-squared for 'all'

        // Clean up memory if necessary
        delete h3all_scaled;
        delete h3all_scaled2;
        delete h3gencsall_scaled;
        delete h3genudall_scaled;
        delete h3genxall_scaled;
        delete h3tagcsall_scaled;
        delete h3tagudall_scaled;
        delete h3tagxall_scaled;
        delete h3gencstagcs_scaled;
        delete h3gencstagud_scaled;
        delete h3gencstagx_scaled;
        delete h3genudtagcs_scaled;
        delete h3genudtagud_scaled;
        delete h3genudtagx_scaled;
        delete h3genxtagcs_scaled;
        delete h3genxtagud_scaled;
        delete h3genxtagx_scaled;
    }
    // Now print the chi-squared values in the array format
std::cout << "double cs[n] = {";
for (size_t i = 0; i < chi2_cs_values.size(); ++i) {
    std::cout << chi2_cs_values[i];
    if (i != chi2_cs_values.size() - 1) std::cout << ", ";
}
std::cout << "};" << std::endl;

std::cout << "double ud[n] = {";
for (size_t i = 0; i < chi2_ud_values.size(); ++i) {
    std::cout << chi2_ud_values[i];
    if (i != chi2_ud_values.size() - 1) std::cout << ", ";
}
std::cout << "};" << std::endl;

std::cout << "double x[n] = {";
for (size_t i = 0; i < chi2_x_values.size(); ++i) {
    std::cout << chi2_x_values[i];
    if (i != chi2_x_values.size() - 1) std::cout << ", ";
}
std::cout << "};" << std::endl;

std::cout << "double all[n] = {";
for (size_t i = 0; i < chi2_all_values.size(); ++i) {
    std::cout << chi2_all_values[i];
    if (i != chi2_all_values.size() - 1) std::cout << ", ";
}
std::cout << "};" << std::endl;
}