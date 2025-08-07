#ifndef JERHANDLER_H
#define JERHANDLER_H

#include <string>
#include <iostream>
#include <random>
#include <memory>
#include <cmath>
#include <algorithm>  // for std::max
#include "JetResolution.h"
#include "JetResolutionObject.h"

class JERHandler {
private:
    std::shared_ptr<JME::JetResolution> resolution;
    std::shared_ptr<JME::JetResolutionScaleFactor> scaleFactor;
    std::mt19937 randomGenerator;
    std::normal_distribution<> gauss;

    std::string currentPeriod;
    bool initialized;

    std::map<std::string, std::string> periodToSFfile = {
        {"Muo16",    "JER_SF/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt"},
        {"Muo16APV", "JER_SF/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt"},
        {"Muo17",    "JER_SF/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt"},
        {"Muo18",    "JER_SF/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt"}
    };

    std::map<std::string, std::string> periodToResolutionFile = {
        {"Muo16",    "JER_SF/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt"},
        {"Muo16APV", "JER_SF/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt"},
        {"Muo17",    "JER_SF/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt"},
        {"Muo18",    "JER_SF/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt"}
    };

    std::string basePath = "./";

public:
    JERHandler()
    : randomGenerator(37428479), gauss(0, 1.0), initialized(false)
    {}

    void setBasePath(const std::string& path) {
        basePath = path;
    }

    void update(const std::string& filename) {
        std::string period = determinePeriod(filename);
        if (period.empty()) {
            std::cerr << "Unknown period for file: " << filename << std::endl;
            exit(1);
        }

        if (period != currentPeriod) {
            std::cout << "Processing file: " << filename << std::endl;
            std::string base = basePath;
            std::string sfFile = base + periodToSFfile[period];
            std::string resFile = base + periodToResolutionFile[period];
            std::cout << "Loading scale factor file: " << sfFile << std::endl;
            std::cout << "Loading resolution file: " << resFile << std::endl;
            resolution = std::make_shared<JME::JetResolution>(resFile);
            scaleFactor = std::make_shared<JME::JetResolutionScaleFactor>(sfFile);
            currentPeriod = period;
            initialized = true;
        }
    }

    double smearPt(double pt, double eta, double rho) {
        if (!initialized) {
            std::cerr << "JERHandler not initialized!" << std::endl;
            exit(1);
        }
    
        JME::JetParameters parameters;
        parameters.setJetPt(pt).setJetEta(eta).setRho(rho);
    
        float sigma = resolution->getResolution(parameters);
        float sf = 1.1; //scaleFactor->getScaleFactor(parameters, JME::Variation::NOMINAL);
    
        std::normal_distribution<> smearDist(0.0, sigma); // N(0, σ_JER)
        double fluctuation = smearDist(randomGenerator);  // draw once
        double smearFactor = 1.0 + fluctuation * std::sqrt(sf * sf - 1.0);
        
        if (pt * smearFactor < 1e-2) {
            smearFactor = 1e-2 / pt;
        }
    
        double smearedPt = pt * smearFactor;
        return smearedPt;
    }

    double getScaleFactor(double pt, double eta, double rho) const {
        JME::JetParameters parameters;
        parameters.setJetPt(pt).setJetEta(eta).setRho(rho);
    
        float sf = scaleFactor->getScaleFactor(parameters, JME::Variation::NOMINAL);
        return sf;
    }

private:
    std::string determinePeriod(const std::string& filename) {
        if (filename.find("Muo16APV") != std::string::npos) return "Muo16APV";
        if (filename.find("Muo16")    != std::string::npos) return "Muo16";
        if (filename.find("Muo17")    != std::string::npos) return "Muo17";
        if (filename.find("Muo18")    != std::string::npos) return "Muo18";
        return "";
    }
};

#endif