#pragma once

#include <string>
#include <vector>

/// Configuration parameters loaded from YAML
struct Config {
    // General simulation settings
    int   nEvents;

    // Negative Binomial distribution parameters
    std::vector<double> mu5TeV;
    std::vector<double> k5TeV;
    std::vector<double> NBDLow;
    std::vector<double> NBDHigh;
    std::vector<double> fracLBC;
    // Kinetic freeze-out temperature per centrality bin
    std::vector<double> Tkin;
    // Blast-Wave flow parameters
    std::vector<double> betaT;
    std::vector<double> n;
    std::vector<double> rho2_p;
    std::vector<double> rho2_L;
    std::vector<double> Rx;

    // Particle ratio & debug flag
    double ratioProtonLambda;
    double ratioProtonInclusive;
    bool   isDebug;

    /// Load configuration from a YAML file
    static Config Load(const std::string& filename);
};
