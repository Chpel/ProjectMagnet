#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "Ising.h"
#include <stdlib.h>

void logger(int L, double J, bool OpenQ, double CoefXY, std::string logText) {
    std::cout << "For size = " << L << ", J = " << J << ", OpenQ = " << OpenQ << " , CoefXY = " << CoefXY << ": " << logText << std::endl;
}

void researcher(int L, double J, bool OpenQ, double CoefXY, int steps, std::vector<double>& m, std::vector<double>& m2, std::vector<double>& m4, std::vector<double>& e, std::vector<double>& e2) {
    double avM, errM, avM2, errM2, avM4, errM4, avE, errE, avE2, errE2;
    bool cov;
    std::tie(avM, errM, cov) = detail::Analysis(m);
    std::tie(avM2, errM2, cov) = detail::Analysis(m2);
    std::tie(avM4, errM4, cov) = detail::Analysis(m4);
    std::tie(avE, errE, cov) = detail::Analysis(e);
    std::tie(avE2, errE2, cov) = detail::Analysis(e2);
    std::ofstream file(std::string("results_") + std::to_string(L) + "_" + std::to_string(J) + "_" + std::to_string(CoefXY) + "_" + std::to_string(OpenQ) + std::string(".txt"), std::ios_base::trunc);
    file << "L  J  CoefXY  OpenQ  avM  errM  avM2  errM2  avM4  errM4  avE  errE  avE2  errE2  U4  steps" << std::endl;
    file << L << " " << J << " " << CoefXY << " " << OpenQ << " " << avM << " " << errM << " " << avM2 << " " << errM2 << " " << avM4 << " " << errM4 << " " << avE << " " << errE << " " << avE2 << " " << errE2 << " " << 1 - avM4 / (3 * pow(avM2, 2)) << " " << steps << std::endl;
    file.close();
}

int main(int argc, char* argv[]) {
    int size = atoi(argv[1]);
    double T = 1 / atof(argv[2]);
    bool OpenQ;
    std::string lType;
    double CoefXY;
    if (argc < 4) {
        lType = "open";
        OpenQ = true;
    }
    else if (std::string(argv[3]) == "pbc") {
        lType = "period";
        OpenQ = false;
    }
    else {
        lType = "open";
        OpenQ = true;
    }
    if (argc < 5) {
        CoefXY = 1;
    }
    else CoefXY = atof(argv[4]);
    //int size = 50;
    //double CoefXY = 0.5, T = 2.269;
    int kIt0 = 1000000, kIt1 = 100000, kIt2 = 100;
    // 
    // 
    // 
    //int kIt0 = 10, kIt1 = 200, kIt2 = 1;
    std::vector<double> m;
    std::vector<double> m2;
    std::vector<double> m4;
    std::vector<double> e;
    std::vector<double> e2;
    Ising ising(size, CoefXY);
    double p = 1 - exp(-(1 / T) * 2);
    logger(size, 1 / T, OpenQ, CoefXY, "program has started!");
    for (int i = 0; i < kIt0; i++) {
        ising.cluster_ising(p, lType);
    }
    logger(size, 1 / T, OpenQ, CoefXY, "presimulations are done!");
    int checker = 0;
    for (int i = 0; i < kIt1; i++) {
        for (int j = 0; j < kIt2; j++) {
            ising.cluster_ising(p, lType);
        }
        e.push_back(ising.energy(1 / T, lType));
        e2.push_back(pow(ising.energy(1 / T, lType),2));
        m.push_back(ising.magnetism());
        m2.push_back(pow(ising.magnetism(), 2));
        m4.push_back(pow(ising.magnetism(), 4));
        checker++;
        if (checker == 100) {
            researcher(size, 1 / T, OpenQ, CoefXY, kIt0 + i*kIt2, m, m2, m4, e, e2);
            checker = 0;
        }
    }
    logger(size, 1 / T, OpenQ, CoefXY, "Done!");
    return 0;
}