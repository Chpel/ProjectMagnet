#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "Ising.h"
#include <stdlib.h>

int main(int argc, char* argv[]) {
    int size = atoi(argv[1]);
    double CoefXY = atof(argv[2]), T = 2.269;
    bool OpenQ = true;
    std::string lType;
    if ((argc == 3) || (argc == 1)) {
        lType = "open";
        OpenQ = true;
    }
    else {
        lType = argv[3];
        if ((lType == "obc") || (lType == "pbc")) {
            OpenQ = (lType == "obc");
            if (OpenQ) {
                lType = "open";
            }
            else {
                lType = "period";
            }
        }
        else {
            lType = "open";
            OpenQ = true;
        }
    }
    //int size = 50;
    //double CoefXY = 0.5, T = 2.269;
    int kIt0 = 1000000, kIt1 = 10000, kIt2 = 100;
    int iCheck = 0;
    std::vector<double> m2;
    std::vector<double> m4;
    Ising ising(size, CoefXY);
    double p = 1 - exp(-(1 / T) * 2);
    std::cout << "The program has started: size = " << size << ", CoefXY = " << CoefXY << ", OpenQ = " << OpenQ << std::endl;
    for (int i = 0; i < kIt0; i++) {
        ising.cluster_ising(p, lType);
        iCheck++;
        if (iCheck == kIt0 / 10) {
            std::cout << "#";
            iCheck = 0;
        }
    }
    std::cout << std::endl;
    std::cout << "Pre-simulations are done." << std::endl;
    iCheck = 0;
    for (int i = 0; i < kIt1; i++) {
        for (int j = 0; j < kIt2; j++) {
            ising.cluster_ising(p, lType);
        }
        m2.push_back(pow(ising.magnetism(), 2));
        m4.push_back(pow(ising.magnetism(), 4));
        iCheck++;
        if (iCheck == kIt1 / 10) {
            std::cout << "#";
            iCheck = 0;
        }
    }
    std::cout << std::endl;
    std::cout << "All simulations are done. Analysis of observable measures starts..." << std::endl;
    double avM2, errM2, avM4, errM4;
    bool cov;
    std::tie(avM2, errM2, cov) = detail::Analysis(m2);
    std::tie(avM4, errM4, cov) = detail::Analysis(m4);
    std::ofstream file(std::string("results") + std::string(".txt"), std::ios_base::app);
    file << size << " " << CoefXY << " " << OpenQ << " " << avM2 << " " << errM2 << " " << avM4 << " " << errM4 << " " << 1 - avM4 / (3 * pow(avM2, 2)) << std::endl;
    file.close();
    std::cout << "Done!" << std::endl;
    return 0;
}
