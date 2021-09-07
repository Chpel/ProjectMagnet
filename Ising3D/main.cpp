#include <iostream>
#include <cstdlib>
#include"mcmc.h"
#include <string>

int main(int argc, char* argv[]) {
    std::string par[3] = { "10", "100", "1" };
    //std::cout << "Hello, World!" << std::endl;
    int N = std::atoi(argv[1]);
    double J = 0.01 * std::atof(argv[2]);
    double  h = 0.1*(double)std::stoi(argv[3]);
    std::cout << J << " " << N << " " << h << " started!" << std::endl;

    //int nSim = std::stoi(par[2]);
    time_t start, end;
    time(&start);

    Protein p(N);

    p.MC(J, h);
    //Protein p(20);
    //p.MC(0.8);


    time(&end);

    double dif = difftime(end, start);
    printf("For N = %d, J = %lf, h = %lf, Time = %lf \n", N, J, h, dif);

    /*for (int i = 4; i < 15 ; i++ )
    {
        Protein p(i);
        p.MC(0.5);
    }*/

    //p.MC(0.8);
    return 0;
}