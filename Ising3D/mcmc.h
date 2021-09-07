//
// Created by kamilla on 13.10.2019.
//
#ifndef MC_CPP_MCMC_H
#define MC_CPP_MCMC_H
//
// Created by kamilla on 13.10.2019.
//
//#ifndef MC_CPP_MCMC_H
//#define MC_CPP_MCMC_H
#include <list>
#include <valarray>
#include <algorithm>
#include <random>
#include <tuple>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <math.h>
#include <set>
#include "observable.h"

typedef long int coord_t;

class Lattice {
public:
    long int lattice_side;
    int ndim() { return 2; }
    //int ndim2() { return 4; }
    int ndim3() { return 6; }
    Lattice(long int max_seq_size = 10);
    void create_lattice(long int max_seq_size = 100);
    //private:
    std::valarray<long int> map_of_contacts_int;
    std::valarray<long int> x_coords;
    std::valarray<long int> y_coords;
    std::valarray<long int> z_coords;

};
class Protein {
public:

    Protein();
    Protein(long int length);

    void Reconnect(int j);
    void Reconnect1(int j);
    bool IsEndInStuck();

    void MC(double J_in = 0, double h_in = 0, int nSimulation = 0, long int steps_to_equilibrium = 400000000, long int mc_steps = 5000000000000000, bool bradius = false);
    //void MC(  double J=0, double h=0, int nSimulation = 0, long int steps_to_equilibrium = 1000000, long int mc_steps = 100000000, bool radius = false);

    void save_calcs();
    void radius();
    void radius_gyration();
    void radius_gyration1();

    void count_contacts();
    bool CheckAndFlipNodes(long int& coord, int& sign);

    void calc_bulk();

public:

    Lattice lattice;

    std::valarray<long int> ordered_coords;

    std::valarray<int> directions; //�� n-1;
    //���� � directions[10] ����� 0, ������, ��������� ������ �� ���������� 10

    std::valarray< int> sequence_on_lattice;
    std::valarray<long int> next_monomers;
    std::valarray<long int> previous_monomers;
    long int end_conformation = 0, start_conformation = 0;

    std::queue<long int>  spins_in_cluster;

    mc_stats::ScalarObservable<double> dists;
    mc_stats::ScalarObservable<double> gyration;

    mc_stats::ScalarObservable<long double> energy; //���������� �������
    mc_stats::ScalarObservable<long double> energy_sq;
    mc_stats::ScalarObservable<long double> energy_4;

    mc_stats::ScalarObservable<long double> magnetization; //���������� ���������������
    mc_stats::ScalarObservable<long double> magnetization_sq;
    mc_stats::ScalarObservable<long double> magnetization_4;

    mc_stats::ScalarObservable<double> bulk6;
    mc_stats::ScalarObservable<double> bulk5;
    mc_stats::ScalarObservable<double> bulk4; //доля узлов с 4 соседами
    mc_stats::ScalarObservable<double> bulk3;
    mc_stats::ScalarObservable<double> bulk2;

    long int number_of_monomers = 0;
    long int E = 0; // -1* ����� �������������� ��������� ������� �����������
    long int current_H_counts = 0;
    int bulk6_now = 0, bulk5_now = 0, bulk4_now = 0, bulk3_now = 0, bulk2_now = 0;

    std::map <long int, long long int> count_E;
    std::map <long int, long long int> count_M;

    long double sum_X = 0, sum_Y = 0, sum_Z = 0;
    double P_add = 0;

private:
    double J = 0;
    double h = 0;
    int nSimulation = 0;

};


#endif //MC_CPP_MCMC_H