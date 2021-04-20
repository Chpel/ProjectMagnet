#pragma once

#include <vector>
#include <string>
#include <random>
#include <queue>
#include <set>
#include <cmath>
#include <utility>
#include <tuple>
#include <algorithm>


class Ising {
public:
    Ising();
    Ising(int sz, double CoefXY=1);
    void cluster_ising(double p, std::string lattice_type="open");
    std::vector<std::pair<int, int>> get_neighbors(const std::pair<int, int>& index, std::string lattice_type="open") const;
    std::vector<std::vector<int>> get_lattice() {return lattice;};
    double magnetism();
private:
    int size;
    int Xsize;
    int Ysize;
    std::random_device rd;
    std::vector<std::vector<int>> lattice;
};

namespace detail {
    const size_t N_B_MAX = 1000;
    template<typename T> void collate(std::vector<T>& arr);
    template<typename T> std::pair<T, T> bSTAT(const std::vector<T>& blocks);
    template<typename T> std::tuple< std::vector<T>, std::vector<T>, std::vector<T> > mrg(const std::vector<T>&);
    template<typename T> std::tuple<T, T, bool> block_stats(const std::vector<T>&, const std::vector<T>);
}



namespace detail {

    /********************************************************
     * Merge blocks in-place: 100 -> 50 twice larger blocks *
     ********************************************************/
    template<typename T>
    void
        collate(std::vector<T>& arr) {
        size_t n2 = floor(arr.size() / 2);
        for (size_t j = 0; j < n2; j++) {
            arr[j] = 0.5 * (arr[2 * j] + arr[2 * j + 1]);
        }
        arr.resize(n2);
    }


    /****************************************
     * Block statistics at fixed block size *
     ****************************************/
    template<typename T>
    std::pair<T, T>
        bSTAT(const std::vector<T>& blocks) {
        T av = 0, av2 = 0;
        size_t n_b = blocks.size();
        for (size_t j = 0; j < n_b; j++) {
            av += blocks[j] / n_b;
            av2 += blocks[j] * blocks[j] / n_b;
        }
        T dif = av2 - av * av;
        if (dif < 0) { dif = 0.; }
        T err = sqrt(dif) / sqrt(1.0 * n_b);
        return std::make_pair(av, err);
    }


    /******************************************
     * Merge blocks, analyze block statistics *
     ******************************************/
    template<typename T>
    std::tuple< std::vector<T>, std::vector<T>, std::vector<T> >
        mrg(const std::vector<T>& blocks) {
        std::vector<T> arr;
        std::copy(blocks.begin(), blocks.end(), back_inserter(arr));

        std::vector<T> v_av, v_err, v_size;

        do {

            std::pair<T, T> av_err = bSTAT(arr);
            v_av.push_back(av_err.first);
            v_err.push_back(av_err.second);
            v_size.push_back(arr.size());

            collate(arr);

        } while (arr.size() > 4);

        return std::tie(v_av, v_err, v_size);
    }


    template<typename T>
    std::tuple<T, T, bool>
        block_stats(const std::vector<T>& v_av, const std::vector<T> v_err) {

        // block mean is the second-to last entry
        T av = v_av.size() > 2 ? v_av[v_av.size() - 2]
            : v_av[0];
        T err = *std::max_element(v_err.begin(), v_err.end());

        // TODO: convergence check
        bool conv = false;

        return std::tie(av, err, conv);
    }

    template<typename T>
    std::tuple<T, T, bool>
        Analysis(const std::vector<T>& _blocks) {
        T av, err;
        bool conv;
        std::vector<T> v_av, v_err, v_size;
        std::tie(v_av, v_err, v_size) = mrg(_blocks);
        return block_stats(v_av, v_err);
    }
}  // namespace detail

