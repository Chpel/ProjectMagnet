#include "Ising.h"

Ising::Ising() {
    Xsize = 30;
    Ysize = 30;
    for (int i = 0; i <Ysize; i++) {
        lattice.emplace_back(Xsize, 1);
    }

}

Ising::Ising(int sz, double CoefXY) {
    size = sz;
    Xsize = sz;
    Ysize = (int)(sz * CoefXY);
    for (int i = 0; i < Ysize; i++) {
        lattice.emplace_back(Xsize, 1);
    }
}

std::vector<std::pair<int, int>> Ising::get_neighbors(const std::pair<int, int>& index, std::string lattice_type) const{
    std::vector<std::pair<int, int>> neig;
    if (lattice_type == "period") {
        if ((index.first + 1) < Ysize)
            neig.emplace_back(index.first + 1, index.second);
        else
            neig.emplace_back(0, index.second);
        if ((index.first - 1) >= 0)
            neig.emplace_back(index.first - 1, index.second);
        else
            neig.emplace_back(Ysize - 1, index.second);
        if ((index.second + 1) < Xsize)
            neig.emplace_back(index.first, index.second + 1);
        else
            neig.emplace_back(index.first, 0);
        if ((index.second - 1) >= 0)
            neig.emplace_back(index.first, index.second - 1);
        else
            neig.emplace_back(index.first, Xsize - 1);
    }
    if (lattice_type == "open") {
        if ((index.first + 1) < Ysize)
            neig.emplace_back(index.first + 1, index.second);
        if ((index.first - 1) >= 0)
            neig.emplace_back(index.first - 1, index.second);
        if ((index.second + 1) < Xsize)
            neig.emplace_back(index.first, index.second + 1);
        if ((index.second - 1) >= 0)
            neig.emplace_back(index.first, index.second - 1);
    }
    return std::move(neig);
}

void Ising::cluster_ising(double p, std::string lattice_type) {
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis1(0, Ysize - 1);
    std::uniform_int_distribution<> dis2(0, Xsize - 1);
    std::pair<int, int> curr, j;

    j.first = dis1(gen);
    j.second = dis2(gen);

    std::vector<std::vector<bool>> used;
    for (int i = 0; i < Ysize; i++) {
        used.emplace_back(Xsize, false);
    }
    std::vector<std::pair<int, int>> neighbors;
    std::queue<std::pair<int, int>> q;

    q.push(j);
    used[j.first][j.second] = true;

    while (!q.empty()) {
        curr = q.front();
        q.pop();
        if (lattice_type == "open") {
            neighbors = get_neighbors(curr, "open");
        } else {
            neighbors = get_neighbors(curr, "period");
        }
        for (auto neighbor: neighbors) {
            auto rand = std::generate_canonical<double, 30>(gen);
            //???
            if (lattice[j.first][j.second] == lattice[neighbor.first][neighbor.second] && rand < p &&
            !used[neighbor.first][neighbor.second]) {
                q.push(neighbor);
                used[neighbor.first][neighbor.second] = true;
                lattice[neighbor.first][neighbor.second] *= -1;
            }
        }
    }
    lattice[j.first][j.second] *= -1;

}


double Ising::magnetism() {
    int s = 0;
    for (auto row : lattice) {
        for (auto spin : row) {
            s += spin;
        }
    }
    return s*1.0 / (Xsize * Ysize);
}