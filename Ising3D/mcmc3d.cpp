#include "mcmc.h"
#include <iostream>
#include <fstream>
#include <cassert>

//Lattice::Lattice() {};
//std::random_device generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);
std::mt19937 generator(1234);
//std::random_device generators1;
std::uniform_int_distribution<int> distribution1(0, 5);
std::mt19937 generators1(123);
//std::random_device generators2;
std::uniform_int_distribution<int> distribution2(0, 1);
std::mt19937 generators3(12377);
//std::random_device generators3;

int ndim3() { return 6; }

Lattice::Lattice(long int max_seq_size) {
    lattice_side = max_seq_size;
    //создается одномерный массив соседей на квадратной решетке
    long int x, y, z;
    ldiv_t n;
    map_of_contacts_int.resize(lattice_side * lattice_side * lattice_side * ndim3());
    for (long int i = 0; i < lattice_side * lattice_side * lattice_side; i++) {
        map_of_contacts_int[ndim3() * i] = i + 1;
        map_of_contacts_int[ndim3() * i + 1] = i - 1;
        map_of_contacts_int[ndim3() * i + 2] = i + lattice_side;
        map_of_contacts_int[ndim3() * i + 3] = i - lattice_side;
        //добавим к этому соседей из смежных плоскостей:
        map_of_contacts_int[ndim3() * i + 4] = i + lattice_side * lattice_side;
        map_of_contacts_int[ndim3() * i + 5] = i - lattice_side * lattice_side;
        //у нас появляется z: i = z * l_s^2 + y * l_s + x
        //поэтому:
        n = div(i, lattice_side * lattice_side);
        z = n.quot;
        n = div(n.rem, lattice_side);
        x = n.rem;
        y = n.quot;
        for (int j = 0; j < ndim3(); j++) {
            if (x == 0) {
                map_of_contacts_int[6 * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[6 * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                map_of_contacts_int[6 * i + 3] = lattice_side * (lattice_side - 1) + x;
            }
            if (y == (lattice_side - 1)) {
                map_of_contacts_int[6 * i + 2] = x;
            }
            //и случай для крайних плоскостей:
            if (z == 0) {
                map_of_contacts_int[6 * i + 5] = x + lattice_side * y + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                map_of_contacts_int[6 * i + 4] = x + lattice_side * y;
            }
        }
    }
}

void Lattice::create_lattice(long int max_seq_size) {

    lattice_side = max_seq_size;
    long int x, y, z;
    ldiv_t n;
    map_of_contacts_int.resize(lattice_side * lattice_side * lattice_side * ndim3());
    x_coords.resize(lattice_side * lattice_side * lattice_side * ndim3());
    y_coords.resize(lattice_side * lattice_side * lattice_side * ndim3());
    z_coords.resize(lattice_side * lattice_side * lattice_side * ndim3());
    for (long int i = 0; i < lattice_side * lattice_side * lattice_side; i++) {
        map_of_contacts_int[ndim3() * i] = i + 1;
        map_of_contacts_int[ndim3() * i + 1] = i - 1;
        map_of_contacts_int[ndim3() * i + 2] = i + lattice_side;
        map_of_contacts_int[ndim3() * i + 3] = i - lattice_side;
        //добавим к этому соседей из смежных плоскостей:
        map_of_contacts_int[ndim3() * i + 4] = i + lattice_side * lattice_side;
        map_of_contacts_int[ndim3() * i + 5] = i - lattice_side * lattice_side;
        //у нас появляется z: i = z * l_s^2 + y * l_s + x
        //поэтому:
        n = div(i, lattice_side * lattice_side);
        z = n.quot;
        n = div(n.rem, lattice_side);
        x = n.rem;
        y = n.quot;
        for (int j = 0; j < ndim3(); j++) {
            if (x == 0) {
                map_of_contacts_int[6 * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[6 * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                map_of_contacts_int[6 * i + 3] = lattice_side * (lattice_side - 1) + x;
            }
            if (y == (lattice_side - 1)) {
                map_of_contacts_int[6 * i + 2] = x;
            }
            //и случай для крайних плоскостей:
            if (z == 0) {
                map_of_contacts_int[6 * i + 5] = x + lattice_side * y + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                map_of_contacts_int[6 * i + 4] = x + lattice_side * y;
            }
        }
        x_coords[i] = x;
        y_coords[i] = y;
        z_coords[i] = z;


    }

}

Protein::Protein() {}

Protein::Protein(long int n) {

    /** type | previous | next   **/
    lattice.create_lattice(n + 5); //создание решетки, на 5 больше, чем длина цепочки
    // 0 - отстуттвие элемента на решетке
    sequence_on_lattice.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, 0); //последовательность мономеров
    next_monomers.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //номер-ссылка на следующий узел
    previous_monomers.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //номер-ссылка на предыдующий узел

    directions.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //направления из {0,1,2,3}

    ordered_coords.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //пока что для последовательной нумерации

    number_of_monomers = n;
    start_conformation = 0;
    end_conformation = n - 1;

    sum_X = 0;
    sum_Y = 0;
    sum_Z = 0;
    for (int i = 1; i < n - 1; i++)
    {
        previous_monomers[i] = i - 1;
        sequence_on_lattice[i] = 1;
        next_monomers[i] = i + 1;

        ordered_coords[i] = i;

        sum_X += i;
    }
    sum_X = sum_X + n - 1;
    ordered_coords[0] = 0;
    ordered_coords[n - 1] = 0;

    sequence_on_lattice[0] = 1;
    sequence_on_lattice[end_conformation] = 1; //начальная последовательность

    next_monomers[0] = 1;
    previous_monomers[n - 1] = n - 2;

    current_H_counts = n;
    E = -(n - 1);


    //сначала все направления - движение вправо
    for (int i = 0; i < n - 1; i++)
    {
        directions[i] = 0;
    }




}


bool Protein::IsEndInStuck()
{
    int hh = 0;
    coord_t  step;
    for (int j = 0; j < lattice.ndim3(); j++) {
        step = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + j];
        if (sequence_on_lattice[step] != 0) {
            hh = hh + 1;
        }
    }
    return hh == lattice.ndim3();
}

void Protein::Reconnect1(int j) {
    int inverse_steps[6] = { 1, 0, 3, 2, 5, 4 };
    int reflect_directions[7][6] =
    { {3, 2, 0, 1, 4, 5}, //90, 0
     {1, 0, 3, 2, 4, 5}, //180, 0
     {2, 3, 1, 0, 4, 5}, //270, 0
     {0, 1, 2, 3, 4, 5}, //360, 0
     {5, 4, 2, 3, 0, 1}, //0, 90
     {1, 0, 2, 3, 5, 4}, //0, 180
     {4, 5, 2, 3, 1, 0} //0, 270
        //возможны ещё, но пока недописал
    };

    //int j = directions[previous_monomers[end_conformation]]; //направление последнего ребра

    long int step = lattice.map_of_contacts_int[lattice.ndim2() * end_conformation + j];
    long int new_end = next_monomers[step];

    next_monomers[step] = end_conformation;

    directions[step] = inverse_steps[j];
    directions[new_end] = -1;

    long int c = end_conformation;
    long int new_c;
    while (c != new_end)
    {
        new_c = previous_monomers[c];

        next_monomers[c] = previous_monomers[c];
        //previous_monomers[new_c]=c;
        directions[c] = inverse_steps[directions[new_c]];
        c = new_c;
    }
    long int temp_prev_next = next_monomers[new_end];

    previous_monomers[end_conformation] = step;
    c = end_conformation;

    while (c != new_end)
    {
        new_c = next_monomers[c];
        previous_monomers[new_c] = c;
        c = new_c;
    }

    end_conformation = new_end;


    previous_monomers[new_end] = temp_prev_next;
    next_monomers[new_end] = -1;

}



void Protein::Reconnect(int j) {
    int inverse_steps[6] = { 1, 0, 3, 2, 5, 4 };
    //    int reflect_directions[4][4] =
    //            {{3, 2, 0, 1}, //90
    //             {1, 0, 3, 2}, //180
    //             {2, 3, 1, 0}, //270
    //             {0, 1, 2, 3}
    //            };
    long int c;

    long int step = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + j];
    long int new_end = next_monomers[step];

    next_monomers[step] = end_conformation;

    directions[step] = inverse_steps[j];

    c = end_conformation;

    //std::cout << " end_conformation " <<  end_conformation << std::endl;
    if (end_conformation > lattice.lattice_side * lattice.lattice_side * lattice.lattice_side)
    {
        std::cout << " HZ why " << std::endl;
        return;

    }
    long int new_c;
    while (c != new_end)
    {
        // std::cout << "c  " << c << std::endl;
        assert(c >= 0 && c < lattice.lattice_side * lattice.lattice_side * lattice.lattice_side);
        if (c > lattice.lattice_side * lattice.lattice_side * lattice.lattice_side)
        {
            std::cout << " HZ why " << std::endl;
            return;
        }

        if (previous_monomers[c] > lattice.lattice_side * lattice.lattice_side * lattice.lattice_side)
        {
            std::cout << "nor  prev " << std::endl;
        }

        new_c = previous_monomers[c];

        if (new_c < 0 || c > lattice.lattice_side * lattice.lattice_side * lattice.lattice_side)
        {
            std::cout << " we have problems " << std::endl;
        }

        next_monomers[c] = previous_monomers[c];


        directions[c] = inverse_steps[directions[new_c]];
        c = new_c;
    }
    long int temp_prev_next = next_monomers[new_end];


    previous_monomers[end_conformation] = step;
    c = end_conformation;


    while (c != new_end)
    {
        new_c = next_monomers[c];
        previous_monomers[new_c] = c;
        c = new_c;
    }

    end_conformation = new_end;


    //std::cout <<"mmmm" << " " << step << std::endl;



    previous_monomers[new_end] = temp_prev_next;
    next_monomers[new_end] = -1;
    directions[new_end] = -1;


}


void Protein::MC(double J_in, double h_in, int nSumulation, long int steps_to_equilibrium, long int mc_steps, bool bradius)
{
    J = J_in;
    h = h_in;
    //[i,j]: i-поворот, j-направление (против часовой)
    int reflect_directions[4][4] =
    { {3, 2, 0, 1 }, //90
     {1,0,3,2 }, //180
     {2,3, 1, 0 }, //270
    { 0, 1, 2, 3}
    };

    int inverse_steps[6] = { 1,0,3,2,5,4 }; //для сохранения направлений в апдейте "перенести конец в начало"
    //double step_rd; //Для выбора апдейта: обычный или реконнект
    double q_rd, p1, p_metropolis; //Для вероятности принятия шага
    int rand_path; // = distribution1(generators1); //выбирается направление: 0 - переставляем начало в конец
    double typeOfUpdate; //0 - простой; 1 - реконнект
    int step;
    int step_on_lattice;//выбор одного из соседей
    long int new_point;
    long int new_E, new_H;
    int hh;
    long int temp, del, oldspin;

    double p_for_reconnect_update = 0.5;

    long int all_steps = steps_to_equilibrium + mc_steps;

    for (long int i = 0; i < all_steps; i++) {
        //std::cout << "STEP : " << i << std::endl;
        typeOfUpdate = distribution(generator);
        if (typeOfUpdate < p_for_reconnect_update) {
            hh = 0;
            rand_path = distribution2(generators3);

            if (rand_path == 0) {//переставляем начало в конец

                step_on_lattice = distribution1(generators1);
                new_point = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[start_conformation];

                if (sequence_on_lattice[new_point] == 0) { //проверка, что в узле нет мономеров

                    //делаем апдейт

                    //добавляем в конец
                    next_monomers[end_conformation] = new_point;
                    sequence_on_lattice[new_point] = 3 * distribution2(generators3) - 1; //выбор спина //!!!!!!!!!
                    previous_monomers[new_point] = end_conformation;
                    end_conformation = new_point;

                    //удаляем начало
                    temp = start_conformation;
                    start_conformation = next_monomers[start_conformation];
                    next_monomers[temp] = -1;
                    previous_monomers[start_conformation] = -1;

                    //смотрим потери
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * temp + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh - sequence_on_lattice[temp] * sequence_on_lattice[step];
                        }
                    }

                    //смотрим выигрыш
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh + sequence_on_lattice[end_conformation] * sequence_on_lattice[step];
                        }
                    }

                    new_E = E + hh;
                    //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[start_conformation];
                    new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];

                    p1 = exp(-(-(new_E - E) * J - (new_H - current_H_counts) * h));
                    p_metropolis = std::min(1.0, p1);
                    q_rd = distribution(generator);
                    if (q_rd < p_metropolis) {
                        E = new_E;
                        current_H_counts = new_H;
                        sequence_on_lattice[temp] = 0; //делаю здесь, так как проще считать энергию(!!!)
                        sum_X = sum_X + lattice.x_coords[end_conformation] - lattice.x_coords[temp];
                        sum_Y = sum_Y + lattice.y_coords[end_conformation] - lattice.y_coords[temp];
                        sum_Z = sum_Z + lattice.z_coords[end_conformation] - lattice.z_coords[temp];
                        //корректируем информацию о направлениях
                        directions[temp] = -1;
                        directions[previous_monomers[end_conformation]] = step_on_lattice;

                    }
                    else {//отменяем изменения
                        //удаляем конец
                        del = end_conformation;
                        end_conformation = previous_monomers[end_conformation];
                        next_monomers[end_conformation] = -1;
                        previous_monomers[del] = -1;
                        sequence_on_lattice[del] = 0;

                        //возвращаем начало
                        previous_monomers[start_conformation] = temp;
                        next_monomers[temp] = start_conformation;
                        start_conformation = temp;
                        sequence_on_lattice[start_conformation] = oldspin;

                    }

                }
                else {//места нет, выходим из шага
                    //continue;
                }

            }
            else {//переставляем конец в начало

                step_on_lattice = distribution1(generators1);
                new_point = lattice.map_of_contacts_int[lattice.ndim3() * start_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[end_conformation];

                if (sequence_on_lattice[new_point] == 0) { //проверка, что в узле нет мономеров

                    //делаем апдейт
                    //добавляем в начало
                    previous_monomers[start_conformation] = new_point;
                    sequence_on_lattice[new_point] = 3 * distribution2(generators3) - 1; //выбор спина //!!!!!!!!!
                    next_monomers[new_point] = start_conformation;
                    start_conformation = new_point;

                    //удаляем конец
                    temp = end_conformation;
                    end_conformation = previous_monomers[end_conformation];
                    if (previous_monomers[end_conformation] < 0) {
                        std::cout << "problem update " << std::endl;
                    }
                    previous_monomers[temp] = -1;
                    next_monomers[end_conformation] = -1;

                    //смотрим потери
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * temp + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh - sequence_on_lattice[temp] * sequence_on_lattice[step];
                        }
                    }

                    //смотрим выигрыш
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * start_conformation + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh + sequence_on_lattice[start_conformation] * sequence_on_lattice[step];
                        }
                    }

                    new_E = E + hh;
                    new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];

                    //p1 = exp(-(new_E - E) * J - (new_H - current_H_counts) * h);

                    p1 = exp(-(-(new_E - E) * J - (new_H - current_H_counts) * h));
                    p_metropolis = std::min(1.0, p1);
                    q_rd = distribution(generator);

                    if (q_rd < p_metropolis) {
                        E = new_E;
                        current_H_counts = new_H;
                        sequence_on_lattice[temp] = 0; //делаю здесь, так как проще считать энергию(!!!)
                        sum_X = sum_X + lattice.x_coords[start_conformation] - lattice.x_coords[temp];
                        sum_Y = sum_Y + lattice.y_coords[start_conformation] - lattice.y_coords[temp];
                        sum_Z = sum_Z + lattice.z_coords[start_conformation] - lattice.z_coords[temp];
                        //sum_X = sum_X + lattice.x_coords[start_conformation];
                        //sum_Y = sum_Y + lattice.y_coords[start_conformation];

                        //корректируем информацию о направлениях
                        directions[end_conformation] = -1;
                        directions[start_conformation] = inverse_steps[step_on_lattice];

                    }
                    else {//отменяем изменения
                     //удаляем начало
                        del = start_conformation;
                        start_conformation = next_monomers[start_conformation];
                        previous_monomers[start_conformation] = -1;
                        next_monomers[del] = -1;
                        sequence_on_lattice[del] = 0;

                        //возвращаем конец
                        next_monomers[end_conformation] = temp;
                        previous_monomers[temp] = end_conformation;
                        end_conformation = temp;
                        sequence_on_lattice[end_conformation] = oldspin;

                        if (previous_monomers[temp] < 0) {
                            std::cout << "problem return " << std::endl;
                        }
                        if (temp < 0) {
                            std::cout << "problem return temp" << std::endl;
                        }
                    }
                }
                else {
                    //некуда идти
                }

            }
        }
        else {
            step_on_lattice = distribution1(generators1);

            new_point = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + step_on_lattice];

            //проверка, что проверенный узел занят спином
            if (sequence_on_lattice[new_point] != 0 && next_monomers[new_point] != -1 && new_point != previous_monomers[end_conformation])
            {
                Reconnect(step_on_lattice);
            }

        }

        if (i > steps_to_equilibrium) {
            save_calcs();
            radius();

        }

        if (i % (number_of_monomers * number_of_monomers) == 0)
        {
            radius_gyration();
            //std::cout << i << std::endl;
        }

        if (i > steps_to_equilibrium && i % 5000000000 == 0)
            //if ( i> steps_to_equilibrium && i%1000==0 )
        {
            std::string filename;
            std::ofstream out_result;

            filename = "Canonical_Ising_" + std::to_string(J) + "_" + std::to_string(h) + "_" + std::to_string(number_of_monomers) + ".txt";
            //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

            out_result.open(filename);
            //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";
            out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq ";
            out_result << "mean_e err_mean_e mean_e_sq err_mean_e_sq mean_e_fourth err_mean_e_fourth ";
            out_result << "mean_m err_mean_m mean_m_sq err_mean_m_sq mean_m_fourth err_mean_m_fourth " << std::endl;

            out_result << number_of_monomers << " " << J << " " << h << " ";
            out_result << dists.mean() << " " << dists.errorbar() << " " << gyration.mean() << " " << gyration.errorbar() << " ";

            out_result << energy.mean() << " " << energy.errorbar() << " ";
            out_result << energy_sq.mean() << " " << energy_sq.errorbar() << " ";
            out_result << energy_4.mean() << " " << energy_4.errorbar() << " ";

            out_result << magnetization.mean() << " " << magnetization.errorbar() << " ";
            out_result << magnetization_sq.mean() << " " << magnetization_sq.errorbar() << " ";
            out_result << magnetization_4.mean() << " " << magnetization_4.errorbar() << " ";
            out_result << i << std::endl; //neede i!!!!


            out_result.close();

            filename = "Counts_E_Ising_" + std::to_string(J) + "_" + std::to_string(number_of_monomers) + ".txt";
            //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

            out_result.open(filename);

            out_result << "N J h steps " << std::endl;
            out_result << number_of_monomers << " " << J << " " << h << " ";
            out_result << i << std::endl;
            for (auto counts : count_E)
            {
                out_result << counts.first << " " << counts.second << std::endl;
            }
            out_result.close();


            filename = "Counts_M_Ising_" + std::to_string(J) + "_" + std::to_string(number_of_monomers) + ".txt";
            //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

            out_result.open(filename);

            out_result << "N J h steps " << std::endl;
            out_result << number_of_monomers << " " << J << " " << h << " ";
            out_result << i << std::endl;
            for (auto counts : count_M)
            {
                out_result << counts.first << " " << counts.second << std::endl;
            }
            out_result.close();


        }


    }

}




void Protein::radius()
{
    long int point1x = end_conformation % lattice.lattice_side;
    long int point1z = end_conformation / (lattice.lattice_side * lattice.lattice_side);
    long int point1y = (end_conformation % (lattice.lattice_side * lattice.lattice_side)) / lattice.lattice_side;
    long int point1xs = start_conformation % lattice.lattice_side;
    long int point1zs = start_conformation / (lattice.lattice_side * lattice.lattice_side);
    long int point1ys = (start_conformation % (lattice.lattice_side * lattice.lattice_side)) / lattice.lattice_side;
    //расстояние на торе
    long int xdiff = abs(point1x - point1xs);
    if (xdiff > (lattice.lattice_side / 2))
        xdiff = lattice.lattice_side - xdiff;

    long int ydiff = abs(point1y - point1ys);
    if (ydiff > (lattice.lattice_side / 2))
        ydiff = lattice.lattice_side - ydiff;

    long int zdiff = abs(point1z - point1zs);
    if (zdiff > (lattice.lattice_side / 2))
        zdiff = lattice.lattice_side - zdiff;

    long int r = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
    dists << r;

}


void Protein::radius_gyration()
{


    long double r_g = 0;
    long int current = start_conformation;
    long double y = 0, x = 0, z = 0;
    long double point1x = 0, point1y = 0, point1z = 0;
    //long double point1x = 1.0*sum_X/number_of_monomers;
    //long double point1y = 1.0*sum_Y/number_of_monomers;
    //long double point1z = 1.0*sum_Z/number_of_monomers;
    long double xdiff, ydiff, zdiff;
    //point1x = start_conformation % lattice.lattice_side;
    //point1z = start_conformation / (lattice.lattice_side * lattice.lattice_side);
    //point1y = (start_conformation % (lattice.lattice_side * lattice.lattice_side)) / lattice.lattice_side;
    current = start_conformation;
    long double r;
    long int point1xs, point1ys, point1zs;

    long second_current;
    for (int e = 0; e < number_of_monomers; e++)
    {
        second_current = start_conformation;

        point1x = lattice.x_coords[current];
        point1y = lattice.y_coords[current];
        point1z = lattice.z_coords[current];

        for (int e1 = 0; e1 < number_of_monomers; e1++) {

            point1xs = lattice.x_coords[second_current];
            point1ys = lattice.y_coords[second_current];
            point1zs = lattice.z_coords[second_current];

            //расстояние на торе
            xdiff = abs(point1x - point1xs);
            if (xdiff > (lattice.lattice_side / 2))
                xdiff = lattice.lattice_side - xdiff;

            ydiff = abs(point1y - point1ys);
            if (ydiff > (lattice.lattice_side / 2))
                ydiff = lattice.lattice_side - ydiff;

            zdiff = abs(point1z - point1zs);
            if (zdiff > (lattice.lattice_side / 2))
                zdiff = lattice.lattice_side - zdiff;

            r = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;

            r_g = r_g + r;

            second_current = next_monomers[second_current];

        }
        current = next_monomers[current];

        //std::cout << current << " ";

    }
    //std::cout << std::endl;
    gyration << 0.5 * r_g / number_of_monomers / number_of_monomers;

}

void Protein::radius_gyration1()
{


    long double r_g = 0;
    long int current = start_conformation;
    long double y = 0, x = 0, z = 0;
    long double point1x = 0, point1y = 0, point1z;
    //long double point1x = 1.0*sum_X/number_of_monomers;
    //long double point1y = 1.0*sum_Y/number_of_monomers;
    //long double point1z = 1.0*sum_Z/number_of_monomers;
    long double xdiff, ydiff, zdiff;
    point1x = start_conformation % lattice.lattice_side;
    point1z = start_conformation / (lattice.lattice_side * lattice.lattice_side);
    point1y = (start_conformation % (lattice.lattice_side * lattice.lattice_side)) / lattice.lattice_side;

    for (int e = 0; e < number_of_monomers; e++)
    {

        long int point1xs = lattice.x_coords[current];
        long int point1ys = lattice.y_coords[current];
        long int point1zs = lattice.z_coords[current];


        //расстояние на торе
        xdiff = abs(point1x - point1xs);
        if (xdiff > (lattice.lattice_side / 2))
            xdiff = lattice.lattice_side - xdiff;

        ydiff = abs(point1y - point1ys);
        if (ydiff > (lattice.lattice_side / 2))
            ydiff = lattice.lattice_side - ydiff;

        zdiff = abs(point1z - point1zs);
        if (zdiff > (lattice.lattice_side / 2))
            zdiff = lattice.lattice_side - zdiff;

        x = x + xdiff;
        y = y + ydiff;
        z = z + zdiff;

        //r = xdiff *xdiff  + ydiff*ydiff + zdiff * zdiff;

        current = next_monomers[current];

        //std::cout << current << " ";

    }

    point1x = x / number_of_monomers;
    point1y = y / number_of_monomers;
    point1z = z / number_of_monomers;

    current = start_conformation;
    long double r;
    long int point1xs, point1ys, point1zs;
    for (int e = 0; e < number_of_monomers; e++)
    {

        point1xs = lattice.x_coords[current];
        point1ys = lattice.y_coords[current];
        point1zs = lattice.z_coords[current];

        //расстояние на торе
        xdiff = abs(point1x - point1xs);
        if (xdiff > (lattice.lattice_side / 2))
            xdiff = lattice.lattice_side - xdiff;

        ydiff = abs(point1y - point1ys);
        if (ydiff > (lattice.lattice_side / 2))
            ydiff = lattice.lattice_side - ydiff;

        zdiff = abs(point1z - point1zs);
        if (zdiff > (lattice.lattice_side / 2))
            zdiff = lattice.lattice_side - zdiff;

        r = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;

        r_g = r_g + r;


        current = next_monomers[current];

        //std::cout << current << " ";

    }
    //std::cout << std::endl;
    gyration << 1.0 * r_g / number_of_monomers;

}

void Protein::save_calcs()
{

    energy << 1.0 * (E) / number_of_monomers;
    energy_sq << 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers;
    energy_4 << 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers;

    magnetization << 1.0 * abs(current_H_counts) / number_of_monomers;
    magnetization_sq << 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers;
    magnetization_4 << 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers;


    count_E[E] = count_E[E] + 1;
    count_M[current_H_counts] = count_M[current_H_counts] + 1;

    //radius();

    /*
    energy << 1.0*(-E) ;
    energy_sq << 1.0*E*E;
    energy_4 << 1.0*E*E*E*E;
    magnetization << 1.0*current_H_counts/number_of_monomers;
    magnetization_sq << 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers;
    magnetization_4 << 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers* 1.0*current_H_counts/number_of_monomers;
*/
}