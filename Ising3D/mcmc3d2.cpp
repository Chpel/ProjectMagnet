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
    //                                                         
    long int x, y, z;
    div_t n;
    map_of_contacts_int.resize(lattice_side * lattice_side * lattice_side * ndim3());
    for (long int i = 0; i < lattice_side * lattice_side * lattice_side; i++) {
        map_of_contacts_int[ndim3() * i] = i + 1;
        map_of_contacts_int[ndim3() * i + 1] = i - 1;
        map_of_contacts_int[ndim3() * i + 2] = i + lattice_side;
        map_of_contacts_int[ndim3() * i + 3] = i - lattice_side;
        //                                             :
        map_of_contacts_int[ndim3() * i + 4] = i + lattice_side * lattice_side;
        map_of_contacts_int[ndim3() * i + 5] = i - lattice_side * lattice_side;
        //                 z: i = z * l_s^2 + y * l_s + x
        //       :
        n = div(i, lattice_side * lattice_side);
        z = i / (lattice_side * lattice_side);
        //n = div( (long)n.rem, lattice_side);
        x = n.rem % lattice_side;
        y = n.rem / lattice_side;
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
            //                               :
            if (z == 0) {
                //map_of_contacts_int[6 * i + 5] = x + lattice_side * y + lattice_side * lattice_side * (lattice_side - 1);
                map_of_contacts_int[6 * i + 5] = i + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                //map_of_contacts_int[6 * i + 4] = x + lattice_side * y;
                map_of_contacts_int[6 * i + 4] = i - lattice_side * lattice_side * (lattice_side - 1);
            }
        }
    }
}

void Lattice::create_lattice(long int max_seq_size) {

    lattice_side = max_seq_size;
    long int x, y, z;
    div_t n;
    map_of_contacts_int.resize(lattice_side * lattice_side * lattice_side * ndim3());
    x_coords.resize(lattice_side * lattice_side * lattice_side * ndim3());
    y_coords.resize(lattice_side * lattice_side * lattice_side * ndim3());
    z_coords.resize(lattice_side * lattice_side * lattice_side * ndim3());
    for (long int i = 0; i < lattice_side * lattice_side * lattice_side; i++) {
        map_of_contacts_int[ndim3() * i] = i + 1;
        map_of_contacts_int[ndim3() * i + 1] = i - 1;
        map_of_contacts_int[ndim3() * i + 2] = i + lattice_side;
        map_of_contacts_int[ndim3() * i + 3] = i - lattice_side;
        //                                             :
        map_of_contacts_int[ndim3() * i + 4] = i + lattice_side * lattice_side;
        map_of_contacts_int[ndim3() * i + 5] = i - lattice_side * lattice_side;
        //                 z: i = z * l_s^2 + y * l_s + x
        //       :
        n = div(i, lattice_side * lattice_side);
        z = n.quot;
        n = div(n.rem, lattice_side);
        x = n.rem;
        y = n.quot;
        //x = n.rem % lattice_side; //n.rem;
        //y = n.rem / lattice_side; //n.quot;
        for (int j = 0; j < ndim3(); j++) {
            if (x == 0) {
                map_of_contacts_int[6 * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[6 * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                //map_of_contacts_int[6 * i + 3] = lattice_side * (lattice_side - 1) + x;
                map_of_contacts_int[6 * i + 3] = lattice_side * (lattice_side - 1) + i;
            }
            if (y == (lattice_side - 1)) {
                //map_of_contacts_int[6 * i + 2] = x;
                map_of_contacts_int[6 * i + 2] = i - lattice_side * (lattice_side - 1);
            }
            //                               :
            if (z == 0) {
                //map_of_contacts_int[6 * i + 5] = x + lattice_side * y + lattice_side * lattice_side * (lattice_side - 1);
                map_of_contacts_int[6 * i + 5] = i + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                //map_of_contacts_int[6 * i + 4] = x + lattice_side * y;
                map_of_contacts_int[6 * i + 4] = i - lattice_side * lattice_side * (lattice_side - 1);
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
    lattice.create_lattice(n + 5); //                ,    5       ,                  
    // 0 -                               
    sequence_on_lattice.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, 0); //                            
    next_monomers.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //     -                        
    previous_monomers.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //     -                          

    directions.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //               {0,1,2,3}

    ordered_coords.resize(lattice.lattice_side * lattice.lattice_side * lattice.lattice_side, -1); //                                       

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
    sequence_on_lattice[end_conformation] = 1; //                            

    next_monomers[0] = 1;
    previous_monomers[n - 1] = n - 2;

    current_H_counts = n;
    E = -(n - 1);


    //                        -                
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
        //            ,                  
    };

    //int j = directions[previous_monomers[end_conformation]]; //                            

    long int step = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + j];
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

void Protein::calc_bulk()
{
    bulk6_now = 0, bulk5_now = 0, bulk4_now = 0, bulk3_now = 0, bulk2_now = 0;

    long int current = start_conformation;
    long int step;
    int k = 0;
    for (int e = 0; e < number_of_monomers; e++)
    {
        k = 0;
        for (int j = 0; j < lattice.ndim3(); j++) {
            step = lattice.map_of_contacts_int[lattice.ndim3() * current + j];
            if (sequence_on_lattice[step] != 0) {
                k += 1;

            }
        }

        if (k == 2) {
            bulk2_now += 1;
        }
        if (k == 3) {
            bulk3_now += 1;
        }
        if (k == 4) {
            bulk4_now += 1;
        }
        if (k == 5) {
            bulk5_now += 1;
        }
        if (k == 6) {
            bulk6_now += 1;
        }

        current = next_monomers[current];

    }


    //std::cout << current << " ";

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
        assert(c >= 0 && c < lattice.lattice_side* lattice.lattice_side* lattice.lattice_side);
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
    //[i,j]: i-       , j-            (              )
    int reflect_directions[4][4] =
    { {3, 2, 0, 1 }, //90
     {1,0,3,2 }, //180
     {2,3, 1, 0 }, //270
    { 0, 1, 2, 3}
    };

    int inverse_steps[6] = { 1,0,3,2,5,4 }; //                                     "                        "
    //double step_rd; //                  :                      
    double q_rd, p1, p_metropolis; //                             
    int rand_path; // = distribution1(generators1); //                      : 0 -                            
    double typeOfUpdate; //0 -        ; 1 -          
    int step;
    int step_on_lattice;//                       
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

            if (rand_path == 0) {//                           

                step_on_lattice = distribution1(generators1);
                new_point = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[start_conformation];

                if (sequence_on_lattice[new_point] == 0) { //        ,                         

                    //             

                    //                 
                    next_monomers[end_conformation] = new_point;
                    sequence_on_lattice[new_point] = 2 * distribution2(generators3) - 1; //            //!!!!!!!!!
                    previous_monomers[new_point] = end_conformation;
                    end_conformation = new_point;

                    //              
                    temp = start_conformation;
                    start_conformation = next_monomers[start_conformation];
                    next_monomers[temp] = -1;
                    previous_monomers[start_conformation] = -1;

                    //              
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * temp + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh - sequence_on_lattice[temp] * sequence_on_lattice[step];
                            //std::cout << sequence_on_lattice[temp] * sequence_on_lattice[step] <<std::endl;
                        }
                    }

                    //               
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh + sequence_on_lattice[end_conformation] * sequence_on_lattice[step];
                            //std::cout << sequence_on_lattice[end_conformation] * sequence_on_lattice[step] <<std::endl;
                        }
                    }

                    new_E = E + hh;
                    //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[start_conformation];
                    new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];


                    p1 = exp(-(-(new_E - E) * J - (new_H - current_H_counts) * h));


                    //std::cout << p1 << std::endl;

                    p_metropolis = std::min(1.0, p1);
                    q_rd = distribution(generator);
                    if (q_rd < p_metropolis) {
                        //std::cout << E << " " << new_E << " " << hh <<std::endl;
                        E = new_E;
                        current_H_counts = new_H;
                        sequence_on_lattice[temp] = 0; //           ,                              (!!!)
                        sum_X = sum_X + lattice.x_coords[end_conformation] - lattice.x_coords[temp];
                        sum_Y = sum_Y + lattice.y_coords[end_conformation] - lattice.y_coords[temp];
                        sum_Z = sum_Z + lattice.z_coords[end_conformation] - lattice.z_coords[temp];
                        //                                      
                        directions[temp] = -1;
                        directions[previous_monomers[end_conformation]] = step_on_lattice;

                    }
                    else {//                  
                        //             
                        del = end_conformation;
                        end_conformation = previous_monomers[end_conformation];
                        next_monomers[end_conformation] = -1;
                        previous_monomers[del] = -1;
                        sequence_on_lattice[del] = 0;

                        //                 
                        previous_monomers[start_conformation] = temp;
                        next_monomers[temp] = start_conformation;
                        start_conformation = temp;
                        sequence_on_lattice[start_conformation] = oldspin;

                    }

                }
                else {//         ,                
                    //continue;
                }

            }
            else {//                           

                step_on_lattice = distribution1(generators1);
                new_point = lattice.map_of_contacts_int[lattice.ndim3() * start_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[end_conformation];

                if (sequence_on_lattice[new_point] == 0) { //        ,                         

                    //             
                    //                  
                    previous_monomers[start_conformation] = new_point;
                    sequence_on_lattice[new_point] = 2 * distribution2(generators3) - 1; //            //!!!!!!!!!
                    next_monomers[new_point] = start_conformation;
                    start_conformation = new_point;

                    //             
                    temp = end_conformation;
                    end_conformation = previous_monomers[end_conformation];
                    if (previous_monomers[end_conformation] < 0) {
                        std::cout << "problem update " << std::endl;
                    }
                    previous_monomers[temp] = -1;
                    next_monomers[end_conformation] = -1;

                    //              
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * temp + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh - sequence_on_lattice[temp] * sequence_on_lattice[step];
                            //std::cout << sequence_on_lattice[temp] * sequence_on_lattice[step] <<std::endl;
                        }
                    }

                    //               
                    for (int j = 0; j < lattice.ndim3(); j++) {
                        step = lattice.map_of_contacts_int[lattice.ndim3() * start_conformation + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh + sequence_on_lattice[start_conformation] * sequence_on_lattice[step];
                            //std::cout << sequence_on_lattice[start_conformation] * sequence_on_lattice[step] <<std::endl;
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
                        sequence_on_lattice[temp] = 0; //           ,                              (!!!)
                        sum_X = sum_X + lattice.x_coords[start_conformation] - lattice.x_coords[temp];
                        sum_Y = sum_Y + lattice.y_coords[start_conformation] - lattice.y_coords[temp];
                        sum_Z = sum_Z + lattice.z_coords[start_conformation] - lattice.z_coords[temp];
                        //sum_X = sum_X + lattice.x_coords[start_conformation];
                        //sum_Y = sum_Y + lattice.y_coords[start_conformation];

                        //                                      
                        directions[end_conformation] = -1;
                        directions[start_conformation] = inverse_steps[step_on_lattice];

                    }
                    else {//                  
                     //              
                        del = start_conformation;
                        start_conformation = next_monomers[start_conformation];
                        previous_monomers[start_conformation] = -1;
                        next_monomers[del] = -1;
                        sequence_on_lattice[del] = 0;

                        //                
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
                    //           
                }

            }
        }
        else {
            step_on_lattice = distribution1(generators1);

            new_point = lattice.map_of_contacts_int[lattice.ndim3() * end_conformation + step_on_lattice];

            //        ,                                  
            if (sequence_on_lattice[new_point] != 0 && next_monomers[new_point] != -1 && new_point != previous_monomers[end_conformation])
            {
                Reconnect(step_on_lattice);
            }

        }


        /*long int c = start_conformation;
        std::cout<<start_conformation << ", " ;
        for (int e = 0; e < number_of_monomers; e++)
        {
            std::cout<< next_monomers[c] << ", " ;
            c = next_monomers[c];
        }
        std::cout<<std::endl;
        std::cout << E << std::endl; */
        //if (i>100) break;

        if (i > steps_to_equilibrium && i % 10000000 == 0) {
            save_calcs();
            //radius();

            calc_bulk();
            bulk2 << 1.0 * bulk2_now / number_of_monomers;
            bulk3 << 1.0 * bulk3_now / number_of_monomers;
            bulk4 << 1.0 * bulk4_now / number_of_monomers;
            bulk5 << 1.0 * bulk5_now / number_of_monomers;
            bulk6 << 1.0 * bulk6_now / number_of_monomers;

        }

        if (i % (number_of_monomers * number_of_monomers) == 0)
        {
            radius_gyration();
            //std::cout << i << std::endl;
        }

        if (i > steps_to_equilibrium && i % 100000000 == 0)
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


            filename = "Geometry_Ising_" + std::to_string(J) + "_" + std::to_string(h) + "_" + std::to_string(number_of_monomers) + ".txt";
            out_result.open(filename);
            
            out_result << "N J h b2 err_b2 b3 err_b3 b4 err_b4 b5 err_b5 b6 err_b6";
            out_result << number_of_monomers << " " << J << " " << h << " ";

            out_result << bulk2.mean() << " " << bulk2.errorbar() << " ";
            out_result << bulk3.mean() << " " << bulk3.errorbar() << " ";
            out_result << bulk4.mean() << " " << bulk4.errorbar() << " ";
            out_result << bulk5.mean() << " " << bulk5.errorbar() << " ";
            out_result << bulk6.mean() << " " << bulk6.errorbar() << " ";

            out_result << i << std::endl; //neede i!!!!


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
    //                  
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

            //                  
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


        //                  
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

        //                  
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