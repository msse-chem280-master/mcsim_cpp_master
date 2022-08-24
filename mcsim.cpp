#include <iostream>
#include <string>

#include "util.hpp"
#include "calculate.hpp"


Coordinates run_simulation(Coordinates coordinates, const double box_length, const double cutoff,
    const double reduced_temperature, const int num_steps, const double max_displacement, const int freq = 1000)
{

    std::vector<int> steps;
    std::vector<double> energies;
    std::vector<Coordinates> all_coordinates;

    double beta = 1.0 / reduced_temperature;

    int num_particles = coordinates.size();

    double total_energy = calculate_total_energy(coordinates, box_length, cutoff);
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);

    for (int i = 0; i < num_steps; i++) {

        int random_particle = random_integer(0, num_particles);

        double current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;

        double proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double delta_energy = proposed_energy - current_energy;

        bool accept = accept_or_reject(delta_energy, beta);

        if (accept) {
            total_energy += delta_energy;
        } else {
            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
        }

        if (i % freq == 0) {
            std::cout << i << ' ' << total_energy / num_particles << std::endl;
            steps.push_back(i);
            energies.push_back(total_energy / num_particles);
            all_coordinates.push_back(coordinates);
        }
    }

    return coordinates;
}

int main(int argc, char ** argv)
{
    if(argc < 2)
    {
        std::cout << "Missing number of steps as the argument" << std::endl;
        return 1;
    }

    double reduced_temperature = 1.5;
    double max_displacement = 0.1; 
    double cutoff = 3.0;

    int num_steps = std::atoi(argv[1]);

    // Initialize the global random number generator
    initialize_rng();

    std::pair<Coordinates, double> xyz_info = read_xyz("lj_sample_configurations/lj_sample_config_periodic1.txt");
    Coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;

    Coordinates newcoords = run_simulation(coords, box_length, cutoff, reduced_temperature, num_steps, max_displacement);

    return 0;
}
