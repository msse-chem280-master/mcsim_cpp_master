#include <cmath>

#include "types.hpp"
#include "util.hpp"


double calculate_LJ(double r_ij)
{

    double inv_r_ij = 1.0 / r_ij;
    double r6_term = pow(inv_r_ij, 6.0);
    double r12_term = r6_term * r6_term;
    double pairwise_energy = 4.0 * (r12_term - r6_term);
    return pairwise_energy;
}

double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length = -1.0)
{

    double distance = 0;

    for (int i = 0; i < 3; i++) {
        double dim_dist = (coord1[i] - coord2[i]);

        if (box_length != -1.0) {
            dim_dist = dim_dist - box_length * round(dim_dist / box_length);
        }
        dim_dist = dim_dist * dim_dist;
        distance += dim_dist;
    }

    distance = sqrt(distance);

    return distance;
}

double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff)
{

    double total_energy = 0.0;

    for (size_t i = 0; i < coordinates.size(); i++) {
        for (size_t j = i + 1; j < coordinates.size(); j++) {

            double dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length);

            if (dist_ij < cutoff)
                total_energy += calculate_LJ(dist_ij);
        }
    }

    return total_energy;
}

double calculate_tail_correction(int num_particles, double box_length, double cutoff)
{

    double const1 = (8.0 * M_PI * num_particles * num_particles) / (3 * pow(box_length, 3));
    double const2 = (1.0 / 3.0) * pow((1.0 / cutoff), 9) - pow((1.0 / cutoff), 3);
    return const1 * const2;
}

double calculate_pair_energy(Coordinates coordinates, int i_particle, double box_length, double cutoff)
{

    double e_total = 0.0;
    int num_atoms = coordinates.size();

    for (int j = 0; j < num_atoms; j++) {
        if (j != i_particle) {
            double dist_ij = calculate_distance(coordinates[i_particle], coordinates[j], box_length);

            if (dist_ij < cutoff) {
                double interaction_energy = calculate_LJ(dist_ij);
                e_total += interaction_energy;
            }
        }
    }

    return e_total;
}

bool accept_or_reject(double delta_e, double beta)
{

    bool accept = false;
    if (delta_e <= 0) {
        accept = true;
    } else {
        double random_number = random_double(0, 1);
        double p_acc = exp(-beta * delta_e);
        if (random_number < p_acc) {
            accept = true;
        } else {
            accept = false;
        }
    }

    return accept;
}
