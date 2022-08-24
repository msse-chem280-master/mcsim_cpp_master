#pragma once

#include "types.hpp"

double calculate_LJ(double r_ij);

double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length = -1.0);

double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff);

double calculate_tail_correction(int num_particles, double box_length, double cutoff);

double calculate_pair_energy(Coordinates coordinates, int i_particle, double box_length, double cutoff);

bool accept_or_reject(double delta_e, double beta);
