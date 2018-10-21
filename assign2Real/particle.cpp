#include <Windows.h>
#include "particle.h"
#include <cmath>

float Particle::get_distance(const Particle& second_particle) const {
    // Get distances
    float dx = x - second_particle.x;
    float dy = y - second_particle.y;

    // Square of distances
    float distance_square = dx * dx + dy * dy;

    // Keep a minimum square of distance
    if (distance_square < 1)
        distance_square = 1;

    return sqrt(distance_square);
}

void Particle::computeVelocity() {
    vel_x += acc_x;
    vel_y += acc_y;
}

void Particle::calculatePosition() {
    x += vel_x;
    y += vel_y;
}

void Particle::add_acceleration(float total_mass, float center_of_mass_x, float center_of_mass_y) {

}