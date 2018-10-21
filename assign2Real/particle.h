#ifndef ASSIGN1CPP_PARTICLE_H
#define ASSIGN1CPP_PARTICLE_H

#include "color.h"

class Particle {
public:
    double x;
    double y;
    double mass;
    double vel_x, vel_y;
    double acc_x, acc_y;
    Color color;

    // Default constructor
    Particle() {
        x = 0.0;
        y = 0.0;
        mass = 0.0;
        vel_x = 0.0;
        vel_y = 0.0;
        acc_x = 0.0;
        acc_y = 0.0;
    }

    // No-forces constructor
    Particle(double x, double y, double mass): x(x), y(y), mass(mass) {
        vel_x = 0.0;
        vel_y = 0.0;
        acc_x = 0.0;
        acc_y = 0.0;
    };

    // Color constructor
    Particle(double x, double y, double mass, float red, float green, float blue): x(x), y(y), mass(mass) {
        vel_x = 0.0;
        vel_y = 0.0;
        acc_x = 0.0;
        acc_y = 0.0;
        color.red = red;
        color.green = green;
        color.blue = blue;
    };

    // Velocity constructor
    Particle(double x, double y, double mass, double vel_x, double vel_y): x(x), y(y), mass(mass), vel_x(vel_x), vel_y(vel_y) {
        acc_x = 0.0;
        acc_y = 0.0;
    };

    // Full constructor
    Particle(double x, double y, double mass, double vel_x, double vel_y, double acc_x, double acc_y)
    : x(x), y(y), mass(mass), vel_x(vel_x), vel_y(vel_y), acc_x(acc_x), acc_y(acc_y) {};

    float get_distance(const Particle& second_particle) const;
    void add_acceleration(float total_mass, float center_of_mass_x, float center_of_mass_y);
    void computeVelocity();
    void calculatePosition();
};

#endif //ASSIGN1CPP_PARTICLE_H
