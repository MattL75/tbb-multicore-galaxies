#include <iostream>
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <tbb/tbb.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GLFW/glfw3.h>
#include "TreeNode.h"
#include "constants.h"
#include <random>
#include <chrono>
#include <time.h>
#include <math.h>
#include <cmath>
#include <mutex>
#include "mingw.thread.h"
#include "ComputeTBB.h"

#define PI 3.14159265

GLint windowWidth = 800;                                    // Width of our window
GLint windowHeight = 800;                                   // Height of our window
const int PARTICLE_NUMBER = 10000;                          // Number of particles
const int GALAXY_NUMBER = 2;                                // Number of galaxies
const int TOTAL_PARTICLES = PARTICLE_NUMBER + GALAXY_NUMBER;
const double STAR_MASS = 10;                                // Mass of a star
const int GALAXY_MASS = 1000000;                            // Mass of a galaxy
const bool debug = false;                                   // Debug mode boolean
TreeNode *root;                                             // Reference to the root tree node
Particle *particles[TOTAL_PARTICLES];                       // Array of particle pointers
const double NUM_SECONDS = 0.033333333;                     // Number of seconds per iteration
const int NUM_THREADS = 4;                                  // Total number of threads
std::thread myThreads[NUM_THREADS];                         // Array of threads
std::mutex m;

// Method used to set viewport on window resize
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    glViewport(0, 0, width, height);
}

// Method used for key inputs
void processInput(GLFWwindow *window) {

    // Exit window if escape key is pressed
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// Method to build the barnes-hut quad tree
void buildTree() {
    for (int i = 0; i < PARTICLE_NUMBER; ++i) {
        // Ignore particles that are out of bounds
        if (!(particles[i]->x > windowWidth || particles[i]->x < 0 || particles[i]->y > windowHeight || particles[i]->y < 0)) {

            // Insert to root, recursive
            root->insertToNode(particles[i]);
        }
    }
}

void ParallelTBB() {
	// tbb::parallel_for(tbb::blocked_range<size_t>(0, TOTAL_PARTICLES), ComputeTBB(*particles, root));
	tbb::parallel_for(tbb::blocked_range<size_t>(0, TOTAL_PARTICLES - 1),
		[&](tbb::blocked_range<size_t>& r)
	{
		for (int i = r.begin(); i < r.end(); ++i)
		{
			Vec2D force = root->calculateTreeForce(*particles[i]);
			particles[i]->acc_x = force.x / particles[i]->mass;
			particles[i]->acc_y = force.y / particles[i]->mass;
			particles[i]->computeVelocity();
			particles[i]->calculatePosition();
		}
	});
}

// Thread method to calculate particle forces
void forcesParallelSetup(int threadid) {

    // Fraction of particles per thread
    int fraction = (PARTICLE_NUMBER) / NUM_THREADS;

    // Loop start
    int start = threadid * fraction;

    // Loop end
    int end = (threadid >= NUM_THREADS ? TOTAL_PARTICLES : fraction * (threadid + 1));

    // Loop to calculate positions
    for (int j = start; j < end; ++j) {
        Vec2D force = root->calculateTreeForce(*particles[j]);
        particles[j]->acc_x = force.x / particles[j]->mass;
        particles[j]->acc_y = force.y / particles[j]->mass;
        particles[j]->computeVelocity();
        particles[j]->calculatePosition();
    }
}

// Method called to calculate forces parallel
void calcForcesParallel() {
    int argument = 0;   // Necessary for some reason - bug in thread library?
    for (int i = 0; i < NUM_THREADS; ++i) {
        myThreads[i] = std::thread(forcesParallelSetup, argument++);
    }
    int count = 0;      // Also necessary
    for (int i = 0; i < NUM_THREADS; ++i) {
        myThreads[count++].join();
    }
}

// Sequential method to calculate forces
void calcForces() {
    for (int i = 0; i < TOTAL_PARTICLES; ++i) {
        Vec2D force = root->calculateTreeForce(*particles[i]);
        particles[i]->acc_x = force.x / particles[i]->mass;
        particles[i]->acc_y = force.y / particles[i]->mass;
        particles[i]->computeVelocity();
        particles[i]->calculatePosition();
    }
}

// Debug method that draws the quadrants
void drawTree() {
    root->drawTree();
}

// Initial state for random particles
void initRandomParticles() {

    // Random number generation code
    std::default_random_engine generator;
    std::uniform_real_distribution<double> width(0, windowWidth);
    std::uniform_real_distribution<double> height(0, windowHeight);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    // Generate particles
    for (int i = 0; i < PARTICLE_NUMBER - 1; ++i) {
        particles[i] = new Particle(width(generator), height(generator), STAR_MASS);
    }

    // Generate galaxy
    particles[PARTICLE_NUMBER - 1] = new Particle(400, 400, GALAXY_MASS);
}

// Method to get necessary orbital velocity around a particle
void GetOrbitalVelocity(Particle *p1, Particle *p2) {

    // Setup distance constants
    double r[2], dist;
    r[0] = p1->x - p2->x;
    r[1] = p1->y - p2->y;

    // Calculate distance
    dist = sqrt(r[0] * r[0] + r[1] * r[1]);

    // Calculate velocity
    double v = sqrt(GRAV * p1->mass / dist) * 0.9;

    // Apply velocity
    p2->vel_x = (r[1] / dist) * v;
    p2->vel_y = (-r[0] / dist) * v;
}

// Method to generate one galaxy in the middle
void spawnGalaxy(double radius) {

    // Set galaxy particles
    particles[0] = new Particle(windowWidth / 2, windowHeight / 2, GALAXY_MASS);

    // Prepare random number generation
    std::default_random_engine generator;
    std::uniform_real_distribution<double> length(0, radius);
    std::uniform_real_distribution<double> angle(0, 2 * PI);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    // First galaxy
    for (int i = 1; i < PARTICLE_NUMBER; ++i) {

        // Get random numbers for angle and length
        double angler = angle(generator);
        double lengther = length(generator);

        // Set new particle
        particles[i] = new Particle(windowWidth / 2 + (cos(angler) * lengther),
                                    windowHeight / 2 + (sin(angler) * lengther), STAR_MASS);

        // Generate orbital velocity around galaxy
        GetOrbitalVelocity(particles[0], particles[i]);
    }
}

// Demo initial state for a collision between galaxies
void spawnParallelCollision() {
    double g1x = 350;
    double g1y = 350;
    double g2x = 650;
    double g2y = 650;
    int g1r = 100;
    int g2r = 100;

    // Set galaxy particles (Maybe color black for black hole?)
    // Initial velocity so they circle each other, not collide
    particles[0] = new Particle(g1x, g1y, GALAXY_MASS, 0.5, -0.5);
    particles[1] = new Particle(g2x, g2y, GALAXY_MASS, -0.5, 0.5);

    // Prepare random number generation
    std::default_random_engine generator;
    std::uniform_real_distribution<double> length(0, g1r);
    std::uniform_real_distribution<double> angle(0, 2 * PI);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    // First galaxy
    for (int i = 2; i < TOTAL_PARTICLES / 2; ++i) {

        // Get random numbers for angle and length
        double angler = angle(generator);
        double lengther = length(generator);

        // Set new particle
        particles[i] = new Particle(g1x + (cos(angler) * lengther), g1y + (sin(angler) * lengther), STAR_MASS, 0.75f, 1,
                                    1);

        // Generate orbital velocity around galaxy
        GetOrbitalVelocity(particles[0], particles[i]);
    }

    // Second galaxy
    std::uniform_real_distribution<double> length2(0, g2r);
    for (int i = TOTAL_PARTICLES / 2; i < TOTAL_PARTICLES; ++i) {

        // Get random numbers for angle and length
        double angler = angle(generator);
        double lengther = length2(generator);

        // Set new particle
        particles[i] = new Particle(g2x + (cos(angler) * lengther), g2y + (sin(angler) * lengther), STAR_MASS, 1.0f,
                                    0.75f, 0.75f);

        // Generate orbital velocity around galaxy
        GetOrbitalVelocity(particles[1], particles[i]);
    }
}

// Demo initial state for a collision between galaxies
void spawnCollision(double g1x, double g1y, double g2x, double g2y, int g1r, int g2r) {

    // Set galaxy particles (Maybe color black for black hole?)
    // Initial velocity so they circle each other, not collide
    particles[0] = new Particle(g1x, g1y, GALAXY_MASS, 0.5, -0.5);
    particles[1] = new Particle(g2x, g2y, GALAXY_MASS, -0.5, 0.5);
//    particles[0] = new Particle(g1x, g1y, GALAXY_MASS, 1.00, 0.00);
//    particles[1] = new Particle(g2x, g2y, GALAXY_MASS, -1.00, 0.00);

    // Prepare random number generation
    std::default_random_engine generator;
    std::uniform_real_distribution<double> length(0, g1r);
    std::uniform_real_distribution<double> angle(0, 2 * PI);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    // First galaxy
    for (int i = 2; i < TOTAL_PARTICLES / 2; ++i) {

        // Get random numbers for angle and length
        double angler = angle(generator);
        double lengther = length(generator);

        // Set new particle
        particles[i] = new Particle(g1x + (cos(angler) * lengther), g1y + (sin(angler) * lengther), STAR_MASS, 0.75f, 1,
                                    1);

        // Generate orbital velocity around galaxy
        GetOrbitalVelocity(particles[0], particles[i]);
    }

    // Second galaxy
    std::uniform_real_distribution<double> length2(0, g2r);
    for (int i = TOTAL_PARTICLES / 2; i < TOTAL_PARTICLES; ++i) {

        // Get random numbers for angle and length
        double angler = angle(generator);
        double lengther = length2(generator);

        // Set new particle
        particles[i] = new Particle(g2x + (cos(angler) * lengther), g2y + (sin(angler) * lengther), STAR_MASS, 1.0f,
                                    0.75f, 0.75f);

        // Generate orbital velocity around galaxy
        GetOrbitalVelocity(particles[1], particles[i]);
    }
}

// Main draw method
void updateAndDraw() {
    // Setup new positions
    delete (root);
    root = new TreeNode(0, windowWidth, windowHeight, 0, nullptr);
    buildTree();
    root->computeMassDistribution();
    // calcForces();
    // calcForcesParallel();

	tbb::task_group g;
	g.run([&] {ParallelTBB(); });
	g.wait();

    // Viewport/buffer clear
    glViewport(0, 0, windowWidth, windowHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Setup matrix and projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, windowWidth, 0.0, windowHeight, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);

    // Debug mode to draw quadrants
    if (debug) {
        drawTree();
    }

    // Actual drawing
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
	for (int i = 0; i < TOTAL_PARTICLES; ++i) {

		// Check if particle is outside the bounding box
		if (!(particles[i]->x > windowWidth || particles[i]->x < 0 || particles[i]->y > windowHeight ||
			particles[i]->y < 0)) {

			// Set color as specified in Particle object
			glColor3f(particles[i]->color.red, particles[i]->color.green, particles[i]->color.blue);

			// Draw particle
			glVertex3d(particles[i]->x, particles[i]->y, 0.0);
		}
	}
	//    tbb::task_group r;
	//    r.run([&] {
	//        for (int i = 0; i < TOTAL_PARTICLES; ++i) {
	//
	//            // Check if particle is outside the bounding box
	//            if (!(particles[i]->x > windowWidth || particles[i]->x < 0 || particles[i]->y > windowHeight ||
	//                  particles[i]->y < 0)) {
	//
	//                // Set color as specified in Particle object
	//                glColor3f(particles[i]->color.red, particles[i]->color.green, particles[i]->color.blue);
	//
	//                // Draw particle
	//                glVertex3d(particles[i]->x, particles[i]->y, 0.0);
	//            }
	//        }
	//    });
	//    r.wait();
    glEnd();
    glFlush();
}

int main() {

    // Initialise glfw and opengl
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

    // Create window
    GLFWwindow *window = glfwCreateWindow(windowWidth, windowHeight, "Galaxy Simulation", NULL, NULL);
    if (window == NULL) {
        std::cout << "Failed to create GLFW window." << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Initialize OpenGL
    glViewport(0, 0, windowWidth, windowHeight);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // Generate initial state sequentially
    // spawnCollision(350, 350, 650, 650, 100, 100);

    // Prepare the initial state in parallel
	spawnParallelCollision();
    //std::thread initialState = std::thread(spawnParallelCollision);
    //initialState.join();

    // Timer base set
    double time_counter = 0;

    // Keep track of time
    clock_t this_time = clock();
    clock_t last_time = this_time;

    // Main drawing loop
    while (!glfwWindowShouldClose(window)) {

        // Manage time
        this_time = clock();
        time_counter += (double) (this_time - last_time);
        last_time = this_time;

        // Only go through if the set amount of time has passed
        if (time_counter > (double) (NUM_SECONDS * CLOCKS_PER_SEC)) {
            time_counter -= (double) (NUM_SECONDS * CLOCKS_PER_SEC);
            processInput(window);
            updateAndDraw();

            glfwSwapBuffers(window);
            glfwPollEvents();
        }
    }

    // Terminate procedure
    glfwTerminate();
    return 0;
}