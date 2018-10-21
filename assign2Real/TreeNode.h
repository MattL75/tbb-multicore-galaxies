#ifndef ASSIGN1CPP_TREENODE_H
#define ASSIGN1CPP_TREENODE_H

#include "particle.h"
#include "vector.h"

// Class that acts as a quad tree
class TreeNode {
public:
    double box_x1, box_y1, box_x2, box_y2;  // Boundaries of the node quadrant. 1 is top left, 2 is bottom right;
    double center_x, center_y;              // Center of the quadrant
    int particleNum = 0;                    // Number of particles in this node
    Particle *particle;                      // Only carry one particle per node
    TreeNode *children[4];                  // Children nodes, max 4
    double cm_x, cm_y, mass;                // Center of mass positions

    // Default constructor
    TreeNode() {
        box_x1 = 0.0;
        box_x2 = 0.0;
        box_y1 = 0.0;
        box_y2 = 0.0;
        center_x = 0.0;
        center_y = 0.0;
    }

    // No-particle constructor
    TreeNode(double box_x1, double box_x2, double box_y1, double box_y2)
    : box_x1(box_x1), box_x2(box_x2), box_y2(box_y2), box_y1(box_y1) {
        center_x = (box_x1 + box_x2) / 2;   // MAYBE WRONG
        center_y = (box_y1 + box_y2) / 2;
        children[0] = nullptr;
        children[1] = nullptr;
        children[2] = nullptr;
        children[3] = nullptr;
    }

    // Particle constructor
    TreeNode(double box_x1, double box_x2, double box_y1, double box_y2, Particle* particle)
    : box_x1(box_x1), box_x2(box_x2), box_y2(box_y2), box_y1(box_y1), particle(particle) {
        center_x = (box_x1 + box_x2) / 2;
        center_y = (box_y1 + box_y2) / 2;
        children[0] = nullptr;
        children[1] = nullptr;
        children[2] = nullptr;
        children[3] = nullptr;
    }

    ~TreeNode() {
        delete(children[0]);
        delete(children[1]);
        delete(children[2]);
        delete(children[3]);
    }

    // Method to get quadrant of a particle
    int getQuadrant(double x, double y);

    // Use to get the top-left x-value of quadrant
    double getQuadrantBox_X1(int quadrant);

    // Use to get the top-left y-value of quadrant
    double getQuadrantBox_Y1(int quadrant);

    // Use to get the bottom-right x-value of quadrant
    double getQuadrantBox_X2(int quadrant);

    // Use to get the bottom-right y-value of quadrant
    double getQuadrantBox_Y2(int quadrant);

    // Method to compute the mass distribution
    void computeMassDistribution();

    // Method to calculate total force on a particle
    Vec2D calculateTreeForce(Particle &targetParticle);

    // Method to calculate force between two particles
    Vec2D calcForce(Particle &p1, Particle &p2);

    // Inserting to node
    void insertToNode(Particle *particle);

    // Method to reset the tree
    void resetTree();

    // Method to draw the tree
    void drawTree();
};

#endif //ASSIGN1CPP_TREENODE_H
