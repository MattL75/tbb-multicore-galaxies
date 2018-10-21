#include <Windows.h>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>
#include "TreeNode.h"
#include "vector.h"
#include "constants.h"
#include <GL/gl.h>
#include <GL/glu.h>

int TreeNode::getQuadrant(double x, double y) {
    if (x <= center_x && y >= center_y) {
        return 0;   // NW
    }
    else if (x >= center_x && y >= center_y) {
        return 1;   // NE
    }
    else if (x <= center_x && y <= center_y) {
        return 2;   // SW
    }
    else if (x >= center_x && y <= center_y) {
        return 3;   // SE
    }
    else {
        throw std::runtime_error("Can't determine quadrant.");
    }
}

void TreeNode::insertToNode(Particle *newParticle) {
    if (particleNum > 1) {
        // There are a lot of particles, existing one has already been split
        int quad = getQuadrant(newParticle->x, newParticle->y);
        if (!children[quad]) {
            children[quad] = new TreeNode(getQuadrantBox_X1(quad), getQuadrantBox_X2(quad), getQuadrantBox_Y1(quad), getQuadrantBox_Y2(quad));
        }
        children[quad]->insertToNode(newParticle);
    } else if (particleNum == 1) {
        // If you are not trying to insert, no need to split this into multiple quadrants so
        // here we are splitting the existing one since we are inserting.
        int quad = getQuadrant(particle->x, particle->y);
        if (!children[quad]) {
            children[quad] = new TreeNode(getQuadrantBox_X1(quad), getQuadrantBox_X2(quad), getQuadrantBox_Y1(quad), getQuadrantBox_Y2(quad));
        }
        children[quad]->insertToNode(particle);
        particle = nullptr;

        // Now for the new particle placement
        int newquad = getQuadrant(newParticle->x, newParticle->y);
        if (!children[newquad]) {
            children[newquad] = new TreeNode(getQuadrantBox_X1(newquad), getQuadrantBox_X2(newquad), getQuadrantBox_Y1(newquad), getQuadrantBox_Y2(newquad));
        }
        children[newquad]->insertToNode(newParticle);
    } else if (particleNum == 0) {
        particle = newParticle;
    }
    particleNum++;
}

void TreeNode::drawTree() {
    glColor3ub(255, 0, 0);
    glBegin(GL_LINE_STRIP);
    glVertex2d(box_x1, box_y1);
    glVertex2d(box_x1, box_y2);
    glVertex2d(box_x2, box_y2);
    glVertex2d(box_x2, box_y1);
    glVertex2d(box_x1, box_y1);
    glEnd();

    for (int i = 0; i < 4; ++i) {
        if (children[i])
            children[i]->drawTree();
    }
}

void TreeNode::computeMassDistribution() {
    if (particleNum == 1) {
        cm_x = particle->x;
        cm_y = particle->y;
        mass = particle->mass;
    } else {
        // Need to compute center of mass based on all of the
        // child quadrants, including centers of mass.
        mass = 0;
        cm_x = 0;
        cm_y = 0;

        for (int i = 0; i < 4; ++i) {
            if (children[i]) {
                children[i]->computeMassDistribution();
                mass += children[i]->mass;
                cm_x += children[i]->cm_x * children[i]->mass;
                cm_y += children[i]->cm_y * children[i]->mass;
            }
        }
        cm_x /= mass;
        cm_y /= mass;
    }
}

Vec2D TreeNode::calculateTreeForce(Particle &targetParticle) {
    Vec2D force;
    double r(0), k(0), d(0);

    if (particleNum == 1) {
        force = calcForce(targetParticle, *particle);
    } else {
        r = sqrt( (targetParticle.x - cm_x) * (targetParticle.x - cm_x) + (targetParticle.y - cm_y) * (targetParticle.y - cm_y));
        d = box_x2 - box_x1;
        if ((d/r) <= THETA) {
            k = (GRAV * mass * targetParticle.mass) / ((r + SOFTENER) * (r + SOFTENER) * (r + SOFTENER));
            force.x = k * (cm_x - targetParticle.x);
            force.y = k * (cm_y - targetParticle.y);
        } else {
            Vec2D temp;
            for (int i = 0; i < 4; ++i) {
                if (children[i]) {
                    temp = children[i]->calculateTreeForce(targetParticle);
                    force.x += temp.x;
                    force.y += temp.y;
                }
            }
        }
    }
    return force;
}

Vec2D TreeNode::calcForce(Particle &p1, Particle &p2) {
    Vec2D force;

    if (&p1==&p2) return force;
    double r = sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
    if (r > 0) {
        double k = GRAV * p1.mass * p2.mass / ((r + SOFTENER) * (r + SOFTENER) * (r + SOFTENER));
        force.x += k * (p2.x - p1.x);
        force.y += k * (p2.y - p1.y);
    } else {
        force.x = force.y = 0;
    }
    return force;
}

void TreeNode::resetTree() {
    delete(children[0]);
    delete(children[1]);
    delete(children[2]);
    delete(children[3]);
}

double TreeNode::getQuadrantBox_X1(int quadrant) {
    if (quadrant == 0 || quadrant == 2) {
        return box_x1;
    } else {
        return center_x;
    }
}

double TreeNode::getQuadrantBox_X2(int quadrant) {
    if (quadrant == 0 || quadrant == 2) {
        return center_x;
    } else {
        return box_x2;
    }
}

double TreeNode::getQuadrantBox_Y1(int quadrant) {
    if (quadrant == 0 || quadrant == 1) {
        return box_y1;
    } else {
        return center_y;
    }
}

double TreeNode::getQuadrantBox_Y2(int quadrant) {
    if (quadrant == 0 || quadrant == 1) {
        return center_y;
    } else {
        return box_y2;
    }
}