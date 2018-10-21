#include <windows.h>
#ifndef ASSIGN1CPP_COMPUTETBB_H
#define ASSIGN1CPP_COMPUTETBB_H
#include <numeric>
#include "tbb/tbb.h"
#include "particle.h"
#include "TreeNode.h"

using namespace tbb;

class ComputeTBB {
	Particle *particles[10002];
	TreeNode *root;
public:
	void operator()(const blocked_range<size_t> &r) const {
		for (size_t i = r.begin(); i != r.end(); ++i) {
			Vec2D force = root->calculateTreeForce(*particles[i]);
			particles[i]->acc_x = force.x / particles[i]->mass;
			particles[i]->acc_y = force.y / particles[i]->mass;
			particles[i]->computeVelocity();
			particles[i]->calculatePosition();
		}
	}

	ComputeTBB(Particle a[10002], TreeNode *root) : root(root) {
		*particles = a;
	}
};

#endif //ASSIGN1CPP_COMPUTETBB_H