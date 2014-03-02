/*
 *  sbond.h
 *  stodyn
 *
 *  Created by seto on 10/10/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef sbond_h
#define sbond_h 1
#include "common.h"
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "contable.h"

using namespace std;
class SParticle;

class SBond {
	int bond_number;
	//Particle *p_particle0;
	//Particle *p_particle1;
private:
public:
	//SBond(){}
	SBond(int bond_number_, int d0, int d1, vec3d bond_pos);
	~SBond(){}

	int particle[2];
	vec3d pos;
	int other_side(int i){
		if (i == particle[0])
			return particle[1];
		else 
			return particle[0];
	}
	
	double force_normal;
	double force_sliding;
	double moment_bending;
	double moment_torsion;
	
};
#endif
