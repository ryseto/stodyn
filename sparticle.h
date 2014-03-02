/*
 *  sparticle.h
 *  stodyn
 *
 *  Created by seto on 10/10/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef sparticle_h
#define sparticle_h 1
#include <vector>
#include "sbond.h"
#include "vec3d.h"
//#include "rigidClusterAnalysis.h"
using namespace std;
class SBond;
extern vector<SBond *> sbond;

class SParticle{	
private:
	int cn_size;
	vector<int> neighbor;
	vector<int> bonds;
	//	System *sy;
public:
	SParticle(int particle_number_, vec3d &p_){
		p = p_;
		particle_number = particle_number_;
	}	
	~SParticle(){}	
	vec3d p;
	int particle_number;
	void setbond(int bondnumber_){
		bonds.push_back(bondnumber_);		
	}
	void output() {
		cerr << p.x << ' ' << p.y << ' ' << p.z << endl;		
	}
	
	int num_bond(){ return bonds.size(); }
	int neighbor_particle(int i, int cutbond){
		if ( cutbond == bonds[i] )
			return -1;
		else
			return sbond[ bonds[i]]->other_side(particle_number); 
		
	}
};

#endif
