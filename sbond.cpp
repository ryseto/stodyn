/*
 *  sbond.cpp
 *  stodyn
 *
 *  Created by seto on 10/10/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "sbond.h"

SBond::SBond(int bond_number_, int d0, int d1, vec3d bond_pos){
	bond_number = bond_number_;
//	cerr << d0 << ' ' << d1 << endl;
	particle[0] = d0;
	particle[1] = d1;
	pos  = bond_pos;

}
