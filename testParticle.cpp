//
//  testParticle.cpp
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "testParticle.h"

TestParticle::TestParticle(int particle_number_, 
                           TestSystem &tsys_,
                           SDsystem &sd_sys_   ){
	tsys = &tsys_;
    sd_sys = &sd_sys_;
    particle_number = particle_number_;
	cn_size = 0;
    
	setVelocityZero();
    orientation.set(1.,0.,0.,0.);
    resetForce();
}

void TestParticle::generateBond(){
    int n_neighbor = neighbor.size();
	for (int i=0; i < n_neighbor; i++){
		if (sq_dist(p, tsys->particle[neighbor[i]]->p) <= tsys->sq_dist_generate ){
//			dem->bond.push_back(new Bond ( neighbor[i], particle_number, *dem));            
			neighbor[i] = neighbor.back();
			neighbor.pop_back();
            n_neighbor --;
			i --; // because neighbor[i] must be checked again.
		}
	}    
}

void TestParticle::makeNeighbor(){
    vector< vector<int> *> neighbor_cells;
	tsys->grid->get_neighbor_list_pointer( p, neighbor_cells);	
	neighbor.clear();
	foreach(vector< vector<int> *>, neighbor_cells, iter){
		foreach( vector<int>, *(*iter), i){
			if ( tsys->ct->connect(particle_number, *i ) == false  
				&& *i < particle_number ){
				neighbor.push_back( *i );
			}
		}
	}

}







