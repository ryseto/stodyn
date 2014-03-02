/*
 *  particle.cpp
 *  CCN_3D
 *
 *  Created by seto on 09/08/11.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "particle.h"

Particle::Particle(int particle_number_, 
                   DEMsystem &dem_,
                   SDsystem &sd_sys_   ){
	dem = &dem_;
    sd_sys = &sd_sys_;
	setInitial(particle_number_);
	setVelocityZero();
    orientation.set(1.,0.,0.,0.);
    resetForce();
}

//Particle::Particle(int particle_number_){
//
//    particle_number = particle_number_;
//    setVelocityZero();
//    orientation.set(1.,0.,0.,0.);
//    resetForce();
//    cn_size = 0;
//}

//Particle::Particle(int particle_number_){
//    setInitial(particle_number_);
//}

void Particle::setInitial(int particle_number_){
	particle_number = particle_number_;
	cn_size = 0;
}

void Particle::setPosition(const vec3d & position){
	p = position;		
}

void Particle::setPositionFromArray(const double *pos){
	p.x = *pos;
	p.y = *(pos+1);
	p.z = *(pos+2);
}

void  Particle::makeNeighbor(){
	vector< vector<int> *> neighbor_cells;
	dem->grid->get_neighbor_list_pointer( p, neighbor_cells);	
	neighbor.clear();
	foreach(vector< vector<int> *>, neighbor_cells, iter){
		foreach( vector<int>, *(*iter), i){
			if ( dem->ct->connect(particle_number, *i ) == false  
				&& *i < particle_number ){
				neighbor.push_back( *i );
			}
		}
	}
}

void Particle::generateBond(){
    int n_neighbor = neighbor.size();
	for (int i=0; i < n_neighbor; i++){
		if (sq_dist(p, dem->particle[neighbor[i]]->p) <= dem->sq_dist_generate ){
			dem->bond.push_back(new Bond ( neighbor[i], particle_number, *dem));            
			neighbor[i] = neighbor.back();
			neighbor.pop_back();
            n_neighbor --;
			i --; // because neighbor[i] must be checked again.
		}
	}    
}

void Particle::delConnectPoint(int bond_number){
	for (int i=0; i < cn_size-1 ; ++i){
		if ( cn[i].bond == bond_number ){			
			cn[i] = cn[ cn_size-1 ];
			dem->bond[ cn[i].bond ]->chPointer(i, particle_number);
			break;
		}
	}
	cn_size --;
	return;
}

void Particle::move_with_velocity(){
    // velocity is given by unit of SD, i.e. 
    p += velocity * dem->timestep;
    d_rotation = omega * dem->timestep;
    orientation.infinitesimalRotation( d_rotation );
    for (int i = 0; i < cn_size ; ++i){
        cn[i].u.rotateInfinitesimal( d_rotation );
        cn[i].tor_angle += dot(d_rotation, cn[i].u);
    }
		
}

//void Particle::output(ofstream &fout){
	//fout << "c " << p.x - dem->pos_center_of_mass.x ;
	//fout << ' '  << p.y - dem->pos_center_of_mass.y ;
	//fout << ' '  << p.z - dem->pos_center_of_mass.z << endl;
    //fout << "c " << p.x - dem->pos_center_of_mass.x ;
	//fout << ' '  << p.y - dem->pos_center_of_mass.y ;
	//fout << ' '  << p.z - dem->pos_center_of_mass.z << endl;
    
//}
