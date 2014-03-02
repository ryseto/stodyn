/*
 *  rigidClusterAnalysis.h
 *  stodyn
 *
 *  Created by seto on 10/10/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include "system.h"
#include "stokes.h"
#include "ewald-3fts.h"
#include "libiter.h"
#include "vec3d.h"
#include "my_utilities.h"
#include "common.h"
#include "sparticle.h"
#include "sbond.h"
class SBond;
class SParticle;
vector<SParticle *> sparticle;
vector<SBond *> sbond;

void splitCluster(int num_particle, int cutbond, vector<int> &div_cluster){
	for (int j=0; j < div_cluster.size() ; j++){
		if ( div_cluster[j] == num_particle ){
			return;
		}
	}
	div_cluster.push_back(num_particle);
	
	for (int i=0; i < sparticle[num_particle]->num_bond(); i++){
		int next = sparticle[num_particle]->neighbor_particle(i, cutbond);
		if ( next != -1 && next != num_particle )
			splitCluster( next, cutbond, div_cluster );
	}
}

void rigidClusterAnalysis(char** argv){
    System sy;    
    sy.type_simu = argv[1][0];	
	vector<vec3d> pos;
	vector<vec3d> force;
	vector<vec3d> torque;
    /* I need to rewrite the function to import data.*/
	sy.importTBMfile(argv[2], pos, force, torque);
    vec3d pcm(50, 50, 50);
	for (int i=0; i < sy.np; i++){
        sparticle.push_back(new SParticle(i, pos[i]));
    }
    int bond_number ;
    double delta = 1e-5;
    double d_sq_init = 4 + delta;
    double d_sq = d_sq_init;
    while (true){
        bond_number = 0;
        for (int i=0; i < sy.np-1; i++){
            for (int j=i+1; j < sy.np; j++){
                double d = sq_dist( sparticle[i]->p, sparticle[j]->p);
                if ( d < d_sq ){
                    vec3d bond_pos = 0.5*(sparticle[i]->p + sparticle[j]->p);
                    sbond.push_back( new SBond(bond_number, i, j, bond_pos) );
                    sparticle[i]->setbond(bond_number);
                    sparticle[j]->setbond(bond_number);
                    bond_number ++;
                }
            }
        }
        cerr << bond_number << endl;
        if ( bond_number == sy.np - 1){
            break;
        }
        if ( bond_number < sy.np - 1){
            d_sq += 0.1*delta;
        } else{
            d_sq -= 0.1*delta;    
        }
    }
	
	for (int i=0; i < pos.size(); i++){
		if ( sparticle[i]->num_bond() == 0 ){
			cerr << "particle which has no bond is found. " << endl;
		}
	}
	//split_cluster
	cerr << pos.size() << ' ' << sbond.size() << ' ' << bond_number << endl;	
	for(int cut_bond=0; cut_bond < sbond.size(); cut_bond++){
		vector<int> clusterA;
		vector<int> clusterB;
		vec3d forceA(0,0,0);
		vec3d momentA(0,0,0);
		splitCluster(sbond[cut_bond]->particle[0], cut_bond, clusterA);
		for (int j=0; j < clusterA.size(); j++){
			forceA+= force[clusterA[j]] ;
			vec3d dr = sparticle[clusterA[j]]->p - sbond[cut_bond]->pos;
			momentA+= cross(dr,	force[clusterA[j]]) + torque[clusterA[j]];	
		}
		vec3d forceB(0,0,0);
		vec3d momentB(0,0,0);
		splitCluster(sbond[cut_bond]->particle[1], cut_bond, clusterB);
		for (int j=0; j < clusterB.size(); j++){
			forceB+= force[clusterB[j]] ;
			vec3d dr = sparticle[clusterB[j]]->p - sbond[cut_bond]->pos;
			momentB+= cross(dr,	force[clusterB[j]])+ torque[clusterB[j]];
		}

		if( clusterA.size() + clusterB.size() != sy.np){
			cerr << sy.np << endl;
			cerr <<  clusterA.size() << " + " << clusterB.size() << endl;
			cout << "error splitting" << endl;
            exit(1);
		}		
		//for (int i=0; i<sy.np ; i++){
		//			cout << "c ";
		//			cout << sparticle[i]->p.x -50<< ' ';
		//			cout << sparticle[i]->p.y -50<< ' ';
		//			cout << sparticle[i]->p.z -50<< endl;
		//		}
		bool stop = false;
		if (abs(forceA.norm() - forceB.norm())/forceA.norm() >0.05){
			cerr << "force does not balance." << endl;
			cerr << abs(forceA.norm() - forceB.norm())/forceA.norm() << endl;
			cerr << forceA.norm() << ' ' << forceB.norm() << endl;
			cerr << clusterA.size() << ' ' <<  clusterB.size() << endl;
			stop = true;
		}
		if (abs(momentA.norm() - momentB.norm())/momentA.norm() >0.05){
			cerr << "moment does not balance." << endl;
			cerr << abs(momentA.norm() - momentB.norm())/momentA.norm() << endl;
			cerr << momentA.norm() << ' ' << momentB.norm() << endl;
			cerr << clusterA.size() << ' ' <<  clusterB.size() << endl;
			stop = true;
		}
		//if ( stop) exit(1);
		int pA = sbond[cut_bond]->particle[0];
		int pB = sbond[cut_bond]->particle[1];
		
		vec3d normal_vector = sparticle[pA]->p - sparticle[pB]->p;
		normal_vector.unitvector();

		if ( clusterA.size() < clusterB.size() ){
			double f_normal = dot(forceA, normal_vector);
			sbond[cut_bond]->force_normal = f_normal;
			sbond[cut_bond]->force_sliding = (forceA - f_normal*normal_vector).norm();
			double m_torsion = dot(momentA, normal_vector);
			sbond[cut_bond]->moment_torsion = abs(m_torsion);
			sbond[cut_bond]->moment_bending = (momentA - m_torsion*normal_vector).norm();
		} else {
			double f_normal = dot(forceB, normal_vector);
			sbond[cut_bond]->force_normal = -f_normal;
			sbond[cut_bond]->force_sliding = (forceB - f_normal*normal_vector).norm();
			double m_torsion = dot(momentB, normal_vector);
			sbond[cut_bond]->moment_torsion = abs(m_torsion);
			sbond[cut_bond]->moment_bending = (momentB - m_torsion*normal_vector).norm();
		}
		//forceA.cerr();
		//forceB.cerr();
		//momentA.cerr();
		//momentB.cerr();
		//sd->pos[i*3 + 0] 
	}
	
	/*

	cout << "y 6" << endl;
	cout << "@ 3" << endl;
	for (int i = 0; i < sbond.size(); i++){
		int p0 = sbond[i]->particle[0];
		int p1 = sbond[i]->particle[1];
		cout << "r " << 0.01*abs(sbond[i]->force_normal) << endl;
		cout << "s " ;
		cout << pos[p0].x -50<< ' ';
		cout << pos[p0].y -50<< ' ';
		cout << pos[p0].z -50<< ' ';
		cout << pos[p1].x -50<< ' ';
		cout << pos[p1].y -50<< ' ';
		cout << pos[p1].z -50<< ' ';
		cout << endl;
	}
	cout << "y 7" << endl;
	cout << "@ 4" << endl;
	for (int i = 0; i < sbond.size(); i++){
		int p0 = sbond[i]->particle[0];
		int p1 = sbond[i]->particle[1];
		cout << "r " << 0.01*sbond[i]->moment_bending << endl;
		cout << "s " ;
		cout << pos[p0].x -50<< ' ';
		cout << pos[p0].y -50<< ' ';
		cout << pos[p0].z -50<< ' ';
		cout << pos[p1].x -50<< ' ';
		cout << pos[p1].y -50<< ' ';
		cout << pos[p1].z -50<< ' ';
		cout << endl;
	}
*/
	double max_force_elong = 0;
	double max_force_comp = 0;
	double max_force_s = 0;
	double max_moment_b = 0;
	double max_moment_t = 0;
	
	for (int i = 0; i < sbond.size(); i++){
		if ( sbond[i]->force_normal > 0 ){
			if (max_force_elong < sbond[i]->force_normal)
				max_force_elong = sbond[i]->force_normal;
		} else {
			if (max_force_comp < abs(sbond[i]->force_normal))
				max_force_comp = abs(sbond[i]->force_normal);
		}            
		if (max_force_s < sbond[i]->force_sliding )
			max_force_s = sbond[i]->force_sliding;
		if (max_moment_b < sbond[i]->moment_bending )
			max_moment_b = sbond[i]->moment_bending;
		
		if (max_moment_t < sbond[i]->moment_torsion )
			max_moment_t = sbond[i]->moment_torsion;
	}
	cout << sy.np << ' ';
	cout << max_force_elong  << ' ';
	cout << max_force_comp  << ' ';
	cout << max_force_s  << ' ';
	cout << max_moment_b  << ' ';
	cout << max_moment_t  << ' ';
	cout << sy.gyration_radius << endl;
    
}
