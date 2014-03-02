/*
 *  bond.cpp
 *  CCN_3D
 *
 *  Created by seto on 09/08/12.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "bond.h"

Bond::Bond(int d0, int d1, DEMsystem &dem_){
    dem = &dem_;
	status = 1;
    initial_bond = dem->initialprocess;
	bond_number = dem->n_bond;
	dem->n_bond ++ ;
	d[0] = d0, d[1] = d1;  
	dem->ct->on_connect( d[0], d[1] );
	p_particle0 = dem->particle[d[0]];
	p_particle1 = dem->particle[d[1]];

	(*p_particle1).setBond(bond_number, d[0]);
	(*p_particle0).setBond(bond_number, d[1]);
	
    pp[0] = (*p_particle0).p_pos(); // pointer to postion of particle 0.
	pu[0] = (*p_particle0).pu_back(); // pointer to u vector of particle 0.

    pp[1] = (*p_particle1).p_pos();
	pu[1] = (*p_particle1).pu_back();

    r_vec = (*pp[1]) - (*pp[0]);

	p_tor[0] = (*p_particle0).p_tor_angle_back();
    p_tor[1] = (*p_particle1).p_tor_angle_back();
	(*p_tor[0]) = 0;
	(*p_tor[1]) = 0;

	r = r_vec.norm();
	e_normal = r_vec / r;

	(*pu[0]) = e_normal;
	(*pu[1]) = - e_normal;
    
    u_vector[0] = (*p_particle0).orientation.ori_backward(*pu[0]);
    u_vector[1] = (*p_particle1).orientation.ori_backward(*pu[1]);

    if (initial_bond){
        u_vector_initial[0] = u_vector[0];
        u_vector_initial[1] = u_vector[1];
    }

    if (initial_bond)
        para = dem->bond0;
    else
        para = dem->bond1;
    
    
    sq_fsc = sq(para.fsc);
    sq_mbc = sq(para.mbc);
    sq_mtc = sq(para.mtc);
    
    calcForce();
}

Bond::~Bond(){
    cerr << "deleted bond" << endl;
}

void Bond::whichparticle(int &i, int &j){
    i = (*p_particle0).particle_number;
    j = (*p_particle1).particle_number;
}

void Bond::addContactForce(){
	calcForce();

	(*p_particle0).stackForce( force0, torque0);
//    (*p_particle1).stackForce(force1, torque1);
	(*p_particle1).stackForce(-force0, torque1);
}

void Bond::calcForce(){
	/* calc basic information */
	r_vec = (*pp[1]) - (*pp[0]);
	r = r_vec.norm();
    e_normal = r_vec/r;

    u01 = *pu[1] - *pu[0];
    q = r - 2.;
    
    d_slid = u01 - dot(u01, e_normal)*e_normal;
    ang_bend = - cross(*pu[0], *pu[1]);
    ang_tort = - (*p_tor[0]) - (*p_tor[1]);
    
	/* calc bond stress */
	force_normal   = para.kn*q;
    force_sliding  = para.ks*d_slid;
    
	moment_bending = para.kb*ang_bend;
	moment_torsion = para.kt*ang_tort;

    
    
	/** compose **/
	/* force to particle 0 */
	force0 = force_normal*e_normal + force_sliding;
	  // force 1 = - force0  
    moment_sliding = cross(e_normal, force_sliding);
    torque_tmp = moment_bending + moment_torsion*e_normal;
	/* torque to particle 0 */
    torque0 = moment_sliding + torque_tmp;
	/* torque to particle 1 */
	torque1 = moment_sliding - torque_tmp;
}

void Bond::chPointer(int i, int particle_num){
	if (d[0] == particle_num){
		// It is rewritten in d[0].
		pu[0] = (*p_particle0).pu(i);
		p_tor[0] = (*p_particle0).p_tor_angle(i);
	} else {
        // It is rewritten in d[1].
		pu[1] = (*p_particle1).pu(i);
		p_tor[1] = (*p_particle1).p_tor_angle(i);
	}
}

void Bond::rupture(){
    status = 0;
	dem->ct->off_connect( d[0], d[1] );
	(*p_particle0).delConnectPoint(bond_number);
    (*p_particle1).delConnectPoint(bond_number);
	if ( d[1] < d[0] )
		(*p_particle0).addNeighbor( d[1] );
	else
		(*p_particle1).addNeighbor( d[0] );

}

void Bond::regeneration(){
    status = 2;
	r_vec = (*pp[1]) - (*pp[0]);
    r = r_vec.norm();
	e_normal = r_vec/r;
    (*pu[0]) = e_normal;
	(*pu[1]) = - e_normal;
    (*p_tor[0]) = 0.0;
	(*p_tor[1]) = 0.0;
    u_vector[0] = (*p_particle0).orientation.ori_backward(*pu[0]);
    u_vector[1] = (*p_particle1).orientation.ori_backward(*pu[1]);

    calcForce();
	//regeneration_tuzuku = true;
}


void Bond::output_normal(ofstream &out){
	if ( force_normal > 0) {
		out << "@ 2" << endl;
		out << "r " << 0.01*force_normal<< endl;			
	} else {
		out << "@ 3" << endl;
		out << "r " << 0.01*force_normal << endl;			
	}
	drawString(out);
}
void Bond::output_sliding(ofstream &out){
//	/**********************************/
//#ifdef NONEWBOND	
//	if (sticky_bond == false)
//		return;
//#endif
//	/**********************************/

	out << "r " << 0.01*force_sliding.norm() << endl;
	drawString(out);
}
void Bond::output_bending(ofstream &out){
	out << "r " << 0.01*moment_bending.norm() << endl;
	drawString(out);
}
void Bond::output_torsion(ofstream &out){
	out << "r " << 0.01* abs(moment_torsion) << endl;
	drawString(out);
}

void Bond::failureCondition(){
	/* Condition of bond failure
     *
     * The bond can be robuster by compression?
     * The normalization for des[0] should not be cp_f_n_max;
     * Now, this effect is neglected.
     */
    if ( q < 0 ){
        // force_normal < 0;
        //des[0] = force_normal * para.robust_bond_compression;
        des[0] = 0;
    } else {
        //des[0] = sq(force_normal/para.fnc);
        des[0] = 0;
    }
	//des[1] = force_sliding.sq_norm()/sq_fsc;
    des[1] = 0;
	des[2] = moment_bending.sq_norm()/sq_mbc;
	des[3] = sq(moment_torsion)/sq_mtc;
	D_function = des[0] + des[1] + des[2] + des[3];
	if ( D_function > 1. ){
        //if ( force_normal > para.fnc){
        //            dem->rupture_bond.push_back(bond_number);                            
        //} else{            
        dem->regeneration_bond.push_back(bond_number);
        //  }
        /*
		int max_des = 0;
		for (int j=1; j<4; j++){
			if (des[ max_des ] < des[j] ){
				max_des = j;
			}
		}
		switch (max_des) {
			case 0:	++ dem->rup_normal;	break;
			case 1:	++ dem->rup_shear; break;
			case 2:	++ dem->rup_bend; break;
			case 3:	++ dem->rup_torsion; break;
		}
         */
	}
	return;
}

void Bond::monitor_state(ofstream &out){
    vec3d u_init[2];
    vec3d r_vec_init;
    vec3d e_normal_init;
    u_init[0] = (*p_particle0).orientation.ori_forward(u_vector[0]);
    u_init[1] = (*p_particle1).orientation.ori_forward(u_vector[1]);
    r_vec_init = *pp[1] - *pp[0];
//	periodicBoundary_rvec(r_vec_init);
    double normal_distance;
    normal_distance = r_vec_init.norm();
	e_normal_init = r_vec_init/normal_distance;
	u01 = u_init[1] - u_init[0];
    
    double slide_distance = (u01 - dot(u01, e_normal_init)*e_normal_init).norm();
    double dot_prod = dot( u_init[0], - u_init[1]);
    double rolling_angle;
    if ( dot_prod > 0.9 ){
        double cross_prod = (cross( u_init[0], - u_init[1])).norm();
        rolling_angle = asin(cross_prod);
    }else{
        rolling_angle = acos(dot_prod);
    }
    double torsional_angle = 0; 
    out << bond_number << ' ' << status  << ' ';
    out << normal_distance << ' ' << slide_distance << ' ' << rolling_angle << ' ' << torsional_angle << endl ;
}


void Bond::drawString(ofstream &out){
    out << "s ";
    out << (*pp[0]).x - dem->lx0 << ' ';
    out << (*pp[0]).y - dem->ly0 << ' ';
    out << (*pp[0]).z - dem->lz0  << ' ';
    out << (*pp[0]).x + (*pu[0]).x - dem->lx0 << ' ';
    out << (*pp[0]).y + (*pu[0]).y - dem->ly0 << ' ';
    out << (*pp[0]).z + (*pu[0]).z - dem->lz0 << endl;
    out << "s ";
    out << (*pp[1]).x - dem->lx0 << ' ';
    out << (*pp[1]).y - dem->ly0 << ' ';
    out << (*pp[1]).z - dem->lz0 << ' ';
    out << (*pp[1]).x + (*pu[1]).x - dem->lx0 << ' ';
    out << (*pp[1]).y + (*pu[1]).y - dem->ly0 << ' ';
    out << (*pp[1]).z + (*pu[1]).z - dem->lz0 << endl;
}


void Bond::output_u(ofstream &out){
    /**********************************************/
    /*
     * The bond can be robuster by compression?
     * The normalization for des[0] should not be cp_f_n_max;
     * Now, this effect is neglected.
     */
    if ( force_normal > 0 )
        des[0] = sq( force_normal / dem->cp_f_n_max);
    else 
        des[0] = force_normal * dem->robust_bond_compression;
    
	des[1] = force_sliding.sq_norm() / sq(dem->cp_f_s_max);
	des[2] = moment_bending.sq_norm() / sq(dem->cp_m_b_max);
	des[3] = sq(moment_torsion / dem->cp_m_t_max);
    
	D_function = des[0] + des[1] + des[2] + des[3];
	int max_des = 0;
	for (int j=1; j<4; j++){
		if (des[ max_des ] < des[j] ){
			max_des = j;
		}
	}
	
	switch (max_des) {
		case 0:	out << "@ 3\n";break;
		case 1:	out << "@ 4\n";break;
		case 2:	out << "@ 5\n";break;
		case 3:	out << "@ 6\n";break;
	}
	out << "r " << 0.3*des[max_des] << endl;
	drawString(out);
}
