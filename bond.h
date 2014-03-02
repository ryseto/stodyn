/*
 *  bond.h
 *  CCN_3D
 *
 *  Created by seto on 09/08/12.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef bond_h
#define bond_h 1
#include "common.h"
#include <cstdlib>
#include <cmath>
#include "contable.h"
#include "particle.h"
#include "system.h"
#include "demsystem.h"
#include "grid.h"
class Particle;
class DEMsystem;

class Bond {
	int bond_number;
    BondParameter para;
    DEMsystem *dem;
    double q;
	double sq_rv;
	vec3d r_vec;
	vec3d force0;
	vec3d torque0;
	vec3d torque1;
	double force_normal;
	vec3d moment_bending;
	vec3d force_sliding;
	vec3d moment_sliding;
	double moment_torsion;
    vec3d torque_tmp;
    vec3d e_normal;
	vec3d u01;
    vec3d d_slid;
    vec3d ang_bend;
    double ang_tort;
    
	int d[2];
	double *p_tor[2];
	vec3d *pp[2];
	vec3d *pu[2];
	Particle *p_particle0;
	Particle *p_particle1;
    vec3d u_vector[2]; 
    vec3d u_vector_initial[2]; 
	double des[4];
    double sq_fsc;
    double sq_mbc;
    double sq_mtc;
private:
	void setting();
	void checkRuptureType(vec3d &, vec3d &, double);
	void periodicBoundary_rvec();
public:
	Bond();
	~Bond();
	Bond(int d0, int d1, DEMsystem &dem);
    /*
     * distance between two particles 
     */
    double r; 
    
    /* status; 
     * 0: The bond is ruptured.
     * 1: Only elastic deformation after the generation.
     * 2: Regenerations are experienced. 
     */
    int status; 
    bool initial_bond;
    double D_function;
    
	void addContactForce();
	void calcForce();
    void monitor_state(ofstream &out);
	void outputCompression(ofstream &out, double f_max);
	void outputTraction(ofstream &out, double f_max);
	void output_u(ofstream &out);
	void output_normal(ofstream &out);
	void output_sliding(ofstream &out);
	void output_bending(ofstream &out);
	void output_torsion(ofstream &out);
	void drawString(ofstream &out);
	void rupture();	
	void chPointer(int i, int particle_num);	
	void regeneration();
	void failureCondition();
    inline double val_F_norm(){ return force_normal;} 
    inline double val_F_slid(){ return force_sliding.norm();}
    inline double val_M_bend(){ return moment_bending.norm();}
    inline double val_M_tors(){ return abs(moment_torsion);}
    
    
    inline double valBendingAngle(){return ang_bend.norm();}
    inline double valTorsionalAngle(){return abs(ang_tort);}
    inline double valSlidingDisplacement(){return d_slid.norm();}
	inline double strenghOfForce(){return force0.norm();}
    void whichparticle(int &i, int &j);
};
#endif
