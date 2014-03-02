//
//  testParticle.h
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#ifndef stodyn_testParticle_h
#define stodyn_testParticle_h

#include "testSystem.h"
#include "sdsystem.h"
#include "common.h"

#include "grid.h"
#include "bond.h"
#include "vec3d.h"
#include "quaternion.h"
#include "contable.h"
class TestSystem;
class Grid;
class Bond;
class SDsystem;

class TestParticle{	
private:
	int cn_size;
	vec3d velocity_heun1;
	vec3d omega_heun1;
	vec3d d_rotation;
	vector<int> neighbor;
	ConnectPoint cn[12];
    SDsystem *sd_sys;
    TestSystem *tsys;
	int i_array[3];

protected:
	void (vec3d::*p_change)(double, double, double);
	void remove_list_neighbor_all();
public:
	TestParticle(int particle_number_, 
                 TestSystem &tsys_, 
                 SDsystem &sd_sys_);
    
    TestParticle(int particle_number_, 
             TestSystem &tsys_);
	TestParticle(int paritcle_number_);
	~TestParticle(){
	}
    // positions of particles
	vec3d p;
    quaternion orientation;
    void set_qu();
	vec3d velocity;
	vec3d omega;
    double stresslet[5];
    /* This force means the total forces from
     * the contacting particles.
     */
    vec3d force; 
	vec3d torque;
	double sq_force;
	void set_force_SD();
	void set_SD();
	int particle_number;
	void makeNeighbor();
	void addNeighbor( int neighbor_particle ){
		neighbor.push_back( neighbor_particle );
	}
	inline void setBond(int bond_number, int next_particle){
		cn[cn_size].next = next_particle;
		cn[cn_size].bond = bond_number;
		cn[cn_size].tor_angle = 0;
		++ cn_size;
	}
	void addFluidForce();
	void addFluidForce_stillfluid();
	void addGravityForce();
	double tmp_u_abs(){
		return cn[0].u.norm();
		
	}
    void move_with_velocity();
    
	inline vec3d vectorForce(){
		return force;
	}
	double abs_sq_velocity(){
		return velocity.sq_norm();
	}
	
	double abs_sq_omega(){
		return omega.sq_norm();
	}
    
	inline void resetForce(){
		force.reset();
		torque.reset();
	}
	inline void resetFz(){
		force.z = 0;
	}
	
	inline void resetTorque(){
		torque.reset();
	}
	inline double valForceZ(){
		double fz = force.z;
		force.reset();
		return  fz;
	}
	
    /* This is called from contacting bonds.
     * The total contact force acting on this particle
     * will be calculated.
     */
	inline void stackForce(const vec3d &force_, const vec3d &torque_){
		force += force_;
		torque += torque_;
	}
	
	void setNorm_u(){
		for (int i = 0; i < cn_size ; ++i){
			cn[i].u.unitvector();
		}
	}
    
    void setVelocityZero(){
		velocity.reset();
		omega.reset();
	}
	void delConnectPoint(int bond_number);	
	inline double x(){return p.x;}
	inline double y(){return p.y;}
	inline double z(){return p.z;}	
    
	void z_shift( double dz ){ p.z += dz; }
	inline vec3d pos(){return p;}	
	void setInitial(int disk_number_);
	void setPosition(const vec3d &position);
	void setPositionFromArray(const double *pos);
	inline double distOverlap(const vec3d &pp){
		return dist(p, pp) - 2.0; // ro = 2a = 1.0
	}
	//vec3d *p_pos(){ return &p;}
	vec3d *pu_back(){return &( cn[cn_size-1].u );}
	vec3d *pu(int i){return &( cn[i].u ); }
	double *p_tor_angle_back(){return &( cn[cn_size-1].tor_angle ); }
	double *p_tor_angle(int i){return &( cn[i].tor_angle ); }
	
	void generateBond();
    
	void cerr(){
		std::cerr << "c " << p.x << ' ' << p.y << ' ' << p.z << std::endl;
		//fout << particle_number << ' ' << p.x << ' ' << p.y << ' ' << p.z << endl;
	}
    
	double valForce(){ return force.norm(); }
	double valTorque(){ return torque.norm(); }
	double valVelocity(){ return velocity.norm(); }
	double valOmega(){ return omega.norm(); }
	double valOmegaY(){ return omega.y; }
	
	int valCn_size(){ return cn_size; }
	
	
	
};


#endif
