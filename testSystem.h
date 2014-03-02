//
//  testSystem.h
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#ifndef stodyn_testSystem_h
#define stodyn_testSystem_h
#include <string>
#include "vec3d.h"
#include "sdsystem.h" // System for Stokesian dynamics.
#include "common.h"
#include "bond.h"
#include "grid.h"
#include "Particle.h"

using namespace std;

class Particle;
class SDsystem;
class Bond;
class Grid;

class TestSystem {    
private:
    SDsystem *sd_sys;
    int num_of_particle;
    int n_bond;
    vector<Bond *> bond;

//    double *mov; // Mobility matrix of SD
//    double *vos; // velocity, omega, stresslet
//    double *fte; // force, torque, rate-of-strain

    /* Simulation box (particles) */
    double l_box[3]; // (lx, ly, lz)
    /* Simulation box (fluid) */
    double l_box_sd[3]; // (lx, ly, lz)


    double L_box_particle;
    double L_box_sd;

    double dist_bond_generate0; //distnce to generate initial bond;
    double dist_bond_generate1; //distnce to generate bond in dynamics.
    
    double time_simulation;
    
    char parameters_file[128];
    char simu_name[128];
    string rootdir;
    string bond0_file;
    string bond1_file;
    BondParameter bond0;
    BondParameter bond1;

    
    double shear_rate;
    double dt;
    double tmax;

    vector<int> regeneration_bond;
    vector<int> rupture_bond;
    
    void readParameterKey(const string &codeword,                               
                          const string &value);
    void setBondParameter(BondParameter &bondparameter, 
                          string &bond_file);
    void prepareBond(BondParameter & _bond);

    void makeInitialBond();

    void makeNeighbor();
    
    void generateBond();


public:
	TestSystem(SDsystem &sd_sys_);
	~TestSystem();
    
    
    
    void setParameterFile(char *parameter_file_);
    void importCluster(char* importfilename, int skipline);

    void readParameterFile();

    void readBondParameter();

    void setSimulationBoxs();
    
    void initDEM();

    vector<Particle *> particle;
    Grid *grid;
    ConTable *ct; // connection table 
    vec3d *pos; // positions of all particles 
    double sq_dist_generate; // square of distance for bond generation.
    
//    vector<testParticle *> particle;
//    vector<testBond *> bond;
//    ConTable *ct;
//    Grid *grid;
    
    
};
#endif
