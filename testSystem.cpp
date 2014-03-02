//
//  testSystem.cpp
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "testSystem.h"
#include "vec3d.h"
#include "contable.h"
#include "grid.h"
using namespace std;

TestSystem::TestSystem(SDsystem &sd_sys_){
    sd_sys = &sd_sys_;
    ct = new ConTable;
    grid = new Grid;
	n_bond = 0 ; /* number of bond */
//    restructuring = true;
}

TestSystem::~TestSystem(){
    delete [] pos;
}


void TestSystem::setParameterFile(char *parameter_file_){
    sprintf(parameters_file, "%s", parameter_file_);
    string s_parameters_file = parameters_file;
       int i_backslash = s_parameters_file.find_last_of( "/") + 1;
       int i_extention = s_parameters_file.find( ".txt" );	
        sprintf(simu_name, "%s",
                (s_parameters_file.substr( i_backslash, i_extention-i_backslash)).c_str());
        cerr << " simu_name = " << simu_name << endl;

}


void TestSystem::importCluster(char* importfilename, int skipline){
    vector<vec3d> tmp_import_data;
    
    ////////////////////////////////////////////////////////////
    // Import positions of particles
	ifstream fin;
	fin.open( importfilename );
    // Skip lines
	char buf[1000];
	for (int i = 0; i< skipline; i++){
		fin.getline( buf, 1000 );
	}
    double x, y, z;
	do{
		fin >> x >> y >> z;
        if(fin.fail())
            break; 
        vec3d new_p(x, y, z);
		tmp_import_data.push_back(new_p);
	} while (!fin.eof());
    fin.close();
    //
    ////////////////////////////////////////////////////////////
    /*
     * Set the imported data in object;
     */
    num_of_particle = tmp_import_data.size();
    pos = new vec3d [num_of_particle];
    for (int i = 0; i < num_of_particle ; i++){
        pos[i] = tmp_import_data[i];
    }	
}

void TestSystem::readParameterFile(){
	ifstream fin;
	fin.open(parameters_file);
	string codeword, value, unit;
	while (!fin.eof()){
		fin >> codeword ;
		if ( codeword == "#") {
			char buf[1024];
            fin.get(buf, 1024);
		} else if (codeword == "!"){
			break;
		} else {
			fin >> value;
            if (fin.eof()){
                break;   
            }
			readParameterKey(codeword, value);
		}
	}
    fin.close();
    cerr << "completed to read parameters." << endl;
}

void TestSystem::readParameterKey(const string &codeword,                               
                                 const string &value){
	map<string,int> keylist;
    keylist["rootdir:"]=1; const int _rootdir = 1;
    keylist["bond0_file:"]=2; const int _bond0_file = 2;
    keylist["bond1_file:"]=3; const int _bond1_file = 3;
    keylist["method_hydroint:"]=4; const int _method_hydroint = 4;
    keylist["shear_rate:"]=5; const int _shear_rate = 5;
    keylist["dt:"]=6; const int _dt = 6;
    keylist["tmax:"]=7;const int _tmax = 7;
    keylist["L_box_particle:"]= 8; const int _L_box_particle = 8;
    keylist["L_box_sd:"]= 9; const int _L_box_sd = 9;

	cerr << "codeword: " << codeword << ' ' << value  << endl;
	switch(keylist[codeword]){
        case _rootdir: 
            rootdir = value; 
            break;
        case _bond0_file: 
            bond0_file = value; 
            break;
        case _bond1_file: 
            bond1_file = value; 
            break;
        case _method_hydroint:
            sd_sys->method_hydroint = atoi(value.c_str());
            break;
        case _shear_rate: 
            shear_rate = atof(value.c_str()); 
            break;
        case _dt:
            dt = atof(value.c_str() );
            break; 
        case _tmax:
            tmax = atof(value.c_str());
            break;
        case _L_box_particle:
            L_box_particle = atof(value.c_str());
            break;
        case _L_box_sd:
            L_box_sd = atof(value.c_str());
            break;
		default:
			cerr << "The codeword " << codeword << " is'nt associated with an parameter" << endl;
			exit(1);
	}
}


void TestSystem::readBondParameter(){
    setBondParameter(bond0, bond0_file);
    setBondParameter(bond1, bond1_file);
    dist_bond_generate0 = bond0.dist_generate;
    dist_bond_generate1 = bond1.dist_generate;
}

void TestSystem::setBondParameter(BondParameter &bondparameter, string &bond_file){
    string path_bond = rootdir + "/" + bond_file;
    ifstream fin;
    if(fin.is_open())
        fin.close();
    fin.open(path_bond.c_str());
    string codeword;
    fin >> codeword >> bondparameter.fnc;
    fin >> codeword >> bondparameter.fsc;
    fin >> codeword >> bondparameter.mbc;
    fin >> codeword >> bondparameter.mtc;
    fin >> codeword >> bondparameter.n_max;
    fin >> codeword >> bondparameter.s_max;
    fin >> codeword >> bondparameter.b_max;
    fin >> codeword >> bondparameter.t_max;
    fin >> codeword >> bondparameter.overlap_force_factor;
    fin >> codeword >> bondparameter.robust_bond_compression;
    fin >> codeword >> bondparameter.dist_generate;
    fin.close();
    prepareBond(bondparameter);
    return;
}

void TestSystem::prepareBond(BondParameter & _bond){
    /* calculate spring constants
     * as critical forces / critical lengths.
     */
    double a0 = 1.0;
    _bond.kn = _bond.fnc / _bond.n_max;
    _bond.ks = _bond.fsc / _bond.s_max;
    _bond.kb = _bond.mbc / (_bond.b_max * a0 * a0);
    _bond.kt = _bond.mtc / (_bond.t_max * a0 * a0);
}

void TestSystem::setSimulationBoxs(){
    /* set box size for Stokesian dynamics
     */
    //sd_sys->setBox(L_sd, L_sd, L_sd);
    l_box_sd[0] = l_box_sd[1] = l_box_sd[2] = L_box_sd;
    /* set box size for DEM simulation
     */
    l_box[0] = l_box[1] = l_box[2] = L_box_particle;
    cerr << l_box[0] << endl;
}

void TestSystem::initDEM(){
    /* 
     */
    ct->set(num_of_particle);
    double grid_size = 2.2;

    grid->init(num_of_particle, l_box, grid_size);
    //  vec3d origin_shift(lx0,ly0,lz0);
    
    
//	for (int i=0; i < num_of_particle ; i++){
//		particle.push_back(new Particle(i));
//	}
    for (int i=0; i < num_of_particle ; i++){
		particle[i]->p = pos[i];
	}
    
	rupture_bond.clear();
	regeneration_bond.clear();
  
	makeInitialBond();
    cerr <<"here" << endl;
    time_simulation = 0.; //set the time clock zero.
//    vos = new double [np*11];
//    fte = new double [np*11];
//    pos = new vec3d [np];
//    set_FDA();
//    pos_center_of_mass.set(lx0, ly0, lz0);
}

void TestSystem::makeInitialBond(){
	sq_dist_generate = sq(dist_bond_generate0);
	makeNeighbor();
	generateBond();	
	sq_dist_generate = sq(dist_bond_generate1);
}

void TestSystem::makeNeighbor(){
    grid->clear_vcell();
    for (int i=0; i < num_of_particle; i++){
        grid->remakeOneParticle(i, particle[i]->p);
    }
    cerr << "-3-";
    for (int i=0; i < num_of_particle; i++){
		particle[i]->makeNeighbor();
	}
    cerr << "-4-";
}

void TestSystem::generateBond(){
	for (int i=0; i < num_of_particle; i++){
		particle[i]->generateBond();
	}
}




