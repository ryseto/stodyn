//
//  sdsystem.cpp
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//
#include <iostream>
#include <iomanip>
#include "sdsystem.h"
using namespace std;

SDsystem::SDsystem(){
	/* lubrication ----> set by profile */
    cerr << "SDsystem()" << endl;
    np = 0;
	lubrication = -1;
	dr = NULL;
	velocity = NULL;
	omega = NULL;
	strain_velocity = NULL;
	stresslet = NULL;
	force = NULL;
	torque = NULL;
}


SDsystem::~SDsystem(){
    cerr << "~SDsystem()" << endl;

}

void SDsystem::setBox(double lx_, double ly_, double lz_){
	lx = lx_;
	ly = ly_;
	lz = lz_;
	lx0 = lx/2.0;
	ly0 = ly/2.0;
	lz0 = lz/2.0;
}

void SDsystem::setFlowType(char type_of_flow_){
    /* s : shear flow
     */
    type_of_flow = type_of_flow_;
}

/* Initialize libstokes
 * Basic parameters for the libstorks
 * and pos[] are set from init_aggregate.
 */
void SDsystem::initLibStokes(){
    nm = np;	/* nm : number of mobile particles  */
    cerr << "number of the particles : " << np << endl;
    sd = stokes_init();
	sd->twobody_lub = lubrication;
    sd->version = 2; /* 0 = F, 1 = FT, 2 = FTS  */
    sd->periodic = 0; 	/* 0 = non periodic, 1 = periodic */
	stokes_set_np(sd, np, nm);
    if (lx == 0 || ly == 0 || lz == 0 ){
        cerr << "lx, ly, lz are not given." << endl;exit(1);
    }
    cerr << "=======" << endl;
	stokes_set_l(sd, lx, ly, lz);
    
    int n3=np*3;
    int n5=np*5;
    try{
        velocity = new double [ n3 ];
        omega = new double [ n3 ];
        strain_velocity = new double [ n5 ];
        force = new double [ n3 ];
        torque = new double [ n3 ];
        stresslet = new double [ n5 ];
        dr = new vec3d [np];
    } catch (bad_alloc &){
        cerr << "bad_alloc at System::init()" << endl;
        exit(1);
    }
    
//    /*
//	 * set pos[] for SD
//	 */ 
//    
//    if (!init_aggregate.empty()){
//        double scale = 1.0;
//        for (int i = 0; i < np; i++){
//            setPosition(i, 
//                        scale*init_aggregate[i].x,
//                        scale*init_aggregate[i].y,
//                        scale*init_aggregate[i].z);
//        }
//        set_dr_from_sdpos();
//    }
}


void SDsystem::set_dr_from_sdpos(){    
    for (int i = 0; i < np; i++){
		dr[i].set(sd->pos[i*3 + 0] - lx0,
				  sd->pos[i*3 + 1] - ly0,
				  sd->pos[i*3 + 2] - lz0);
	}  
}

void SDsystem::setPosition(int i, const vec3d &position){
	int j = 3*i;
	sd->pos[j] = position.x;
	sd->pos[j+1] = position.y;
	sd->pos[j+2] = position.z;
}


// (U,O,S) = mov*(F,T,E)
void SDsystem::calcGrandMovMatrix(double* mov){
    calc_mob_3fts_matrix(sd, mov);
}


