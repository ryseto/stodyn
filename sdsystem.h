//
//  sdsystem.h
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#ifndef stodyn_sdsystem_h
#define stodyn_sdsystem_h
#include <string>
#include "stokes.h" //Ichiki-san
#include "ewald-3fts.h" //Ichiki-san
#include "ewald-3fts-matrix.h" //Ichiki-san
#include "ewald-3ft.h" //Ichiki-san
#include "libiter.h" //Ichiki-san
#include "vec3d.h"

class SDsystem {
private:
    char type_of_flow;

public:
	SDsystem();
	~SDsystem();
    void setBox(double lx_, double ly_, double lz_);
    void setFlowType(char type_of_flow_);

    void initLibStokes();
    void set_dr_from_sdpos();
	void setPosition(int, const vec3d &);

    
    /*************************
	 *   Stokesian Dynamics
	 *************************/
    struct stokes *sd;
    void calcGrandMovMatrix(double* mov);

    /* method_hydroint
     * 0: Free-draining approximation
     * 1: Stokesian dynamics without lubrication
     * 2: Stokesian dynamics with lubrication
     */
    int method_hydroint;
    int np;
	int nm;
    int lubrication;
    vec3d *dr;
	double * velocity;
	double * omega;
	double * strain_velocity;
	double * force;
	double * torque;
	double * stresslet;
    double lx;
    double ly;
    double lz;
    double lx0;
    double ly0;
    double lz0;
    
};

#endif
