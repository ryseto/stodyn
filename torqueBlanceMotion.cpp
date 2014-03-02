//
//  torqueBlanceMotion.cpp
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "torqueBlanceMotion.h"
using namespace std;

void torqueBalancedMotion(int argc, char** argv){
    System sy;
    sy.type_simu = argv[1][0];	

    
    vec3d torque_free_U;
    vec3d torque_free_O;
    /*
	 * set parameters
	 */ 
    sy.setBox(100, 100, 100);
	sy.importCluster(argv[2], atoi(argv[3])); 
    sy.setMethodHydrodynamicInteraction(atoi(argv[4]));
	
    sy.type_of_flow = 's';
    sy.set_output_precision(9);
    /*
     * Compose Rigid resistance matrix
     */
    if (sy.method_hydroint >= 0){
        sy.initLibStokes();
    }
    sy.calcResistanceMatrixRigidAggregate();
    
    /*
     * Solve force- and torque-free motion in shear flow.
     */
    sy.setSimpleShearFlow(1.0);
    sy.solveForceFreeRigidMotion(torque_free_U, torque_free_O);
    /*
     * Calculate force/torque/stresslet
     */ 
    cerr << "Uag: "; torque_free_U.cerr();
    cerr << "Oag: "; torque_free_O.cerr();
    sy.calcFTS_RigidAggregate(torque_free_U, torque_free_O);
    sy.outFTS_RigidAggregate();
    sy.calcFTS_GrandResistanceMatrix();
    /*
     *
     */
    sy.calcGyrationRadius();
    export_profile(sy, false);
    //sy.openOutputFileStream();    
    return;
}

void export_profile(System &sy, bool export_GRM){
    string common_name;
	switch (sy.method_hydroint ) {
		case 0:
			common_name = "FD_";
			break;
		case 1:
			common_name = "SD_nolub_";
			break;
		case 2:
			common_name = "SD_lub_";
			break;
	}
    common_name += sy.filename;
    
    string outfilename_TMB;
    outfilename_TMB = "TBM_" + common_name + ".dat";
	ofstream output_tbm;
	output_tbm.open(outfilename_TMB.c_str());
	sy.output_TorqueBalancedMotion(output_tbm);	
	output_tbm.close();
    
    string outfilename_RRA;
    outfilename_RRA = "RRA_" + common_name + ".dat";
    ofstream output_rra;
    output_rra.open(outfilename_RRA.c_str());
    sy.output_ResistanceMatrixRigidAggregate(output_rra);
    output_rra.close();
    
    if (export_GRM){
        string outfilename_GRM;
        outfilename_GRM = "GRM_" + common_name + ".dat";
        fstream output_grm;
        
        output_grm.open(outfilename_GRM.c_str(),
                        ios::out | ios::binary);
        sy.output_GrandResistanceMatrix(output_grm);
        output_grm.close();
    }
    
}


