//
//  calcDragAndTorque.h
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#ifndef stodyn_calcDragAndTorque_h
#define stodyn_calcDragAndTorque_h
#include "system.h"


void calcDragAndTorque(int argc, char** argv){    
    System sy;

    sy.type_simu = argv[1][0];
	/* Shear-rate or typical velocity is set one.
	 *
	 */
	if (argc == 2){
		/* POSITIONS indicates the path to a file including (x,y,z) of particles
		 *
		 */
		cerr << "Usage: stodyn u POSITIONS skipline method" << endl;
		cerr << "Usage: stodyn s POSITIONS skipline method" << endl;
        cerr << "method=1 : Stokesian dynamics without lublication correction" << endl;
		cerr << "method=2 : Stokesian dynamics with lublication correction " << endl;
        cerr << "This calculation will be conducted by dimensionless variables." << endl;
        cerr << "The unit of length is given by the radius of particle." << endl;
        cerr << "The unit of velocity is given by the velocity of uniform flow." << endl;
        cerr << "Or, the unit of velocity is given by the product of the shear rate and the unit of length." << endl;
        cerr << "Here, U=1.0 is used for uniform flow." << endl;
        cerr << "G=1.0 is used for shear flow." << endl;
        cerr << "F=U is expected for one-body solution." << endl;
        cerr << "===================" << endl;
        cerr << "For rescaling:" << endl;
        cerr << "* Uniform flow. " << endl;
        cerr << " a and U are given." << endl;
        cerr << "L0 = a, U0 = U, F0 = 6*M_PI*eta*L0*U0" << endl;
        cerr << "* Shear flow. " << endl;
        cerr << " a and G are given." << endl;
        cerr << "L0 = a, U0 = G*a, F0 = 6*M_PI*eta*L0*U0" << endl;
        cerr << "-------------------" << endl;
        cerr << "F ---> F*F0 " << endl;
        cerr << "T ---> T*F0*L0 " << endl;
        cerr << "-------------------" << endl;
		return;
	}
    
    sy.method_hydroint = atoi(argv[4]); 
    switch( sy.method_hydroint ){
		case 1:
            sy.setLubrication(0);
			break;
		case 2:
			sy.setLubrication(1);
			break;
        default:
            cerr << "!!!!!!! Warning !!!!!!!!!!!!" << endl;
            cerr << "4th argument should be 1 or 2:" << endl;
            cerr << "1: without lubrication" << endl;
            cerr << "2: with lubrication" << endl;
            cerr << "exit" << endl;
            exit(1);
	} 
    
    sy.set_output_precision(12);
	sy.type_of_flow = sy.type_simu;
	/* init() is called in this function*/
    sy.setBox(100, 100, 100);
	sy.importCluster(argv[2], atoi(argv[3]));
    sy.initLibStokes();
    //    sy.setSD_IterationMethod();
	/*
	 * Set imposed flow 
	 */ 
	if (sy.type_of_flow == 'u'){
		vec3d U(1.0, 0.0, 0.0);	
		sy.setImposedFlow_uniform(U);
	} else {
		sy.setSimpleShearFlow(1.0);
	}
	
	/*
	 * Set velocities of particles
	 */ 
    sy.cl_velocity.set(0.0, 0.0, 0.0);
    sy.cl_omega.set(0.0, 0.0, 0.0);
    sy.setRigidVelocities();
    
	sy.setSD_IterationMethod();


	// main calculation
    double t_start,t_end, calc_time;
    //    t_start = gettimeofday_sec();
	sy.StokesianDynamics();

    //    t_end = gettimeofday_sec();
	calc_time = t_end - t_start;
    //output data
	sy.openOutputFileStream();
    
	//sy.calcFreeDrainingForce();
    //sy.FreeDrainingApproximation();
    
	sy.calcTotalForceTorqueStress();


    sy.outputForComparison();
    exit(1);
    //sy.outputSampleData();
	
}



#endif
