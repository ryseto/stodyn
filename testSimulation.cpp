//
//  testSimulation.cpp
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include "testSimulation.h"
#include "sdsystem.h"
using namespace std;
void testSimulation(int argc, char** argv){
    if (argc == 2){
        cerr << "usage:" << endl;
        cerr << "stodyn T parameters cluster skip shearrate version" << endl;
        return;
    }
    SDsystem sd_sys;
    //sd_sys.type_simu = argv[1][0];
    DEMsystem dem(sd_sys);

    dem.initialprocess = true;
    /* parameter file */
	dem.setParameterFileDEM(argv[2]);   
    /* Cluster configuration : num. of skip lines */
	dem.importCluster(argv[3], atoi(argv[4]));
    dem.setVersion(argv[5]);
    /* Method of Hydrodynamic interaciton FDA(0), SD_nolub(1), SD_lub(2)*/
	dem.readParameterFileDEM();
    dem.readBondParameter();
    dem.setSimulationBoxs();
    /* set imposed flow */
    sd_sys.setFlowType('s');
	/* DEM */
	dem.setDEMParameters();
	dem.initDEM();
    
//    
//    SDsystem sd_sys;
//    TestSystem tsimu(sd_sys);
//    
//    tsimu.setParameterFile(argv[2]);   
//    tsimu.importCluster(argv[3], atoi(argv[4]));
//    tsimu.readParameterFile();
//
//    tsimu.readBondParameter();
//    tsimu.setSimulationBoxs();
//
//    sd_sys.setFlowType('s'); // s = shear flow
//
//	tsimu.initDEM();
//
//    /* Stokesian Dynamics */
////    sy.initLibStokes();
//
//    
//    cerr << "here" << endl;
//    
    return;
}








