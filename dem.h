/*
 *  dem.h
 *  stodyn
 *
 *  Created by seto on 10/07/07.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include "demsystem.h"
#include "common.h"
#include "bond.h"
#include "particle.h"
#include "SDsystem.h"
#include <vector>
using namespace std;

void demSimulation(int argc, char** argv){
    SDsystem sd_sys;    
//    sd_sys.type_simu = argv[1][0];	
    
    DEMsystem dem(sd_sys);
    if (argc == 2){
        cerr << "usage:" << endl;
        cerr << "stodyn D parameters cluster skip shearrate version" << endl;
        return;
    }
    dem.initialprocess = true;
    /* parameter file */
	dem.setParameterFileDEM(argv[2]);   
    /* Cluster configuration : num. of skip lines */
	dem.importCluster(argv[3], atoi(argv[4]));
    dem.setVersion(argv[5]);
    	dem.readParameterFileDEM();
    dem.readBondParameter();
    dem.setSimulationBoxs();
    /* set imposed flow */
    sd_sys.setFlowType('s');

	/* DEM */
	dem.setDEMParameters();
	dem.initDEM();
    sd_sys.np = dem.np;
    /* Stokesian Dynamics */
    sd_sys.initLibStokes();
	dem.openOutputFileDEM();
	dem.outputYaplotDEM();
    dem.outputDeformationConf();

	ofstream out_povray ;
    dem.setInitialMotion();
    dem.shiftClusterToCenter();
    dem.getMovMatrix();
    dem.initialprocess = false;
    sd_sys.set_dr_from_sdpos();
    dem.outputConfiguration();
    
    for(int m = 0; m < dem.stepnumber_shear_vs_bond; m++){
        
        if (m == 0){
            dem.shear_vs_bond = dem.shear_vs_bond_min;
        } else {
            dem.shear_vs_bond *= dem.ratio_shear_vs_bond;
        }
        cerr << "shear_vs_bond = "  << dem.shear_vs_bond << endl;

        dem.setTimeStep();
        double time_to_keep_shear = dem.time_simulation + dem.time_each_shear_vs_bond;
        do{
            double nexttime_output = dem.time_simulation + dem.time_interval_output;
            dem.makeNeighbor();
            int i =0;
            while (dem.time_simulation < nexttime_output){
                
                dem.calcInterParticleForce();
                
                if (sd_sys.method_hydroint == 0){
                    dem.freeDrainingApproximation();
                } else {
                    dem.calcVelocityByMov();
                }
                
                dem.checkFailure();
                if (!dem.regeneration_bond.empty()){
                    dem.regeneration_onebyone();
                    dem.regeneration_bond.clear();
                }
                
                dem.generateBond();
                
                if (sd_sys.method_hydroint != 0){
                    //dem.outputYaplotDEM();
                    dem.estimateClusterRotation();

                    if ( dem.step_deformation > dem.critical_deformation_SD ){
                        cerr << "renew matrix" << ' ' << dem.step_deformation << endl;
                        dem.shiftClusterToCenter();                    
                        dem.getMovMatrix();
                        dem.resetDeformation();
                        sd_sys.set_dr_from_sdpos();
                        if (dem.version[0] == 'o'){
                            dem.outputYaplotDEM();
                            dem.outputDeformationConf();
                        }
                        if (dem.step_deformation > 100){
                            cerr << "the simulation failed." << endl;
                            exit(1);
                        }                        
                    }
                }      
                i ++; 
                if ( i % 200 == 0){
                    dem.shiftClusterToCenter(); 
                    dem.q_rot.normalize();
                    for (int i=0; i < sd_sys.np; i++){
                        dem.particle[i]->setNorm_u(); /*********** IMPORTANT ************/
                        dem.particle[i]->orientation.normalize(); /*********** IMPORTANT ************/
                    }
                }
             //   dem.outputYaplotDEM();

            }
            cerr << "t = " << dem.time_simulation << endl;
            dem.shiftClusterToCenter();
            dem.calcGyrationRadius();
            dem.checkState();
            dem.calcTotalFTS();
            dem.calcTotalDeformation();
            dem.outputDeformationConf();
            dem.outputYaplotDEM();
            dem.outputLogDEM();
            dem.outputData();
            dem.outputConfiguration();
        }while (dem.time_simulation  < time_to_keep_shear);    
    }

	return;
}
