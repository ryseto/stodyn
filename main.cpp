#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "my_utilities.h"
#include "system.h"
#include "vec3d.h"
#include "dem.h"
#include "common.h"
#include "rigidClusterAnalysis.h"
#include "torqueBlanceMotion.h"
#include "calcDragAndTorque.h"
#include "testSimulation.h"
using namespace std;
int main (int argc, char** argv) {
	/*********************************************
     * 
     * 
     *
     *********************************************/
	if ( argc <= 1 ){
		cerr << "Usage: stodyn TYPE ..." << endl;
        cerr << "R: Rigid matrix and torque free motion" << endl;
        cerr << "D: DEM simulation" << endl;
		cerr << "U: calcDragAndTorque in uniform flows" << endl;
		cerr << "S: calcDragAndTorque in shear flows" << endl;
        cerr << "A: Max forces" << endl;
        cerr << "T: Test simulation" << endl;

        cerr << "--inactive development--" << endl;
		cerr << "c: compare_fem_and_sd" << endl;
		cerr << "t: testRyuon" << endl;
		cerr << "g: makeProfileData" << endl;
		cerr << "f: sedimentRigidAggregate" << endl;
		return 0;
	}

	/*********************************************/
    switch (argv[1][0]){
        case 'R':
            /* Torque balance motion
             * (force-free aggregates in shear flow)
             * for the paper.
             * "Hydrodynamic stress on small colloidal aggregates 
             *  in shear flow using Stokesian dynamics"
			 *
             */
            torqueBalancedMotion(argc, argv);
            break;
        case 'D':
            /* DEM simulation for an isolated cluster in shear flow.
             * The method is explained in the publification;
             * "Restructuring of colloidal aggregates in shear flow: 
             * Coupling interparticle contact models with Stokesian dynamics".
             */
           
            demSimulation(argc, argv);
            break;
        case 'U':
            /* Drag forces in uniform flow.
             */
            calcDragAndTorque(argc, argv);

            break;
        case 'S':
            /* Drag forces in shear flow.
             */
            calcDragAndTorque(argc, argv);

            break;
        case 'A':
            /*
             * This is for the paper:
             * "Hydrodynamic stress on small colloidal aggregates 
             *  in shear flow using Stokesian dynamics"
             */
            rigidClusterAnalysis(argv);
            break;
        case 'T':
            /*
             * Simple version of particle simulation.
             *  
             *
             */
            testSimulation(argc, argv);
            break;
        default:
            cerr << "R/D/U/S/I/A/" << endl;
    }
    cerr << "-" << endl;
	return 0;
}
