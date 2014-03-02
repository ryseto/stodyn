//
//  torqueBlanceMotion.h
//  stodyn
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012年 __MyCompanyName__. All rights reserved.
//

#ifndef stodyn_torqueBlanceMotion_h
#define stodyn_torqueBlanceMotion_h
#include "system.h"
//#include "my_utilities.h"



/*
 * Find torqueBalancedMotion in shear flow by using the resistantce matrix for rigid body.
 *
 * INPUT:
 * arg1='R' : This analysis is selected already. 
 * arg2 : File of configuration of particles
 * arg3 : skip number of lines for the configuration file
 * arg4 : Method (0: FDA, 1: SD without lub, 2: SD with lub
 */ 
void torqueBalancedMotion(int argc, char** argv);

void export_profile(System &sy, bool export_GRM);



#endif
