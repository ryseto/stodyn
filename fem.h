/*
 *  fem.h
 *  stodyn
 *
 *  Created by seto on 10/05/21.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include "system.h"

string trim_cstr(char *cstring, string from_trim, string to_trim){
	string tmp_string = cstring;
	cerr << tmp_string << endl;
	int i_backslash = tmp_string.find_last_of(from_trim.c_str()) + 1;
	int i_extention = tmp_string.find(to_trim.c_str());
	cerr << "trim " << i_backslash << ' ' << i_extention-i_backslash << endl;
	return tmp_string.substr( i_backslash, i_extention-i_backslash);
	
}

/* read fortran double precision number 
 * XD+Y (X = significand, Y = exponent)
 */
double abuf_to_double(char *moji){
	string s_moji = moji;
    int index = 0;
	double d_value;
	if ( s_moji.find("D") < s_moji.size() ) 	{
		int i = s_moji.find_last_of( "D") ;
		int j = s_moji.size();
		char exponent[5];
		char significand[32];
		sprintf(exponent, "%s", (s_moji.substr(i+1,j)).c_str());
		sprintf(significand, "%s", (s_moji.substr(0,i)).c_str());
		index = atoi(exponent);
		d_value = atof(significand);
	} 	
	return d_value*pow_10(index);
}

void import_FEM(System &sy,
				char *full_filename,
				int &xyz){
	ifstream fin;
	fin.open(full_filename);
	if (! fin.is_open() ){
		cerr << " no such file : " << full_filename << endl;
		exit(1);
	}
	if (fin.fail()) {
        std::cerr << "Error: Could not open\n";
        exit (8);
    }
	string filename;
	filename = trim_cstr(full_filename, "/", ".dat");
	sy.setFilename(filename);
	xyz = 0; // 1=x, 2=y, 3=z;
	if ( filename.find("uni") < filename.size() ) {
		sy.type_of_flow = 'u';		
		if (filename.find("uni01") < filename.size()  
			|| filename.find("_x") < filename.size() 
			|| filename.find("_x2") < filename.size()){
			xyz = 1;			
		} else if (filename.find("_y") < filename.size()
				   || filename.find("uni02") < filename.size()){
			xyz = 2;
		} else if ( filename.find("_z.") < filename.size()
				   || filename.find("uni03") < filename.size()
				   ){
			xyz = 3;
		} else {
			cerr << "file name uni0[123]" << endl;
			exit(1);
		}
	} else if ( filename.find("shear") < filename.size() ){
		sy.type_of_flow = 's';
		if ( filename.find("shear01") < filename.size() ){
			xyz = 1;			
		} else if ( filename.find("shear02") < filename.size()){
			xyz = 2;
		} else if ( filename.find("shear03") < filename.size()){
			xyz = 3;
		} else {
			cerr << "file name shear0[123]" << endl;
			exit(1);
		}
	} else {
		
		cerr << "file name is not correct" << endl;
	}
	char buf[1000];
	fin.getline( buf, 1000 );
	double pos_[3];
	double a;
	double area;
	double f_[3];
	double t_[3];
	double force_unit = 1e-9;
	double torque_unit = 1e-9*1e-6;
	while ( !fin.eof() ) {
		fin >> buf >> pos_[0] >> pos_[1] >> pos_[2] >> a ;
		fin >> buf; area  = abuf_to_double(buf); 
		fin >> buf; f_[0] = force_unit*abuf_to_double(buf);// [nN]
		fin >> buf; f_[1] = force_unit*abuf_to_double(buf);// [nN]
		fin >> buf; f_[2] = force_unit*abuf_to_double(buf);// [nN]
		fin >> buf; t_[0] = torque_unit*abuf_to_double(buf);// [nN.micro_m]
		fin >> buf; t_[1] = torque_unit*abuf_to_double(buf);// [nN.micro_m]
		fin >> buf; t_[2] = torque_unit*abuf_to_double(buf);// [nN.micro_m]
		sy.pos_fem.push_back( vec3d(pos_[0], pos_[1], pos_[2]));
		sy.f_fem.push_back(   vec3d(f_[0], f_[1], f_[2]));
		sy.t_fem.push_back(   vec3d(t_[0], t_[1], t_[2]));
	}
		
	sy.pos_fem.pop_back();
	sy.f_fem.pop_back();
	sy.t_fem.pop_back();
	
	
	
}

void compare_fem_and_sd(char *filename, System &sy)
{
	double L = 200;
	sy.setBox(L, L, L); 
	sy.L0_fem = 1.0e-6;// [m]
	sy.a_fem = 0.735e-6;//		
	//double eta_fem = 0.001;
	int xyz;
	import_FEM(sy, filename, xyz);
	vector< vec3d > pos;
	pos.resize( sy.pos_fem.size() );
	double length_unit_fem = sy.a_fem / sy.L0_fem ;
	for(int i=0; i< sy.pos_fem.size(); i++){
		pos[i] = sy.pos_fem[i] / length_unit_fem;
	}
		
	sy.setPosition(pos);
	sy.init();
	if ( sy.type_of_flow == 'u'){
		sy.U0_fem = 5e-5; // [m/s]
		sy.F0_fem = 6.0*M_PI*sy.eta*sy.U0_fem*sy.a_fem;
		vec3d U;
		switch (xyz) {
			case 1: U.set(-1,0,0); break;
			case 2:	U.set(0,-1,0); break;
			case 3:	U.set(0,0,-1); break;				
			default : cerr << "xyz = " << xyz << endl; exit(1);
		}
		sy.setImposedFlow_uniform(U);
	} else if (sy.type_of_flow == 's'){
		double shearrate = 50;
		//double shearrate_fem = 0.00005;
		sy.U0 = shearrate*sy.a;
		sy.F0 = 6.0*M_PI*sy.eta*sy.a*sy.U0;
		double shearrate_fem = 50;

		//double omega_fem = shearrate_fem / 2;
		//sy.U0_fem = shearrate_fem*(2*sy.a_fem);
		sy.U0_fem = shearrate_fem*sy.a_fem;
		sy.F0_fem = 6.0*M_PI*sy.eta*sy.a_fem*sy.U0_fem;
		//double T0_fem = 8.0*M_PI*sy.eta*sy.a_fem*sy.a_fem*sy.a_fem*omega_fem;
		sy.setImposedFlow_simpleshearflow(-1, xyz);
	} else {exit(1);}
	
	cerr << filename << endl;
	string tmp_string = sy.filename;	
	int i_backslash = tmp_string.find("00a") + 3;
	int i_extention = tmp_string.find("shear");
	cerr << "trim " << i_backslash << ' ' << i_extention-i_backslash << endl;
	string sample =  tmp_string.substr( i_backslash, i_extention-i_backslash);
	char out_configuration[64];
	sprintf(out_configuration, "o%d-%s.dat", sy.np, sample.c_str());
	ofstream fout_init(out_configuration);
	double sum_sq_dist =0;
	for(int i = 0 ; i< sy.np; i++){
		vec3d p(pos[i].x, pos[i].y, pos[i].z);
		sum_sq_dist += p.sq_norm();
	}
	//double Rg = sqrt(sum_sq_dist /sy.np );
	
	fout_init << "# N " << sy.np << endl;
	fout_init << "# a 1.0 "<< endl;
	fout_init << "# gap 0.00136 "<< endl;
	for(int i=0; i< sy.pos_fem.size(); i++){
		fout_init << pos[i].x << ' ' << pos[i].y  << ' ' << pos[i].z << endl;
	}

	
	
	cerr << "a_fem = " << sy.a_fem << endl;
	cerr << "U0_fem = " << sy.U0_fem << endl;
	cerr << "F0_fem = " << sy.F0_fem << endl;
	sy.setStokesianDynamicsCalc();
	vec3d U(0,0,0);
	sy.setAllParticleVelocity_to(U);
	
	/***************************************************************/
	/* main calculation                                            */
	//sy.solve();
	/***************************************************************/
	/***************************************************************/
	/* output data                                                 */
	//sy.calcFreeDrainingForce();
    sy.FreeDrainingApproximation();
	sy.openOutputFileStream();
	sy.outputYaplot();
	
	sy.FEM_normalization();
	//	for(int i=0; i< sy.np; i++){
	//		f_fem[i] = f_fem[i] / sy.F0_fem;
	//		t_fem[i] = t_fem[i] / (sy.F0_fem*length_unit_fem);
	//	}
	sy.outputYaplot_FD();
	sy.outputYaplot_fem();
	
	sy.outputAnalysis();
}
