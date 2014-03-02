//
//  demsystem.cpp
//  stodyn
//
//  Created by SETO Ryohei on 31/07/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#include "demsystem.h"
#include <iostream>
#include <map>

DEMsystem::DEMsystem(SDsystem &sd_sys_){
    sd_sys = &sd_sys_;
    ct = new ConTable;
    grid = new Grid;
    mov_allocated = false;

    sprintf(version, "0");
	n_bond = 0 ; /* number of bond */
    restructuring = true;
}

DEMsystem::~DEMsystem(){
    if (fout_dem.is_open())
        fout_dem.close();
    if (fout_dem_log.is_open())
        fout_dem_log.close();
    if (fout_conf.is_open())
        fout_conf.close();
    if (fout_conf_def.is_open())
        fout_conf_def.close();
    if (fout_trace.is_open())
        fout_trace.close();
    if (fout_data.is_open())
        fout_data.close();
    if (mov_allocated){
        DELETE(mov);
        DELETE(vos);
        DELETE(fte);
    }
    delete grid;
    delete ct;
};

void DEMsystem::setBox(double lx_, double ly_, double lz_){
    lx = lx_;
    ly = ly_;
    lz = lz_;    
    lx0 = lx/2;
    ly0 = ly/2;
    lz0 = lz/2;
}

void DEMsystem::setVersion(char *version_)
{
	sprintf(version, "%s", version_);
}

void DEMsystem::setInitialMotion(){
    /* (1)resting 
     * (2)torque balanced motion
     */
    for (int i=0; i< np; i++){
        particle[i]->velocity.set(0,0,0);
        particle[i]->omega.set(0,0,0);
    }
}

void DEMsystem::calcInterParticleForce(){
    ForAllParticle{ 
        (*p_iter)->resetForce();
    }
    ForAllBond{
        if( (*bond_iter)->status ){
            (*bond_iter)->addContactForce();
        }
    }
}

void DEMsystem::calcVelocityOmega(){
    /* velocity and angular velocity vectors are on the rotating coordinate
     */
    for (int i=0; i < np; i++){
        for (int k=0; k < 6 ; k++){
            int ii = i*11 + k;
            vos[ii] = 0;
            for (int j = 0 ; j < n11 ; j++){
                vos[ii] += mov[ii*n11 + j]*fte[j];
            }
        }
    }
}

void DEMsystem::calcStresslet(){
    /* The stresslets are on the rotating cooridnate
     */
    for (int i=0; i < np ; i++){
        for (int k=6; k < 11 ; k++){
            int ii = i*11 + k;
            vos[ii] = 0;
            for (int j = 0 ; j < n11 ; j++){
                vos[ii] += mov[ii*n11 + j]*fte[j];
            }
        }
    }
}

void DEMsystem::freeDrainingApproximation(){
    count_MD_steps ++;
    time_simulation += timestep;
    double shear_vs_bond_inv = 1 / shear_vs_bond;
    for (int i = 0; i < np ; i ++){
        f_tmp = shear_vs_bond_inv * particle[i]->force;
        t_tmp = shear_vs_bond_inv * particle[i]->torque;
        
        particle[i]->velocity.set(f_tmp.x + particle[i]->p.z - lz0, f_tmp.y, f_tmp.z);
        
        particle[i]->omega.set(t_tmp.x, t_tmp.y + 0.5, t_tmp.z);
    }
    
    ForAllParticle{ 
        (*p_iter)->move_with_velocity();
    }
    
}


void DEMsystem::calcVelocityByMov(){
    count_MD_steps ++;
    time_simulation += timestep;
    /* Stokesian dynamics */
    //////////////////////
    /* for inv_rotation_quickcalc_r and rotation_quickcalc
     */
    q_rot.preparer_rotation_quickcalc(); 
    double tmp_13 = q_rot.q.x * q_rot.q.z;
    double tmp_02 = q_rot.q0  * q_rot.q.y;
    double tmp_00_22 = q_rot.q0  * q_rot.q0 - q_rot.q.y * q_rot.q.y;
    double tmp_11_33 = q_rot.q.x * q_rot.q.x - q_rot.q.z * q_rot.q.z;
    /*
     * E' = B.E.A
     */
    double rot_0  = tmp_00_22 - tmp_11_33;
    double rot_2  = 2*(tmp_13 - tmp_02);
    double rot_3  = 2*(q_rot.q.x*q_rot.q.y - q_rot.q0*q_rot.q.z);
    double rot_5  = 2*(q_rot.q.y*q_rot.q.z + q_rot.q0*q_rot.q.x);
    double irot_2 = 2*(tmp_13 + tmp_02);
    double irot_8 = tmp_00_22 - tmp_11_33;
    e[0] = rot_0*rot_2;
    e[1] = 0.5*(rot_0* rot_5 + rot_2* rot_3);
    e[2] = 0.5*(rot_0*irot_8 + rot_2*irot_2);
    e[3] = 0.5*(rot_3*irot_8 + rot_5*irot_2);
    e[4] = rot_3*rot_5;
    double shear_vs_bond_inv = 1 / shear_vs_bond;
    for (int i = 0; i < np ; i ++){
        int i11 = i*11;
        f_tmp = q_rot.inv_rotation_quickcalc_r( particle[i]->force );
        t_tmp = q_rot.inv_rotation_quickcalc_r( particle[i]->torque );
        //////////////////////////////////////////
        //// 
        fte[i11]   = shear_vs_bond_inv*f_tmp.x;
        fte[i11+1] = shear_vs_bond_inv*f_tmp.y;
        fte[i11+2] = shear_vs_bond_inv*f_tmp.z;
        fte[i11+3] = shear_vs_bond_inv*t_tmp.x;
        fte[i11+4] = shear_vs_bond_inv*t_tmp.y;
        fte[i11+5] = shear_vs_bond_inv*t_tmp.z;
        //////////////////////////////////////////
        fte[i11+6]  = e[0];
        fte[i11+7]  = e[1];
        fte[i11+8]  = e[2];
        fte[i11+9]  = e[3];
        fte[i11+10] = e[4];    
    }
    calcVelocityOmega();
 
    for (int i=0; i < np ; i++){
        int i11= i*11;
        v.set(vos[i11],vos[i11+1], vos[i11+2]);
        o.set(vos[i11+3], vos[i11+4], vos[i11+5]);
        q_rot.rotation_quickcalc( v );
        q_rot.rotation_quickcalc( o );
        /* 
         * adding imposed flow u_inf and o_inf 
         */
        particle[i]->velocity.set(v.x + particle[i]->p.z - lz0, v.y, v.z);
        particle[i]->omega.set(o.x, o.y + 0.5, o.z);
    }
    ForAllParticle{ 
        (*p_iter)->move_with_velocity();
    }
}

/* Positions of particles are set to SD_system
 */
void DEMsystem::pos_from_DEM_to_SD(){
    for( int i=0; i< np; i++){
        sd_sys->setPosition(i, particle[i]->p);
    }        
}  

void DEMsystem::set_FDA(){
    double self_VF = 1.0;
    double self_OT = 0.75;
    double self_ES = 0.9;
    for (int i = 0; i < n11*n11; i++){
        mov[i] = 0.0;
    }
    
    for (int i = 0; i < np; i++){
        int i11 = i*11;
        for (int k = 0; k < 3; k++){
            mov[n11*(i11+k)   + i11 + k] = self_VF ;
            mov[n11*(i11+3+k) + i11 + 3+k] = self_OT ;
        }
        for (int k = 0; k < 5; k++){
            mov[n11*(i11+6+k) + i11 + 6+k] = 1.0/self_ES ;
        }
    }
}

void DEMsystem::getMovMatrix(){
    // position transport from DEMsystem to SDsystem.
    count_SD_calc ++;
    pos_from_DEM_to_SD();
    sd_sys->calcGrandMovMatrix(mov);
}

/* Import the configuration (x,y,z) of particles.
 * The center of mass is set to (0,0,0)
 * INPUT 
 * importfilename: File name
 * skipline: Line number of the header
 */
void DEMsystem::importCluster(char* importfilename, int skipline){
    vector<vec3d> init_config;

	sprintf(init_config_file, "%s", importfilename);
	string s_filename = init_config_file;
	int i_backslash = s_filename.find_last_of( "/") + 1;
	int i_extention = s_filename.find( ".dat" );	
	sprintf(init_config_file, "%s",
			(s_filename.substr(i_backslash,i_extention-i_backslash)).c_str());
	ifstream fin;
    cerr << init_config_file << endl;
	fin.open( importfilename );
	double x, y, z;
	char buf[1000];
	for (int i = 0; i< skipline; i++){
		fin.getline( buf, 1000 );
	}
	do{
		fin >> x >> y >> z;
        if(fin.fail())
            break; 
        vec3d new_p(x, y, z);
		init_config.push_back(new_p);
	} while (!fin.eof());
	shiftCenterOfMass( init_config );
    np = init_config.size();
    n11 = 11*np;
    
    pos_init = new vec3d [np];

    double sum_sq_dist = 0.0;	
    for (int i = 0; i< np ; i++){
        pos_init[i] = init_config[i];
        sum_sq_dist += init_config[i].sq_norm();
    }	
    radius_of_gyration = sqrt(sum_sq_dist/np);    
    init_radius_of_gyration = radius_of_gyration;
}




void DEMsystem::shiftCenterOfMass(vector<vec3d> &p){
	vec3d cm(0,0,0);
	foreach( vector<vec3d>, p, p_iter){
		cm += (*p_iter);
	}
	cm *= 1.0/p.size();
	foreach( vector<vec3d>, p, p_iter){
		*p_iter -= cm;
	}
}

void prepareBond(BondParameter & _bond){
    double a0 = 1.0;
    _bond.kn = _bond.fnc / _bond.n_max;
    _bond.ks = _bond.fsc / _bond.s_max;
    _bond.kb = _bond.mbc / (_bond.b_max * a0 * a0);
    _bond.kt = _bond.mtc / (_bond.t_max * a0 * a0);
    _bond.kn3 = _bond.overlap_force_factor * _bond.kn / sq(_bond.n_max);
    
    _bond.c_norm = 2*sqrt( _bond.kn);
    _bond.c_slid = 2*sqrt( _bond.ks);
    _bond.c_bend = 2*sqrt( _bond.kb * 1.0 * 1.0);
    _bond.c_tort = 2*sqrt( _bond.kt * 1.0 * 1.0);
}

void setBondParameter( BondParameter &bondparameter, string &bond_file, string &dir){
    string path_bond = dir + "/" + bond_file;
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
    fin.close();
    prepareBond(bondparameter);
    return;
}

void DEMsystem::readBondParameter(){
    setBondParameter(bond0, bond0_file, rootdir);
    setBondParameter(bond1, bond1_file, rootdir);
}

void DEMsystem::readParameterKey(const string &codeword,                               
                                 const string &value){
	map<string,int> keylist;
    keylist["rootdir:"]=1; const int _rootdir = 1;
    keylist["bond0_file:"]=2; const int _bond0_file = 2;
    keylist["bond1_file:"]=3; const int _bond1_file = 3;
    keylist["method_hydroint:"]=4; const int _method_hydroint = 4;
    keylist["shear_vs_bond_min:"]=5; const int _shear_vs_bond_min = 5;
    keylist["shear_vs_bond_max:"]=6; const int _shear_vs_bond_max = 6;
    keylist["stepnumber_shear_vs_bond:"]=7;const int _stepnumber_shear_vs_bond = 7;
    keylist["time_each_shear_vs_bond:"]=8; const int _time_each_shear_vs_bond = 8;
    
    keylist["time_interval_output:"]=9; const int _time_interval_output = 9;
    keylist["critical_deformation_SD:"]=10; const int _critical_deformation_SD = 10;
    keylist["L_dem:"]= 30; const int _L_dem = 30;
    keylist["L_sd:"]= 31; const int _L_sd = 31;
	cerr << codeword << ' ' << value  << endl;
	switch(keylist[codeword]){
        case _rootdir: rootdir = value; break;
        case _bond0_file: bond0_file = value; break;
        case _bond1_file: bond1_file = value; break;
        case _method_hydroint: sd_sys->method_hydroint = atoi(value.c_str()); break;
        case _shear_vs_bond_min: shear_vs_bond_min =  atof(value.c_str() ); break;
        case _shear_vs_bond_max: shear_vs_bond_max =  atof(value.c_str() ); break;
//        case _time_interval: time_interval = atof(value.c_str()); break;
        case _stepnumber_shear_vs_bond: stepnumber_shear_vs_bond = atoi(value.c_str()); break;
        case _time_each_shear_vs_bond: time_each_shear_vs_bond = atof(value.c_str() ); break; 
        case _time_interval_output: time_interval_output = atof(value.c_str());
        case _L_dem: L_dem = atof(value.c_str() ); break;
        case _L_sd: L_sd = atof(value.c_str() ); break;
        case _critical_deformation_SD: critical_deformation_SD = atof(value.c_str()); break;
		default:
			cerr << "The codeword " << codeword << " is'nt associated with an parameter" << endl;
			exit(1);
	}
}

void DEMsystem::readParameterFileDEM(){
	//////////////////////////////////////////	
	// input arguments:
	//////////////////////////////////////////	
	// Read parameter file
	ifstream fin;
	fin.open(parameters_file);
	string codeword, value, unit;
	while (!fin.eof()){
		fin >> codeword ;
		if ( codeword == "#") {
			char buf[1024]; fin.get(buf, 1024);
		} else if (codeword == "!"){
			break;
		} else {
			fin >> value; if (fin.eof()) break;
			//set_value(codeword, value, unit);
			readParameterKey(codeword, value);
		}
	}
    cerr << "completed to read parameters." << endl;
}

void DEMsystem::setSimulationBoxs(){
    sd_sys->setBox(L_sd, L_sd, L_sd);
    setBox(L_dem, L_dem, L_dem);
}

void DEMsystem::setTimeStep(){
    if ( sd_sys->method_hydroint == 0){
        timestep = 0.0001*shear_vs_bond;
        if ( timestep > 0.0001)
            timestep = 0.0001;
    } else {
        timestep = 0.0004*shear_vs_bond;
        if ( timestep > 0.0002)
            timestep = 0.0002;
    }
    cerr << "dt = " << timestep << endl;
}

void DEMsystem::initDEM(){
	for (int i=0; i < np ; i++){
		particle.push_back(new Particle(i, *this, *sd_sys) );
	}
	ct->set(np);	
	initGrid();
    
    string str_version = version;
    if (str_version == "nr"){
        // no restructuring
        restructuring = false;
    }
   
    time_simulation = 0.;
    init_continuity = 1;
//    rup_normal = 0;
//    rup_shear = 0;
//    rup_bend = 0;
//    rup_torsion = 0;
	rupture_bond.clear();
	regeneration_bond.clear();

	counterBreak = 0;
	counterRegenerate = 0;
    count_MD_steps = 0;
    count_SD_calc = 0;
    vec3d origin_shift(lx0,ly0,lz0);

    for (int i=0; i < np ; i++){
		particle[i]->p = pos_init[i] + origin_shift;
	}
    double contact_distance = 2.001;
    if (restructuring == false){
        contact_distance = 2.1;
    }
        
	makeInitialBond(contact_distance);
    mov_allocated = true;
    mov = new double [np*np*121];
    vos = new double [np*11];
    fte = new double [np*11];
    pos = new vec3d [np];

    set_FDA();
    pos_center_of_mass.set(lx0, ly0, lz0);
    
    q_rot.reset();
    q_rot_total.reset();
//    delta_grad_method = 1e-3;
}


void DEMsystem::makeInitialBond(double generation_distance)
{
	double tmp = sq_dist_generate;
	sq_dist_generate = sq(generation_distance);
	makeNeighbor();
	generateBond();	
	sq_dist_generate = tmp;
}

void DEMsystem::makeNeighbor(){
	grid->remake(particle);
	ForAllParticle{ (*p_iter)->makeNeighbor(); }
}

void DEMsystem::generateBond(){
	for (int i=0; i < particle.size(); i++){
		particle[i]->generateBond();
	}
}

void DEMsystem::resetDeformation(){
    q_rot_total_step = q_rot*q_rot_total_step;
    q_rot.reset();
}

void DEMsystem::outputDeformationConf(){
    static bool firsttime = true;
	if ( firsttime ){firsttime = false;}else{fout_conf_def << endl;}

    ForAllParticle{
        vec3d p((*p_iter)->p.x - lx0,
                (*p_iter)->p.y - ly0,
                (*p_iter)->p.z - lz0);
        vec3d pp = q_rot_total.inv_rotation(p);
        fout_conf_def << "c ";
        fout_conf_def << pp.x << ' ';
        fout_conf_def << pp.y << ' ';
        fout_conf_def << pp.z << endl;
    }
}

double DEMsystem::evaluateObjFunction(quaternion & q_, vec3d *po){
    double Obj_tmp = 0;
    q_.preparer_rotation_quickcalc();
    for (int i = 0; i < np; i++){
        Obj_tmp += (q_.rotation_quickcalc_r(po[i]) - pos[i]).sq_norm();
    }
    return Obj_tmp / np;
}

double DEMsystem::findObjMinimum(quaternion &q_try, vec3d *po,
                                 double delta_init, double conv_diff){                                 
    double delta_grad_method = delta_init;
    double gradObj[4];
    /* pos[i  ] is the precent position of particles*/
    for (int i = 0; i< np; i++){
        pos[i].set(particle[i]->p.x -lx0,
                   particle[i]->p.y -ly0,
                   particle[i]->p.z -lz0);
    }
    double Obj;
    double Obj_old = evaluateObjFunction(q_try, po);
    vec3d d_pos;
    cnt_grad_method = 0;
    vec3d tmp1;
    vec3d tmp2;
    vec3d tmp3;
    quaternion q_tmp;    
    do{
        gradObj[0] = gradObj[1] = gradObj[2] = gradObj[3] = 0;
        q_try.preparer_rotation_quickcalc();
        for (int i = 0; i< np; i++){
            /* po[i]    :  master position
             * rot_pos[i]   :  \bm{r}'_i
             * pos[i]       :  \bm{r}_i current position
             */
            d_pos = q_try.rotation_quickcalc_r(po[i]) - pos[i];
            double dot_r_q = dot(po[i], q_try.q);
            tmp1.set( dot_r_q         , -q_try.q0*po[i].z,  q_try.q0*po[i].y);
            tmp2.set( q_try.q0*po[i].z,  dot_r_q         , -q_try.q0*po[i].x);
            tmp3.set(-q_try.q0*po[i].y,  q_try.q0*po[i].x,  dot_r_q         );
            gradObj[0] += 4*dot(d_pos,  q_try.q0 *po[i] - cross(po[i] , q_try.q));            
            gradObj[1] += 4*dot(d_pos, -q_try.q.x*po[i] + po[i].x*q_try.q + tmp1);
            gradObj[2] += 4*dot(d_pos, -q_try.q.y*po[i] + po[i].y*q_try.q + tmp2);
            gradObj[3] += 4*dot(d_pos, -q_try.q.z*po[i] + po[i].z*q_try.q + tmp3);
        }
        
        while (true){
            q_tmp.q0  = q_try.q0  - delta_grad_method*gradObj[0];
            q_tmp.q.x = q_try.q.x - delta_grad_method*gradObj[1];
            q_tmp.q.y = q_try.q.y - delta_grad_method*gradObj[2];
            q_tmp.q.z = q_try.q.z - delta_grad_method*gradObj[3];
            Obj = evaluateObjFunction(q_tmp, po);
            if (Obj <= Obj_old ){
                q_try = q_tmp;
                break;
            }
            delta_grad_method *= 0.5;
        }
        if (abs(Obj - Obj_old) < conv_diff){
            break;
        }
        delta_grad_method *= 1.1;
        Obj_old = Obj;
        cnt_grad_method ++;
    } while (true);
    q_try.normalize();
    //Obj = evaluateObjFunction(q_try, po);
    return sqrt(Obj);
}

void DEMsystem::estimateClusterRotation(){
    step_deformation = findObjMinimum(q_rot, sd_sys->dr, 1e-4, 1e-12);
}

/*
 * DEM simulation
 */ 

void DEMsystem::setParameterFileDEM(char *parameters_file_){
	sprintf(parameters_file, "%s", parameters_file_);
	string s_parameters_file = parameters_file;
	int i_backslash = s_parameters_file.find_last_of( "/") + 1;
	int i_extention = s_parameters_file.find( ".txt" );	
	sprintf(parameters, "%s",
			(s_parameters_file.substr( i_backslash, i_extention-i_backslash)).c_str());
	cerr << " parameters = " << parameters << endl;
}

//
//void DEMsystem::setParticlePosition_sdpos_to_dem(){
//	for (int i=0; i < np; i++){
//        sd_sys->copyPosition( particle[i]->p ,i);		
//	}
//}

void DEMsystem::initGrid(){
	if (np == 0)
	{
		cerr << "before setting system box, ";
		cerr << "it is required to import configuration of particles.";
		cerr << endl;
		exit(1); 
	}
	double cell_h = 2.2;
    double l_box[3]={lx, ly, lz};
	grid->init(np, l_box, cell_h);
}

void DEMsystem::checkFailure(){
    ForAllBond{
        if ( (*bond_iter)->status ){
            (*bond_iter)->failureCondition();
        }
    }
}

void DEMsystem::regeneration_onebyone(){

    double D_max = 0.;
    int most_stressed_bond = -1;
    int n_regeneration_bond = regeneration_bond.size();
    for (int i=0; i < n_regeneration_bond; i++){
        if ( D_max < bond[ regeneration_bond[ i ] ]->D_function ){
            D_max =  bond[ regeneration_bond[ i ] ]->D_function;
            most_stressed_bond = regeneration_bond[ i ];
        }
    }
    bond[ most_stressed_bond  ]->regeneration();
    counterRegenerate ++;
}

void DEMsystem::rupture_onebyone(){
    double D_max = 0;
    int most_stressed_bond = -1;
    int n_rupture_bond =  rupture_bond.size();
    for (int i=0; i< n_rupture_bond; i++){
        if ( D_max < bond[ rupture_bond[ i ] ]->D_function ){
            D_max =  bond[ rupture_bond[ i ] ]->D_function;
            most_stressed_bond = rupture_bond[ i ];
        }
    }
    bond[ most_stressed_bond ]->rupture();
    counterBreak ++;
}

void DEMsystem::calcTotalFTS(){
    /*!!!!!!!!!!!!!!!!!!!!!!
     * This calculation can be totally wrong.
     * This is due to the confusion about dimensionless handling.
     * The total stresslet need to be build from
     * the output data.
     *!!!!!!!!!!!!!!!!!!!!!!
     */
    calcStresslet();

    vec3d f, t;
    double s[5];   // stress[ xx, xy, xz, yz, yy ];
    total_torque.reset();
	total_force.reset();    
    for (int i=0; i< 5; i++){
        total_stresslet[i] = 0.0;
    }

    double s_average[5];
    double sp_total[5];
    double so_total[5];
    for( int k=0; k < 5; k++){
        s_average[k] = 0;
        sp_total[k] = 0;
        so_total[k] = 0;
    }

    //    cout << "r 0.1" << endl;
	for (int i = 0; i < np; i ++){
		f = - particle[i]->force;
		t = - particle[i]->torque;
        //      cout << "c " ; 
        //        cout << pos[i].x << ' ' << pos[i].y << ' ' << pos[i].z << endl;
        //        cout << "l ";
        //        cout << pos[i].x << ' ' << pos[i].y << ' ' << pos[i].z << ' ';
        //        cout << pos[i].x + 10*f.x << ' ' << pos[i].y + 10*f.y << ' ' << pos[i].z + 10*f.z << endl;
        
        for (int k=0; k < 5 ; k++){
            s[k] = vos[i*11 + k + 6];
            s_average[k] += s[k];
        }
		total_force += f;
        total_torque += (cross(pos[i],f) + t);
        sp_total[0] += s[0];
        sp_total[1] += s[1];
        sp_total[2] += s[2];
        sp_total[3] += s[3];
        sp_total[4] += s[4];
        so_total[0] += (1./3)*(2*pos[i].x*f.x - pos[i].y*f.y - pos[i].z*f.z);
        so_total[1] += (1./2)*(  pos[i].y*f.x + pos[i].x*f.y);
        so_total[2] += (1./2)*(  pos[i].z*f.x + pos[i].x*f.z);
        so_total[3] += (1./2)*(  pos[i].z*f.y + pos[i].y*f.z);
        so_total[4] += (1./3)*(2*pos[i].y*f.y - pos[i].x*f.x - pos[i].z*f.z);
	}
    for (int k=0; k < 5 ; k++){
        total_stresslet[k] = sp_total[k] + so_total[k];
    }
    norm_total_stresslet = sqrt(0.5*(sq(total_stresslet[0])
                                     +sq(total_stresslet[4])
                                     +sq(total_stresslet[0] +total_stresslet[4]))
                                + sq(total_stresslet[1])+sq(total_stresslet[2])+sq(total_stresslet[3]));
//    cout << endl;
//    cerr << sp_total[0]  << ' ';
//    cerr << sp_total[1]  << ' ';
//    cerr << sp_total[2]  << ' ';
//    cerr << sp_total[3]  << ' ';
//    cerr << sp_total[4]  << endl;
//    cerr << so_total[0]  << ' ';
//    cerr << so_total[1]  << ' ';
//    cerr << so_total[2]  << ' ';
//    cerr << so_total[3]  << ' ';
//    cerr << so_total[4]  << endl;
//    total_torque.cerr();
//    cerr << endl;

    return;
}


void DEMsystem::shiftClusterToCenter(){
    vec3d center_of_mass;
    center_of_mass.reset();
    for (int i=0; i < np ; i++){
        center_of_mass += particle[i]->p;
    }
    center_of_mass = center_of_mass/np;
    center_of_mass.x -= lx0;
    center_of_mass.y -= ly0;
    center_of_mass.z -= lz0;
    for (int i=0; i < np ; i++){
        particle[i]->p.x -= center_of_mass.x ;
        particle[i]->p.y -= center_of_mass.y;
        particle[i]->p.z -= center_of_mass.z;
    }
    pos_center_of_mass.add( -center_of_mass.x, -center_of_mass.y, -center_of_mass.z);
}

void DEMsystem::setDEMParameters(){
    double r_exponent = 1.0 /(stepnumber_shear_vs_bond - 1);
    ratio_shear_vs_bond = pow(shear_vs_bond_max/shear_vs_bond_min , r_exponent);
    dist_generate = 2.0;
	sq_dist_generate = dist_generate*dist_generate;	/* square of distance to generate new bond  */
}

void DEMsystem::openOutputFileDEM(){
	char filenameYaplotDEM[128];
    char filename_conf[128];
	char filenameLog[128];
    char filename_trace[128];
    char filename_conf_def[128];
    char filename_data[128];
    //	char filenameSample[64];
	string s_dem_algorithm;

    switch (sd_sys->method_hydroint){
        case 0:
            s_dem_algorithm = "FD";
            break;
        case 1:
            s_dem_algorithm = "SD_nolub";
            break;
        case 2:
            s_dem_algorithm = "SD_lub";
            break;
    }

    sprintf(filenameYaplotDEM, 
            "DEM_%s_%s_%s_%s.yap", 
            s_dem_algorithm.c_str(),
            parameters, 
            init_config_file, version);
    
	sprintf(filenameLog, "log_DEM_%s_%s_%s_%s.dat", s_dem_algorithm.c_str(),parameters, 
            init_config_file, version);
    
    sprintf(filename_conf,"conf_DEM_%s_%s_%s_%s.dat", 
            s_dem_algorithm.c_str(),parameters, 
            init_config_file, version);

    sprintf(filename_conf_def,"def_DEM_%s_%s_%s_%s.dat", 
            s_dem_algorithm.c_str(),parameters, 
            init_config_file, version);
    
    sprintf(filename_trace,"trace_DEM_%s_%s_%s_%s.dat", 
            s_dem_algorithm.c_str(),parameters, 
            init_config_file, version);
    
    sprintf(filename_data,"data_DEM_%s_%s_%s_%s.dat", 
            s_dem_algorithm.c_str(),parameters, 
            init_config_file, version);
    
	fout_dem.open(filenameYaplotDEM);
	fout_dem_log.open(filenameLog);
	fout_conf.open( filename_conf);
    fout_conf_def.open( filename_conf_def);
    fout_trace.open ( filename_trace);
    fout_data.open( filename_data);
    
	writeHeader(fout_dem);
	writeHeader(fout_dem_log);
}


void DEMsystem::calcLocalStrains(){
    /* The minimum and maximum distances between two particles.
     */
    r_min = 2.;
    r_max = 2.;
    bending_angle_max = 0.;
    torsional_angle_max = 0.;
    sliding_disp_max = 0.;
    moment_bend_max = 0.;
    moment_tors_max = 0.;
    force_trac_max = 0.;
    force_comp_max = 0.;
    force_slid_max = 0.;
    
    for (int i=0; i < bond.size(); i++){
        if (bond[i]->status){
            if (bond[i]->r < r_min){
                r_min = bond[i]->r;
            }
            if (bond[i]->r > r_max){
                r_max = bond[i]->r;
            }
            if (bond[i]->valSlidingDisplacement() > sliding_disp_max ){
                sliding_disp_max = bond[i]->valSlidingDisplacement(); 
            }                
            if (bond[i]->valBendingAngle() > bending_angle_max){
                bending_angle_max = bond[i]->valBendingAngle();
            }
            if (bond[i]->valTorsionalAngle() > torsional_angle_max){
                torsional_angle_max = bond[i]->valTorsionalAngle();
            }
            // force_normal = para.kn*q;
            if (bond[i]->val_F_norm() > 0.){
                if ( bond[i]->val_F_norm() > force_trac_max ){
                    force_trac_max = bond[i]->val_F_norm();
                }
            } else {
                if ( - bond[i]->val_F_norm() > force_comp_max ){
                    force_comp_max = - bond[i]->val_F_norm();
                }
            }
            if ( bond[i]->val_F_slid() > force_slid_max){
                force_slid_max = bond[i]->val_F_slid();
            } 
            if ( bond[i]->val_M_bend() > moment_bend_max){
                moment_bend_max = bond[i]->val_M_bend();
            }
            if ( bond[i]->val_M_tors() > moment_tors_max){
                moment_tors_max = bond[i]->val_M_tors();
            }
        }
    }
    
    if (r_max > 4){
        cerr << "r_max > 4 " << endl;
        exit(1);
    }
    
    return;    
}

void DEMsystem::calcGyrationRadius(){
    /* calculate gyration radius */
    double sum_sq_dist = 0.0;
    for (int i = 0; i< np ; i++){
        pos[i].set(particle[i]->p.x -lx0,
                   particle[i]->p.y -ly0,
                   particle[i]->p.z -lz0);
        sum_sq_dist += pos[i].sq_norm();
    }	
    radius_of_gyration = sqrt(sum_sq_dist/np);
}

void DEMsystem::calcTotalDeformation(){
    q_rot_total = q_rot*q_rot_total_step;
    total_deformation = findObjMinimum(q_rot_total, pos_init, 1e-4, 1e-12);

    if (init_continuity == 1){
        if (cnt_grad_method > 500)
            init_continuity = 0;
    }
}

void DEMsystem::checkState(){
//	for (int i=0; i < n_bond; i++){
//		if (bond[i]->status){
//			bond[i]->addContactForce();
//		}
//	}
    int number_of_bond = n_bond - counterBreak;
	double bforce_max = 0;
	double bforce_ave = 0;
	for (int i=0; i < n_bond; i++){
		if (bond[i]->status){
			double bforce = bond[i]->strenghOfForce();
			bforce_ave += bforce;
			if ( bforce > bforce_max ){
				bforce_max = bforce ;
			}
		}
	}
	ave_bondforce = bforce_ave / number_of_bond ;
	ave_force = 0.;
	max_force = 0.;
	max_velocity = 0.;
	max_ang_velocity = 0.;

    max_x = -1000;
    min_x = 1000;
    
	for (int i=0; i < np; i++){
		double force = particle[i]->valForce();
		double velocity = particle[i]->valVelocity();
		double ang_velocity = particle[i]->valOmega();
		ave_force += force;
		if ( force > max_force )
			max_force = force ;
		if ( velocity > max_velocity)
			max_velocity = velocity;
		if ( ang_velocity > max_ang_velocity )
			max_ang_velocity = ang_velocity;
        if ( max_x < particle[i]->p.x ) max_x = particle[i]->p.x;
        if ( min_x > particle[i]->p.x ) min_x = particle[i]->p.x;
	}
	ave_force = ave_force / np;
    calcLocalStrains();
//
//	for (int i=0; i < np; i++){
//		particle[i]->resetForce();
//	}
//    
    
}


void DEMsystem::writeHeader(ofstream &out){
	time_t tt= time(0);
    tm *time_now = localtime(&tt);
	
	out << "# parameters: " << parameters << endl;
	out << "# initial_cluster: " << init_config_file << endl;
	out << "# time: "  << time_now->tm_year + 1900 \
	<< ":" << time_now->tm_mon + 1\
	<< ":" << time_now->tm_mday  \
	<< ":" << time_now->tm_hour  \
	<< ":" << time_now->tm_min  \
	<< ":" << time_now->tm_sec << endl;
}

void DEMsystem::outputLogDEM(){
    static bool firsttime = true;
	if ( firsttime ){
		firsttime = false;
		fout_dem_log << "#1 time_simulation\n";
        fout_dem_log << "#2 shear_vs_bond\n";
        fout_dem_log << "#3 radius_of_gyration\n";
        fout_dem_log << "#4 relative_radius_of_gyration\n";
        fout_dem_log << "#5 total_deformation\n";
        fout_dem_log << "#6 step_deformation\n";
        fout_dem_log << "#7 init_continuity\n";
		fout_dem_log << "#8 average_contact_number\n";        
        fout_dem_log << "#9 bond number\n";
        fout_dem_log << "#10 counterRegenerate\n";
		fout_dem_log << "#11 counterBreak\n";
        fout_dem_log << "#12 total_force\n"; 
        fout_dem_log << "#13 total_torque\n"; 
        fout_dem_log << "#14 norm_total_stresslet\n"; 
		fout_dem_log << "#15 force_trac_max\n";
		fout_dem_log << "#16 force_comp_max\n";
        fout_dem_log << "#17 force_slid_max\n";
		fout_dem_log << "#18 moment_bend_max\n";			
		fout_dem_log << "#19 moment_tors_max\n";
        fout_dem_log << "#20 r_min\n";
        fout_dem_log << "#21 r_max\n";
        fout_dem_log << "#22 count_MD_steps\n";
        fout_dem_log << "#23 count_SD_calc\n";
        fout_dem_log << "#24 total_stresslet0\n";
        fout_dem_log << "#25 total_stresslet1\n";
        fout_dem_log << "#26 total_stresslet2\n";
        fout_dem_log << "#27 total_stresslet3\n";
        fout_dem_log << "#28 total_stresslet4\n";
        fout_dem_log << "#29 min_x\n";
        fout_dem_log << "#30 max_x\n";            
	}    
	/*1*/
    int number_of_bond = n_bond - counterBreak;
    double average_contact_number = 2.0*number_of_bond / np;

	fout_dem_log << time_simulation << ' '; // 1
    fout_dem_log << shear_vs_bond << ' '; //2
    fout_dem_log << radius_of_gyration << ' ' ; // 3
    fout_dem_log << radius_of_gyration/init_radius_of_gyration << ' '; // 4
    fout_dem_log << total_deformation << ' '; //5
    fout_dem_log << step_deformation << ' '; //6
    fout_dem_log << init_continuity << ' ' ; //7
    fout_dem_log << average_contact_number << ' '; //8
    fout_dem_log << number_of_bond << ' '; //9
    fout_dem_log << counterRegenerate << ' '; //10
    fout_dem_log << counterBreak << ' '; //11
    fout_dem_log << total_force.norm() << ' ';//12
    fout_dem_log << total_torque.norm() << ' ';//13
    fout_dem_log << norm_total_stresslet << ' ';//14
    fout_dem_log << force_trac_max << ' '; // 15
    fout_dem_log << force_comp_max << ' '; // 16
    fout_dem_log << force_slid_max << ' '; // 17
    fout_dem_log << moment_bend_max << ' '; // 18
    fout_dem_log << moment_tors_max << ' '; // 19
    fout_dem_log << r_min << ' '; //20
    fout_dem_log << r_max << ' '; //21
    fout_dem_log << count_MD_steps << ' '; // 22
    fout_dem_log << count_SD_calc << ' '; // 23
    fout_dem_log << total_stresslet[0] << ' ';//24
    fout_dem_log << total_stresslet[1] << ' ';//25
    fout_dem_log << total_stresslet[2] << ' ';//26
    fout_dem_log << total_stresslet[3] << ' ';//27
    fout_dem_log << total_stresslet[4] << ' ';//28
    fout_dem_log << min_x - lx0 << ' '; //29
    fout_dem_log << max_x - lx0 << ' '; //30
	fout_dem_log << endl; 

    fout_trace << time_simulation << ' ' ; // 1
    fout_trace << cnt_grad_method << ' ' ; // 2
    fout_trace << step_deformation << ' ' ; // 3
//    fout_trace << delta_grad_method << ' ';
    fout_trace << q_rot_total.q0 << ' '; //4
    fout_trace << q_rot_total.q.x << ' ';//5
    fout_trace << q_rot_total.q.y << ' ';//6
    fout_trace << q_rot_total.q.z << ' ';//7
    fout_trace << q_rot_total.norm() << endl;
}

void DEMsystem::outputConfiguration(){
    int number_of_live_bonds = n_bond - counterBreak ;
	fout_conf << "# ";
	fout_conf << time_simulation << ' ' << number_of_live_bonds;
	//fout_conf << ' ' << average_contact_number;
	fout_conf << endl;
    fout_conf << "P " << particle.size() << endl;
	ForAllParticle{
		fout_conf << (*p_iter)->p.x - lx0  << ' ';
		fout_conf << (*p_iter)->p.y - ly0  << ' ';
		fout_conf << (*p_iter)->p.z - lz0  << ' ';
		fout_conf << (*p_iter)->orientation.q0 << ' ';
		fout_conf << (*p_iter)->orientation.q.x << ' ';
		fout_conf << (*p_iter)->orientation.q.y << ' ';
        fout_conf << (*p_iter)->orientation.q.z << endl;
	}
    fout_conf << "B " << number_of_live_bonds << endl;   
	ForAllBond{
        if ( (*bond_iter)->status ){
            int i_p0, i_p1;
            int bondtype;
            if ((*bond_iter)->initial_bond){
                bondtype = 0;
            } else {
                bondtype = 1;
            }
            
            (*bond_iter)->whichparticle(i_p0, i_p1);
            fout_conf << i_p0 << ' ' << i_p1 << ' ' ;
            fout_conf << (*bond_iter)->val_F_norm() << ' ';
            fout_conf << (*bond_iter)->val_F_slid() << ' ';
            fout_conf << (*bond_iter)->val_M_bend() << ' ';
            fout_conf << (*bond_iter)->val_M_tors() << endl;
            //fout_bconf << (*bond_iter)->orientation.q[3] << endl;
        }
    }
}

void DEMsystem::revertStresslet(){
    double A[9];
    double B[9];
    q_rot.makeMatrixA_B(A,B);
    /*
     * S= A.S'.B
     */
    for (int i = 0 ; i < np; i++){
        double s_[5];
        for (int kk=0; kk < 5 ; kk++){
            s_[kk] = vos[i*11 + 6 + kk];
        }
        double smat_[9] = {s_[0], s_[1], s_[2], 
            s_[1], s_[4], s_[3], 
            s_[2], s_[3], -s_[0] - s_[4]};
        double smat[9];
        double tmp[9];
        for (int jc = 0; jc < 3 ; jc++){
            for (int jr = 0; jr < 3 ; jr++){
                int j = 3*jc + jr;
                tmp[j] = 0;
                for (int jj = 0; jj < 3 ; jj++){
                    tmp[j] += smat_[3*jc + jj] * B[3*jj + jr];
                }
            }
        }
        for (int jc = 0; jc < 3 ; jc++){
            for (int jr = 0; jr < 3 ; jr++){
                int j = 3*jc + jr;
                smat[j] = 0;
                for (int jj = 0; jj < 3 ; jj++){
                    smat[j] += A[3*jc + jj] * tmp[3*jj + jr];
                }
            }
        }   
        /* 
         * stresslet[5] = {mat[0], mat[1], mat[2], mat[5], mat[3] } 
         */
        
        particle[i]->stresslet[0] = smat[0];
        particle[i]->stresslet[1] = smat[1];
        particle[i]->stresslet[2] = smat[2];
        particle[i]->stresslet[3] = smat[5]; // This is not mistake.
        particle[i]->stresslet[4] = smat[4];
        

        
    }
}


void DEMsystem::outputData(){
    revertStresslet();
    quaternion q_tmp;
    q_tmp = q_rot*q_rot_total_step;
    fout_data << "# t " << time_simulation << endl;
    fout_data << "# shear_vs_bond " << shear_vs_bond << endl;
    fout_data << "# q_rot_init ";
    fout_data << q_rot_total.q0 << ' ';
    fout_data << q_rot_total.q.x << ' ';
    fout_data << q_rot_total.q.y << ' ';
    fout_data << q_rot_total.q.z << endl;
    fout_data << "# q_rot_cntns ";
    fout_data << q_tmp.q0 << ' ';
    fout_data << q_tmp.q.x << ' ';
    fout_data << q_tmp.q.y << ' ';
    fout_data << q_tmp.q.z << endl; 
    for (int i=0; i < np; i++){
		fout_data << particle[i]->p.x - lx0  << ' ';//1
		fout_data << particle[i]->p.y - ly0  << ' ';//2
		fout_data << particle[i]->p.z - lz0  << ' ';//3
		fout_data << particle[i]->orientation.q0 << ' '; //4
		fout_data << particle[i]->orientation.q.x << ' ';//5
		fout_data << particle[i]->orientation.q.y << ' ';//6
        fout_data << particle[i]->orientation.q.z << ' ';//7
        fout_data << particle[i]->velocity.x << ' ';//8
        fout_data << particle[i]->velocity.y << ' ';//9
        fout_data << particle[i]->velocity.z << ' ';//10
        fout_data << particle[i]->omega.x << ' ';//11
        fout_data << particle[i]->omega.y << ' ';//12
        fout_data << particle[i]->omega.z << ' ';//13
        //////////////////////////////////////////////////////
        // This force is the force acting on particles.
        // This is F_P / F_P0
        fout_data << particle[i]->force.x << ' ';//14
        fout_data << particle[i]->force.y << ' ';//15
        fout_data << particle[i]->force.z << ' ';//16
        // This is T_P / (a F_P0)
        fout_data << particle[i]->torque.x << ' ';//17
        fout_data << particle[i]->torque.y << ' ';//18
        fout_data << particle[i]->torque.z << ' ';//19
        // This is S_H / (a F_H0)
        fout_data << particle[i]->stresslet[0] << ' ';//20
        fout_data << particle[i]->stresslet[1] << ' ';//21
        fout_data << particle[i]->stresslet[2] << ' ';//22
        fout_data << particle[i]->stresslet[3] << ' ';//23
        fout_data << particle[i]->stresslet[4] << endl;//24
	}
    
}

//
//void DEMsystem::outputParticleData(){
//
//    fout_pconf << "# " ; 
//	fout_pconf << ' ' << rup_normal;
//	fout_pconf << ' ' << rup_shear;
//	fout_pconf << ' ' << rup_bend;
//	fout_pconf << ' ' << rup_torsion;
//	fout_pconf << endl;
//    ForAllParticle{
//        fout_pconf << (*p_iter)->p.x - lx0 << ' ';
//        fout_pconf << (*p_iter)->p.y - ly0 << ' ';
//        fout_pconf << (*p_iter)->p.z - lz0 << ' ';
//        fout_pconf << (*p_iter)->orientation.q[0] << ' ';
//        fout_pconf << (*p_iter)->orientation.q[1] << ' ';
//        fout_pconf << (*p_iter)->orientation.q[2] << ' ';
//        fout_pconf << (*p_iter)->orientation.q[3] << endl;
//    }
//    fout_pconf << endl;
//    
//}

void DEMsystem::outputYaplotDEM(){
	static bool firsttime = true;
	if ( firsttime ){firsttime = false;}else{fout_dem << endl;}
	fout_dem << "y 9" << endl;
	fout_dem << "@ 2" << endl;
	ForAllParticle{
		//(*p_iter)->output(fout_dem);
        fout_dem << "c " << (*p_iter)->p.x - lx0;
        fout_dem << ' '  << (*p_iter)->p.y - ly0;
        fout_dem << ' '  << (*p_iter)->p.z - lz0 << endl;
	}
	fout_dem << "y 5" << endl;
	ForAllBond{
        if ( (*bond_iter)->status )
            (*bond_iter)->output_normal(fout_dem);
	}
//	fout_dem << "y 6" << endl;
//	fout_dem << "@ 4" << endl;
//	ForAllBond{
//        if ( (*bond_iter)->status )
//            (*bond_iter)->output_sliding(fout_dem);
//	}
	
	fout_dem << "y 7" << endl;
	fout_dem << "@ 5" << endl;
	ForAllBond{
        if ( (*bond_iter)->status )
            (*bond_iter)->output_bending(fout_dem);
	}
    //	fout_dem << "y 8" << endl;
    //	fout_dem << "@ 6" << endl;
    //	ForAllBond{
    //        if ( (*bond_iter)->status )
    //          (*bond_iter)->output_torsion(fout_dem);
    //	}
}
