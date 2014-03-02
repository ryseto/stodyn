/*
 *  system.cpp
 *  stodyn
 *
 *  Created by seto on 10/04/19.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "system.h"
#include <map>
#include <algorithm>
#include <ctime>
#include <fstream>

#include "memory-check.h"
#include "stokes.h"
#include "bench.h"
#include "fts.h"
#include "ft.h"
#include "f.h"
/*******************DEM Ryuon **********/

#include "dgetri_c.h" // lapack_inv_() 
#include "ewald.h" // make_matrix_mob_3all ()
#include "matrix.h"
#include "lub-matrix.h"
#include "ewald-3fts-matrix.h"
#include "ewald-3ft-matrix.h"

#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
					   double *w, double *work, int *lwork, int *info );
#endif
extern vector<Particle *> particle;
extern vector<Bond *> bond;

System::System(){
	/* lubrication ----> set by profile */
    cerr << "System()" << endl;
    np = 0;
	lubrication = -1;
    lx = ly = lz = 0;
	dr = NULL;
	velocity = NULL;
	omega = NULL;
	strain_velocity = NULL;
	stresslet = NULL;
	force = NULL;
	torque = NULL;
    Rfu = NULL;
    grandResistanceMatrix = NULL;
}

System::~System(){
    cerr << "~System()" << endl;
    /*
    if (dr)
        DELETE(dr);
    if (velocity)
        DELETE(velocity);
    if (omega)
        DELETE(omega);
    if (strain_velocity)
        DELETE(strain_velocity);
    if (stresslet)
        DELETE(stresslet);
    if (force)
        DELETE(force);
    if (torque)
        DELETE(torque);
    if (grandResistanceMatrix)
        DELETE(grandResistanceMatrix);
    if (Rfu)
        deleteRmatrix(); 
     */
}

//////////////////////
void System::out_xyz(ofstream &out, double xx, double yy, double zz){
	out << xx - lx0 << ' ';
	out << yy - ly0 << ' ';
	out << zz - lz0 << ' ';
}
void System::shiftCenterOfMass(vector<vec3d> &p){
	vec3d cm(0,0,0);
	foreach( vector<vec3d>, p, p_iter){
		cm += (*p_iter);
	}
	cm *= 1.0/p.size();
	foreach( vector<vec3d>, p, p_iter){
		*p_iter -= cm;
	}
}
//////////////////////

void System::set_output_precision(int value){
    output_precision = value;
    output_width = value + 9;
}

/* Set position of a particle 
 * INPUT
 *  i : number of the particle
 *  x,y,z or position(vec3d) : position 
 * OUTPUT
 */
void System::setPosition(int i ,const double &x, const double &y , const double &z){
	int i3 = 3*i;
	sd->pos[i3] = x + lx0;
	sd->pos[i3+1] = y + ly0;
	sd->pos[i3+2] = z + lz0;
}

void System::setPosition(int i, const vec3d &position){
	int j = 3*i;
	sd->pos[j] = position.x;
	sd->pos[j+1] = position.y;
	sd->pos[j+2] = position.z;
}

void System::setMethodHydrodynamicInteraction(int method_hydroint_){
    method_hydroint =  method_hydroint_;
    switch(method_hydroint){
        case 0: //FDA
            cerr << "set to free-draining approximation" << endl;
            lubrication = -1;
			break;
		case 1: //SD without lubrication
            cerr << "set to Stokesian dynamics without lublication" << endl;
            lubrication = 0;
			break;
		case 2: //SD with lubrication
            cerr << "set to Stokesian dynamics wit lublication" << endl;
            lubrication = 1;
			break;
	}
}


/* Initialize libstokes
 * Basic parameters for the libstorks
 * and pos[] are set from init_aggregate.
 */
void System::initLibStokes(){
    if (init_aggregate.empty() && np ==  0){
        cerr << "configuration of particles are not imported yet.";
        exit(1);
    }
    if (np == 0)
        np = init_aggregate.size(); 	/* np : number of all particles  */
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
    
    /*
	 * set pos[] for SD
	 */ 
    
    if (!init_aggregate.empty()){
        double scale = 1.0;
        for (int i = 0; i < np; i++){
            setPosition(i, 
                        scale*init_aggregate[i].x,
                        scale*init_aggregate[i].y,
                        scale*init_aggregate[i].z);
        }
        set_dr_from_sdpos();
    }
 }

void System::setSD_IterationMethod(){
    /* setup ewald summation */
	double ewald_tr = 60.25;
	double xi = xi_by_tratio(sd, ewald_tr);
	double ewald_eps = 1.0e-12;
	stokes_set_xi(sd, xi, ewald_eps);
	//fprintf (stdout, "xi = %f\n", xi);
	//sys->lubmin = 2.0000000001;
	
	sd->lubmin2 = 4.0000000001;
	sd->lubmax = 4.00;
	
	//sd->lubmax = 2.0;
	/* set iter param
	 * INPUT
	 *   solver : string indicating the solver
	 *            sta, sta2, gpb, otmk , or gmres (default)
	 *   eps and log10_eps
	 *   max (and restart)
	 *   debug = 0 : no debug info
	 *         = 1 : iteration numbs and residue
	 *   out   : FILE * to output debug info.
	 */
	stokes_set_iter(sd, "gmres", 2000, 20, 1.0e-6, 0, stdout);
	//stokes_set_iter(sd, "sta", 2000, 20, 1.0e-6, 0, stdout);
}

void System::setBox(double lx_, double ly_, double lz_){
	lx = lx_;
	ly = ly_;
	lz = lz_;
	lx0 = lx/2.0;
	ly0 = ly/2.0;
	lz0 = lz/2.0;
}

void System::setFilename(string &filename_){
	sprintf( filename, "%s", filename_.c_str());
}

void System::copyInitialConfiguration( vector<vec3d> &init_positions_){
    for (int i=0; i< np; i++){
        init_positions_.push_back(init_aggregate[i]);
    }
}

void System::copyPosition(vec3d &p, int i)
{
    int i3 = i*3;
    p.x = sd->pos[i3];
    p.y = sd->pos[i3+1];
    p.z = sd->pos[i3+2];    
}



/* Set simple shear flow 
 *  U^inf = (G*(z-lz/2), 0, 0)
 * INPUT
 * shear_rate_: shear rate
 */
void System::setSimpleShearFlow(double shear_rate_){
	shear_rate0 = shear_rate_;
    sd->Ui[0] = -shear_rate0*lz/2 ;
    sd->Oi[1] = shear_rate0/2;
    sd->Ei[2] = shear_rate0/2;
}

void System::setImposedFlow_uniform(vec3d &U_){
	U_imposed = U_;
	stokes_set_Ui(sd, U_imposed.x, U_imposed.y, U_imposed.z);
}

void System::setRigidVelocities(){
    for (int i = 0; i < np; i++){
		vec3d v = cross(cl_omega, dr[i]) + cl_velocity;
		velocity[i*3] = v.x;
		velocity[i*3+1] = v.y;
		velocity[i*3+2] = v.z;
		omega[i*3] = cl_omega.x;
		omega[i*3+1] = cl_omega.y;
		omega[i*3+2] = cl_omega.z;
	}
}

void System::set_dr_from_sdpos(){    
    for (int i = 0; i < np; i++){
		dr[i].set(sd->pos[i*3 + 0] - lx0,
				  sd->pos[i*3 + 1] - ly0,
				  sd->pos[i*3 + 2] - lz0);
	}  
}

void System::setLubrication(int lub_){
    lubrication = lub_;
}

/* Import the configuration (x,y,z) of particles.
 * The center of mass is set to (0,0,0)
 * INPUT 
 * importfilename: File name
 * skipline: Line number of the header
 */
void System::importCluster(char* importfilename, int skipline){
	sprintf( filename, "%s", importfilename);
	string s_filename = filename;
	int i_backslash = s_filename.find_last_of( "/") + 1;
	int i_extention = s_filename.find( ".dat" );	
	sprintf(filename, "%s",
			(s_filename.substr(i_backslash,i_extention-i_backslash)).c_str());
	

	ifstream fin;
	fin.open( importfilename );
	double x, y, z;
	char buf[1000];
	for (int i = 0; i< skipline; i++){
		fin.getline( buf, 1000 );
	}
	do{
		fin >> x >> y >> z;
        if( fin.fail() )
            break; 
        vec3d new_p(x, y, z);
		init_aggregate.push_back(new_p);
	} while (!fin.eof());
	shiftCenterOfMass( init_aggregate );
}

void System::importTBMfile(char* importfilename,
						   vector<vec3d> &pos_,
						   vector<vec3d> &force_,
						   vector<vec3d> &torque_){
    
	sprintf( filename, "%s", importfilename);
	string s_filename = filename;
	int i_backslash = s_filename.find_last_of( "/") + 1;
	int i_extention = s_filename.find( ".dat" );	
	sprintf(filename, "%s",
			(s_filename.substr( i_backslash, i_extention-i_backslash)).c_str());
	ifstream fin;
	fin.open( importfilename );
	char buf[1000];
    // This number should be checked to use

	int skipline = 3;
	for (int i = 0; i< skipline; i++){
		fin.getline( buf, 1000 );
	}
    fin >> buf >> np;
    fin >> buf >> gyration_radius;
    double Uag[3];
    double Oag[3];
    double Sag[5];
    fin >> buf >> Uag[0] >> Uag[1] >> Uag[2];
    fin >> buf >> Oag[0] >> Oag[1] >> Oag[2];
    fin >> buf >> Sag[0] >> Sag[1] >> Sag[2] >> Sag[3] >> Sag[4] >> buf;
	cl_velocity.set(Uag[0], Uag[1], Uag[2]);
    cl_omega.set(Oag[0], Oag[1], Oag[2]);
    for (int k=0; k<5; k++)
        cl_stresslet[k] = Sag[k];
    
	double pos[3];
	double f[3];
	double t[3];
    double s[5];
    for (int i = 0; i < np ;i++){
		fin >> pos[0] >> pos[1] >> pos[2];
        fin >> f[0] >> f[1] >> f[2] ;
        fin >> t[0] >> t[1] >> t[2] ;
        fin >> s[0] >> s[1] >> s[2] >> s[3] >> s[4];
		pos_.push_back( vec3d(pos[0],pos[1],pos[2]) );
		force_.push_back( vec3d(f[0],f[1],f[2]) );
		torque_.push_back( vec3d(t[0], t[1], t[2] ));
	}

}


void System::openOutputFileStream(){    
	char filenameYaplot[64];
	char filenameSD[64];
	char filenameSample[64];
    
    if ( type_simu == 'R'){
		sprintf(filenameYaplot, "SD%d_%s_%d_%s.yap", lubrication, filename, (int)shear_rate, id_string.c_str() );
        fout_yaplot.open(filenameYaplot);
    }else{
        // for shear and uniform
		if (type_simu == 'S') {
			sprintf(filenameYaplot, "sd%d_shear_%s.yap",lubrication,  filename );
			sprintf(filenameSD, "sd%d_shear_%s.dat",lubrication, filename);
		} else if (type_simu == 'U'){
			sprintf(filenameYaplot, "sd%d_uni_%s.yap", lubrication, filename);
			sprintf(filenameSD, "sd%d_uni_%s.dat", lubrication, filename);
		}
        fout_SD.open(filenameSD);
        fout_SD.setf( ios::scientific, ios::floatfield ); // set math types
        fout_SD.precision( output_precision ); // set precision
    }
}
void System::outputYaplot(){
	if (firsttime_yap == true ){
		firsttime_yap = false;
	} else{
		fout_yaplot << endl;
	}
	fout_yaplot << "y 5\n";
	fout_yaplot << "@ 2\n";
	for (int i = 0; i < np; i ++){
		fout_yaplot << "c " ;
		out_xyz(fout_yaplot, sd->pos[i*3+0] , sd->pos[i*3+1]  , sd->pos[i*3+2] );
		fout_yaplot << endl;
    }

	fout_yaplot << "y 6\n";
	fout_yaplot << "@ 3\n";
	for (int i = 0; i < np; i ++){
		fout_yaplot << "l ";
		out_xyz(fout_yaplot,
				sd->pos[i*3+0],
				sd->pos[i*3+1],
				sd->pos[i*3+2]);
		out_xyz(fout_yaplot, 
				sd->pos[i*3+0] - force[i*3+0],
				sd->pos[i*3+1] - force[i*3+1],
				sd->pos[i*3+2] - force[i*3+2]);
		fout_yaplot << endl;
    }
		
	
	
	fout_yaplot << "y 7\n";
	fout_yaplot << "@ 3\n";
	for (int i = 0; i < np; i ++){
		fout_yaplot << "l ";
		out_xyz(fout_yaplot,
				sd->pos[i*3],
				sd->pos[i*3+1],
				sd->pos[i*3+2]);
		out_xyz(fout_yaplot, 
				sd->pos[i*3] - torque[i*3],
				sd->pos[i*3+1] - torque[i*3+1],
				sd->pos[i*3+2] - torque[i*3+2]);
		fout_yaplot << endl;
	}
	fout_yaplot << "y 10\n";
	fout_yaplot << "@ 0\n";
	fout_yaplot << "l ";
	out_xyz(fout_yaplot, 0,0,0);out_xyz(fout_yaplot,lx,0,0);fout_yaplot << endl;
	fout_yaplot << "l ";
	out_xyz(fout_yaplot, lx,0,0);out_xyz(fout_yaplot,lx,ly,0);fout_yaplot << endl;
	fout_yaplot << "l " ;
	out_xyz(fout_yaplot, lx,ly,0);out_xyz(fout_yaplot,0,ly,0);fout_yaplot << endl;
	fout_yaplot << "l " ;
	out_xyz(fout_yaplot, 0,ly,0);out_xyz(fout_yaplot,0,0,0);fout_yaplot << endl;
	fout_yaplot << "l ";
	out_xyz(fout_yaplot, 0,0,lz);out_xyz(fout_yaplot,lx,0,lz);fout_yaplot << endl;
	fout_yaplot << "l ";
	out_xyz(fout_yaplot, lx,0,lz);out_xyz(fout_yaplot,lx,ly,lz);fout_yaplot << endl;
	fout_yaplot << "l " ;
	out_xyz(fout_yaplot, lx,ly,lz);out_xyz(fout_yaplot,0,ly,lz);fout_yaplot << endl;
	fout_yaplot << "l " ;
	out_xyz(fout_yaplot, 0,ly,lz);out_xyz(fout_yaplot,0,0,lz);fout_yaplot << endl;
}


void System::calcGyrationRadius(){
    /* calculate gyration radius */
	double sum_sq_dist = 0.0;	
	for (int i = 0; i< init_aggregate.size(); i++){
		sum_sq_dist += init_aggregate[i].sq_norm();
	}	
	gyration_radius = sqrt(sum_sq_dist/np);
}

void System::calcMomentOfInertia(){
	for(int k = 0; k < 9; k++)
		cl_moi_tensor[k] = 0.0;
	/* matrix
	 *  0 3 6
	 *  1 4 7
	 *  2 5 8
	 */
	for (int i = 0; i < np; i ++){
		cl_moi_tensor[0] += mass_DLE*(sq(dr[i].y) + sq(dr[i].z)); // xx
		cl_moi_tensor[4] += mass_DLE*(sq(dr[i].x) + sq(dr[i].z)); // yy
		cl_moi_tensor[8] += mass_DLE*(sq(dr[i].x) + sq(dr[i].y)); // zz
		cl_moi_tensor[1] += - mass_DLE*(dr[i].x)*(dr[i].y); // xy
		cl_moi_tensor[2] += - mass_DLE*(dr[i].x)*(dr[i].z); // xz
		cl_moi_tensor[5] += - mass_DLE*(dr[i].y)*(dr[i].z); // yz
	}
	cl_moi_tensor[3]=cl_moi_tensor[1];
	cl_moi_tensor[6]=cl_moi_tensor[2];
	cl_moi_tensor[7]=cl_moi_tensor[5];
    
	/*--------------------------------------------------*/
	/* Find Principal axes of inertia */
	char jobz, uplo;
	jobz = 'V', uplo = 'L';
	int n = 3, lda = 3;
	double *w=NULL;
	double *work=NULL;
	w = new double [n];
	work = new double [1];
	int	lwork = -1, info;
#ifdef OSX
	dsyev_(&jobz, &uplo, 
           (__CLPK_integer *) &n,
           (__CLPK_doublereal *) cl_moi_tensor, 
           (__CLPK_integer *) &lda,
           (__CLPK_doublereal *) w, 
           (__CLPK_doublereal *) work,
           (__CLPK_integer *) &lwork, 
           ( __CLPK_integer *) &info);
#else
    dsyev_(&jobz, &uplo, &n,
           cl_moi_tensor, &lda, w, 
           work, &lwork, &info);
    
#endif
	lwork= (int)work[0];
	DELETE(work);
	work = new double [9];
	jobz = 'V';
	uplo = 'L';
#ifdef OSX
    dsyev_(&jobz, &uplo, 
           (__CLPK_integer *) &n,
           (__CLPK_doublereal *) cl_moi_tensor, 
           (__CLPK_integer *) &lda,
           (__CLPK_doublereal *) w, 
           (__CLPK_doublereal *) work,
           (__CLPK_integer *) &lwork, 
           ( __CLPK_integer *) &info);
#else
    dsyev_(&jobz, &uplo, &n,
           cl_moi_tensor, &lda, w, 
           work, &lwork, &info);
#endif
    if (info != 0 ){
		cerr << "fail to calculate the moment of inetia" << endl;
		exit(1);
	}
	cl_moi[0] = w[0] ;
	cl_moi[1] = w[1] ;
	cl_moi[2] = w[2] ;
	p_axes_inertia[0].set(cl_moi_tensor[0], cl_moi_tensor[1], cl_moi_tensor[2]);
	p_axes_inertia[1].set(cl_moi_tensor[3], cl_moi_tensor[4], cl_moi_tensor[5]);
	p_axes_inertia[2].set(cl_moi_tensor[6], cl_moi_tensor[7], cl_moi_tensor[8]);
	/*--------------------------------------------------*/
	cl_moi_ave = (cl_moi[0] + cl_moi[1] + cl_moi[2])/3;
    DELETE(work);
    DELETE(w);
}

//double System::calcReynoldsNumber(){
//	return (gyration_radius*U0) / (eta/rho_solvant);	
//}


void System::calcTotalForceTorqueStress(){
    vec3d f, t;
    double s[5];   // stress[ xx, xy, xz, yz, yy ];
    cl_torque.reset();
	cl_force.reset();
    for (int i=0; i< 5; i++){
        cl_stresslet[i] = 0.0;
    }
    
	for (int i = 0; i < np; i ++){
		f.set(- force[i*3],  -force[i*3+1],  -force[i*3+2]);
		t.set(-torque[i*3], -torque[i*3+1], -torque[i*3+2]);
        for (int k=0; k < 5; k++)
            s[k] = - stresslet[i*5+k];
        
		cl_force += f;
        cl_torque += (cross(dr[i], f) + t);
        cl_stresslet[0] += ((1.0/3)*(2*dr[i].x*f.x - dr[i].y*f.y - dr[i].z*f.z) + s[0]); //xx
        cl_stresslet[1] += ((1.0/2)*(  dr[i].y*f.x + dr[i].x*f.y)               + s[1]); //xy
        cl_stresslet[2] += ((1.0/2)*(  dr[i].z*f.x + dr[i].x*f.z)               + s[2]); //xz
        cl_stresslet[3] += ((1.0/2)*(  dr[i].z*f.y + dr[i].y*f.z)               + s[3]); //yz
        cl_stresslet[4] += ((1.0/3)*(2*dr[i].y*f.y - dr[i].x*f.x - dr[i].z*f.z) + s[4]); //yy
	}
    cl_stresslet_norm = sqrt(0.5*( sq(cl_stresslet[0])
                                  +sq(cl_stresslet[4])
                                  +sq(cl_stresslet[0]+cl_stresslet[4]))
                             + sq(cl_stresslet[1])+sq(cl_stresslet[2])+sq(cl_stresslet[3]));
    return;
}



////////////////////////
/*
 *
 * output
 *  Utf: torque-free velocity 
 *  Otf: torque-free angular velocity
 */

void System::solveForceFreeRigidMotion(vec3d &Utf,
                                       vec3d &Otf){
    smallmatrix Mag_uf; Mag_uf.allocate_memory(3,3);
    smallmatrix Mag_ut; Mag_ut.allocate_memory(3,3);
    smallmatrix Mag_us; Mag_us.allocate_memory(3,5);
    //////////////////////////
    smallmatrix Mag_of; Mag_of.allocate_memory(3,3);
    smallmatrix Mag_ot; Mag_ot.allocate_memory(3,3);
    smallmatrix Mag_os; Mag_os.allocate_memory(3,5);
    //////////////////////////
    smallmatrix Mag_ef; Mag_ef.allocate_memory(5,3);
    smallmatrix Mag_et; Mag_et.allocate_memory(5,3);
    smallmatrix Mag_es; Mag_es.allocate_memory(5,5);
    
    for (int i=0; i < 121; i++){
        Mag[i] = Rag[i];
    }
    lapack_inv_ (11, Mag);
    
    for (int l=0; l< 3; l++){
        for (int k=0; k < 3 ; k++){
            Mag_uf.element[l][k] = Mag[k   + 11*l];
            Mag_ut.element[l][k] = Mag[3+k + 11*l];
            Mag_of.element[l][k] = Mag[k   + 11*(3+l)];
            Mag_ot.element[l][k] = Mag[3+k + 11*(3+l)];            
        }
        for (int k=0; k < 5 ; k++){
            Mag_us.element[l][k] = Mag[6+k + 11*l];
            Mag_os.element[l][k] = Mag[6+k + 11*(3+l)];
        } 
    }
    for (int l=0; l< 5; l++){
        for (int k=0; k < 3 ; k++){
            Mag_ef.element[l][k] = Mag[k   + 11*(6+l)];
            Mag_et.element[l][k] = Mag[3+k + 11*(6+l)];
        }
        for (int k=0; k < 5 ; k++){
            Mag_es.element[l][k] = Mag[6+k + 11*(6+l)];
        }
    }
    
    double *array_inv_Mag_es;
    array_inv_Mag_es = new double [25];
    for (int l=0; l< 5; l++){
        for (int k=0; k < 5 ; k++){
            array_inv_Mag_es[k + l*5] = Mag_es.element[l][k];
        }
    }
    
    lapack_inv_ (5, array_inv_Mag_es);
    smallmatrix inv_Mag_es;
    
    inv_Mag_es.allocate_memory(5,5);
    for (int l=0; l< 5; l++){
        for (int k=0; k < 5 ; k++){
            inv_Mag_es.element[l][k] = array_inv_Mag_es[k + l*5];
        }
    }
    
    double *Sag;
    Sag = new double [5];
    for (int l=0; l < 5; l++) {Sag[l] = 0;}
    inv_Mag_es.multiplyVector( sd->Ei, Sag );
    
    double *Uag, *Oag;
    Uag = new double[3];
    Oag = new double[3];
    
    Mag_us.multiplyVector( Sag, Uag );
    for (int l=0; l < 3 ; l++) Uag[l] = - Uag[l];
    Mag_os.multiplyVector( Sag, Oag );
    for (int l=0; l < 3 ; l++) Oag[l] = - Oag[l];
    
    Oag[0] += sd->Oi[0];
    Oag[1] += sd->Oi[1];
    Oag[2] += sd->Oi[2];
    
    Utf.set(Uag[0], Uag[1], Uag[2]);
    Otf.set(Oag[0], Oag[1], Oag[2]);
    
    DELETE(Sag);
    DELETE(Uag);
    DELETE(Oag);
    DELETE(array_inv_Mag_es);
    return;
}

void System::StokesianDynamics(){
	if ( type_of_flow == 's'){
		if (sd->twobody_lub == 0){
            solve_res_3fts(sd, velocity, omega, strain_velocity,
                           force, torque, stresslet);
            //solve_res_3fts_matrix(sd, 
            //                  velocity, omega, strain_velocity,
            //                  force, torque, stresslet);
		} else {
			solve_res_lub_3fts(sd,
							   velocity, omega, strain_velocity,
							   force, torque, stresslet);
		}
		
	} else {
		if (sd->twobody_lub == 0){
            //sd->version = 1;
			//solve_res_3ft(sd, velocity, omega, force, torque);
            //solve_res_3ft_matrix(sd, velocity, omega, force, torque);
            solve_res_3fts(sd, velocity, omega, strain_velocity,
                           force, torque, stresslet);
		} else {
			//solve_res_lub_3ft(sd, velocity, omega, force, torque);
			solve_res_lub_3fts(sd,
							   velocity, omega, strain_velocity,
                                force, torque, stresslet);
		}
	}
    
}

// (U,O,S) = mov*(F,T,E)
void System::calcGrandMovMatrix(double* mov){
    calc_mob_3fts_matrix(sd, mov);
}


void System::FreeDrainingApproximation(){
	vec3d Uoo(sd->Ui[0], sd->Ui[1],	sd->Ui[2]);
	vec3d Ooo(sd->Oi[0], sd->Oi[1],	sd->Oi[2]);
    double Eoo[5];
    for (int i=0; i < 5; i++){
        Eoo[i] = sd->Ei[i];
    }
	for (int i = 0; i < np; i ++){
		vec3d p(sd->pos[i*3+0], sd->pos[i*3+1], sd->pos[i*3+2]);  
		vec3d u = Uoo + cross(Ooo, dr[i]) + dr[i].product_rate_of_strain( sd->Ei );
		force[i*3+0] = -(u.x - velocity[i*3+0]);
		force[i*3+1] = -(u.y - velocity[i*3+1]);
		force[i*3+2] = -(u.z - velocity[i*3+2]);
		torque[i*3+0] = -(4.0/3)*(Ooo.x - omega[i*3+0])/shear_rate0;
		torque[i*3+1] = -(4.0/3)*(Ooo.y - omega[i*3+1])/shear_rate0;
		torque[i*3+2] = -(4.0/3)*(Ooo.z - omega[i*3+2])/shear_rate0;
        stresslet[i*5+0] = -(10.0/9)*Eoo[0]/shear_rate0;
        stresslet[i*5+1] = -(10.0/9)*Eoo[1]/shear_rate0;
        stresslet[i*5+2] = -(10.0/9)*Eoo[2]/shear_rate0;
        stresslet[i*5+3] = -(10.0/9)*Eoo[3]/shear_rate0;
        stresslet[i*5+4] = -(10.0/9)*Eoo[4]/shear_rate0;
	}
}

/*
 * Resistance matrix + rigid aggregate
 */ 

void System::memoryForRMatrix(){
    Rfu = new smallmatrix * [np];
    Rfo = new smallmatrix * [np];
    Rfe = new smallmatrix * [np];
    Rtu = new smallmatrix * [np];
    Rto = new smallmatrix * [np];
    Rte = new smallmatrix * [np];
    Rsu = new smallmatrix * [np];
    Rso = new smallmatrix * [np];
    Rse = new smallmatrix * [np];
    for (int i = 0; i < np; i++){
        Rfu[i] = new smallmatrix [np];
        Rfo[i] = new smallmatrix [np];
        Rfe[i] = new smallmatrix [np];
        Rtu[i] = new smallmatrix [np];
        Rto[i] = new smallmatrix [np];
        Rte[i] = new smallmatrix [np];
        Rsu[i] = new smallmatrix [np];
        Rso[i] = new smallmatrix [np];
        Rse[i] = new smallmatrix [np];
        for (int j = 0; j < np; j++){
            Rfu[i][j].allocate_memory(3,3);
            Rfo[i][j].allocate_memory(3,3);
            Rfe[i][j].allocate_memory(3,5);
            Rtu[i][j].allocate_memory(3,3);
            Rto[i][j].allocate_memory(3,3);
            Rte[i][j].allocate_memory(3,5);
            Rsu[i][j].allocate_memory(5,3);
            Rso[i][j].allocate_memory(5,3);
            Rse[i][j].allocate_memory(5,5);
        }
    }
}

void System::makeSubMatrix(){
    int n11 = np*11;
    for (int i = 0 ; i < np ;i++ ){
        for (int j = 0 ; j < np ; j++ ){
            for(int l=0; l < 3; l++){ 
                for (int k=0; k < 3 ; k ++) {
                    Rfu[i][j].element[l][k] = grandResistanceMatrix[11*j+k     + n11*(11*i+l)]; // ok
                    Rfo[i][j].element[l][k] = grandResistanceMatrix[11*j+3+k   + n11*(11*i+l)]; // ok
                    Rtu[i][j].element[l][k] = grandResistanceMatrix[11*j+k     + n11*(11*i+3+l)];//ok
                    Rto[i][j].element[l][k] = grandResistanceMatrix[11*j+3+k   + n11*(11*i+3+l)];//ok
                }           
                for (int k=0; k < 5 ; k ++) {
                    Rfe[i][j].element[l][k] = grandResistanceMatrix[11*j+6+k   + n11*(11*i+l)]; // ok;
                    Rte[i][j].element[l][k] = grandResistanceMatrix[11*j+6+k   + n11*(11*i+3+l)];//ok;
                }
            }
            for(int l=0; l < 5; l++){ 
                for (int k=0; k < 3 ; k ++) {
                    Rsu[i][j].element[l][k] = grandResistanceMatrix[11*j+k     + n11*(11*i+6+l)]; // ok
                    Rso[i][j].element[l][k] = grandResistanceMatrix[11*j+3+k   + n11*(11*i+6+l)]; // ok
                }            
                for (int k=0; k < 5 ; k ++) {
                    Rse[i][j].element[l][k] = grandResistanceMatrix[11*j+6+k   + n11*(11*i+6+l)]; // ok
                }
            }                        
        }
    }
}

void System::calcFTS_GrandResistanceMatrix(){
    if (grandResistanceMatrix == NULL){
        cerr << "Obtain the grand R-matrix first" << endl;
        exit(1);
    }
    setRigidVelocities();   
    int n11 = np*11;
    double *voe;
    voe = new double[n11];
    double *fts;
    fts = new double [n11];
    for (int j = 0 ; j < np ; j++){
        for (int k=0; k < 3 ; k++){
            // @@@@@@@@@@@@@@@@@ REWRITE
            double uinf[3];
            uinf[0] = shear_rate0*dr[j].z;
            uinf[1] = 0;
            uinf[2] = 0;
            voe[j*11+k] = (velocity[j*3+k] - uinf[k]);
            voe[j*11+3+k] = (omega[j*3+k] - sd->Oi[k]);
        }
        for (int k=0; k < 5 ; k++){
            voe[j*11+6+k] = (0.0 -sd->Ei[k]); 
        }
    }
    
    for(int i = 0; i < n11; i++){
        fts[i] = 0.0;
        for (int j = 0 ; j < n11 ; j++){
            fts[i] += grandResistanceMatrix[i*n11 + j]*voe[j];
        }
    }
    
    for (int i = 0 ; i < np ; i++){
        for (int k = 0; k < 3; k++){
            force[3*i+k] =  fts[11*i +k];
            torque[3*i+k] = fts[11*i +3+k];
        }
        for (int k=0; k < 5; k++){
            stresslet[5*i+k] = fts[11*i +6+k];
        }
    }
    
    DELETE(fts);
    DELETE(voe); 
}



void System::calcResistanceMatrixRigidAggregate(){
    int n11 = np*11;
    try{
        grandResistanceMatrix = new double [n11*n11];
    }catch  (bad_alloc &){
        cerr << "bad_alloc at System::init()" << endl;
        exit(1);
    }
    if (method_hydroint == 0){
        cerr << "FDA" << endl;
        calc_res_3fts_matrix_FDA();
    } else {
        cerr << "SD" << endl;
        calc_res_3fts_matrix(sd, grandResistanceMatrix );
    }
    
    memoryForRMatrix();
    makeSubMatrix();     
    /////////////////////////////////////////////////////
    Rag_fu.allocate_memory(3,3);Rag_fo.allocate_memory(3,3);Rag_fe.allocate_memory(3,5);
    Rag_tu.allocate_memory(3,3);Rag_to.allocate_memory(3,3);Rag_te.allocate_memory(3,5);
    Rag_su.allocate_memory(5,3);Rag_so.allocate_memory(5,3);Rag_se.allocate_memory(5,5);
    
    smallmatrix *til_Rfu, *til_Rfo, *til_Rfe;
    smallmatrix *hat_Rsu, *hat_Rso, *hat_Rse;
    til_Rfu = new smallmatrix [np];
    til_Rfo = new smallmatrix [np];
    til_Rfe = new smallmatrix [np];
    hat_Rsu = new smallmatrix [np];
    hat_Rso = new smallmatrix [np];
    hat_Rse = new smallmatrix [np];
    for (int i=0; i<np; i++){
        til_Rfu[i].allocate_memory(3, 3);
        til_Rfo[i].allocate_memory(3, 3);
        til_Rfe[i].allocate_memory(3, 5);
        hat_Rsu[i].allocate_memory(5, 3);
        hat_Rso[i].allocate_memory(5, 3);
        hat_Rse[i].allocate_memory(5, 5);
    }
    smallmatrix tmpmat3x3,tmpmat3x5,tmpmat5x3,tmpmat5x5;
    tmpmat3x3.allocate_memory(3,3);tmpmat3x5.allocate_memory(3,5);
    tmpmat5x3.allocate_memory(5,3);tmpmat5x5.allocate_memory(5,5); 
    
    /////////////////////////////////////////////////////
    /* Calculation for Rag_fu  */
    Rag_fu.set_zero();
    for (int i = 0; i< np; i++){
        til_Rfu[i].set_zero();
        for (int j = 0; j< np; j++){
            til_Rfu[i] += Rfu[i][j];
        }
        Rag_fu += til_Rfu[i];
    }
    /* Calculation for Rag_fo  */
    Rag_fo.set_zero();
    for (int i = 0; i< np; i++){
        til_Rfo[i].set_zero();
        for (int j = 0; j< np; j++){
            til_Rfo[i] += Rfo[i][j];
            (Rfu[i][j]).matprod1(tmpmat3x3, dr[j].x, dr[j].y, dr[j].z);
            til_Rfo[i] += tmpmat3x3;
        }
        Rag_fo += til_Rfo[i];
    }
    /* Calculation for Rag_fe  */
    Rag_fe.set_zero();
    for (int i = 0; i< np; i++){
        til_Rfe[i].set_zero();
        for (int j = 0; j< np; j++){
            til_Rfe[i] += Rfe[i][j];
            (Rfu[i][j]).matprod2(tmpmat3x5, dr[j].x, dr[j].y, dr[j].z);
            til_Rfe[i] += tmpmat3x5;
        }
        Rag_fe += til_Rfe[i];
    }
    
    /* Calculation for Rag_tu  */
    Rag_tu.set_zero();
    for (int i = 0; i< np; i++){ 
        for (int j = 0; j< np; j++){
            Rag_tu += Rtu[i][j];
        }
        (til_Rfu[i]).mat_vector_prod(tmpmat3x3, dr[i].x, dr[i].y, dr[i].z);        
        Rag_tu += tmpmat3x3;
    }
    
    /* Calculation for Rag_to  */
    Rag_to.set_zero();
    for (int i = 0; i< np; i++){
        for (int j = 0; j< np; j++){
            Rag_to += Rto[i][j];
            (Rtu[i][j]).matprod1(tmpmat3x3, dr[j].x, dr[j].y, dr[j].z);
            Rag_to += tmpmat3x3;    
        }
        (til_Rfo[i]).mat_vector_prod(tmpmat3x3, dr[i].x, dr[i].y, dr[i].z);
        Rag_to += tmpmat3x3;
    }
    
    /* Calculation for Rag_te  */
    Rag_te.set_zero();
    for (int i = 0; i< np; i++){
        for (int j = 0; j< np; j++){
            Rag_te += Rte[i][j];
            (Rtu[i][j]).matprod2(tmpmat3x5, dr[j].x, dr[j].y, dr[j].z);
            Rag_te += tmpmat3x5;
        }
        (til_Rfe[i]).mat_vector_prod(tmpmat3x5, dr[i].x, dr[i].y, dr[i].z);
        Rag_te += tmpmat3x5;
    }
    
    /* Calculation for Rag_su  */
    for (int i = 0; i< np; i++){
        for (int k=0; k < 3 ; k++){
            hat_Rsu[i].element[0][k] 
            = (1.0/3)*(2*dr[i].x*til_Rfu[i].element[0][k]
                       - dr[i].y*til_Rfu[i].element[1][k] - dr[i].z*til_Rfu[i].element[2][k]);
            hat_Rsu[i].element[1][k] 
            = (1.0/2)*(  dr[i].x*til_Rfu[i].element[1][k] + dr[i].y*til_Rfu[i].element[0][k]);
            hat_Rsu[i].element[2][k] 
            = (1.0/2)*(  dr[i].x*til_Rfu[i].element[2][k] + dr[i].z*til_Rfu[i].element[0][k]);
            hat_Rsu[i].element[3][k] 
            = (1.0/2)*(  dr[i].y*til_Rfu[i].element[2][k] + dr[i].z*til_Rfu[i].element[1][k]);
            hat_Rsu[i].element[4][k] 
            = (1.0/3)*(2*dr[i].y*til_Rfu[i].element[1][k]
                       - dr[i].x*til_Rfu[i].element[0][k] - dr[i].z*til_Rfu[i].element[2][k]);
        }
    }
    Rag_su.set_zero();
    for (int i = 0; i < np; i++){
        for (int j = 0; j < np; j++){
            Rag_su +=  Rsu[i][j];
        }
        Rag_su += hat_Rsu[i];
    }
    
    /* Calculation for Rag_so  */
    for (int i = 0; i< np; i++){
        for (int k=0; k < 3 ; k++){
            hat_Rso[i].element[0][k] 
            = (1.0/3)*(2*dr[i].x*til_Rfo[i].element[0][k]
                       - dr[i].y*til_Rfo[i].element[1][k] - dr[i].z*til_Rfo[i].element[2][k]); // xx
            hat_Rso[i].element[1][k] 
            = (1.0/2)*(dr[i].x*til_Rfo[i].element[1][k] + dr[i].y*til_Rfo[i].element[0][k]); // xy
            hat_Rso[i].element[2][k] 
            = (1.0/2)*(dr[i].x*til_Rfo[i].element[2][k] + dr[i].z*til_Rfo[i].element[0][k]); // xz
            hat_Rso[i].element[3][k] 
            = (1.0/2)*(dr[i].y*til_Rfo[i].element[2][k] + dr[i].z*til_Rfo[i].element[1][k]); // yz 
            hat_Rso[i].element[4][k] 
            = (1.0/3)*(2*dr[i].y*til_Rfo[i].element[1][k]
                       - dr[i].x*til_Rfo[i].element[0][k] - dr[i].z*til_Rfo[i].element[2][k]); // yy       
        }
    }
    Rag_so.set_zero();
    for (int i = 0; i < np; i++){
        for (int j = 0; j< np; j++){
            Rag_so += Rso[i][j];
            (Rsu[i][j]).matprod1(tmpmat5x3, dr[j].x, dr[j].y, dr[j].z);
            Rag_so += tmpmat5x3;
        }
        Rag_so += hat_Rso[i];
    }
    
    /* Calculation for Rag_se  */
    for (int i = 0; i< np; i++){
        for (int k=0; k < 5 ; k++){
            hat_Rse[i].element[0][k] 
            = (1.0/3)*(2*dr[i].x*til_Rfe[i].element[0][k]
                       - dr[i].y*til_Rfe[i].element[1][k] - dr[i].z*til_Rfe[i].element[2][k]);
            hat_Rse[i].element[1][k] 
            = (1.0/2)*(dr[i].x*til_Rfe[i].element[1][k] + dr[i].y*til_Rfe[i].element[0][k]);
            hat_Rse[i].element[2][k] 
            = (1.0/2)*(dr[i].x*til_Rfe[i].element[2][k] + dr[i].z*til_Rfe[i].element[0][k]);
            hat_Rse[i].element[3][k] 
            = (1.0/2)*(dr[i].y*til_Rfe[i].element[2][k] + dr[i].z*til_Rfe[i].element[1][k]);
            hat_Rse[i].element[4][k] 
            = (1.0/3)*(2*dr[i].y*til_Rfe[i].element[1][k] 
                       - dr[i].x*til_Rfe[i].element[0][k] - dr[i].z*til_Rfe[i].element[2][k]);            
        }
    }
    
    Rag_se.set_zero();
    for (int i = 0; i < np; i++){
        for (int j = 0; j < np; j++){
            Rag_se += Rse[i][j];
            (Rsu[i][j]).matprod2(tmpmat5x5, dr[j].x, dr[j].y, dr[j].z);
            Rag_se += tmpmat5x5;
        }
        Rag_se += hat_Rse[i];
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    for (int l=0; l < 3; l++){
        for (int k=0; k < 3 ; k++){
            Rag[11*l      + k]   = Rag_fu.element[l][k];
            Rag[11*l      + 3+k] = Rag_fo.element[l][k];
            Rag[11*(l+3)  + k]   = Rag_tu.element[l][k];
            Rag[11*(l+3)  + 3+k] = Rag_to.element[l][k];
            
        }
        for (int k=0; k < 5 ; k++){
            Rag[11*l      + 6+k] = Rag_fe.element[l][k];
            Rag[11*(l+3)  + 6+k] = Rag_te.element[l][k];
        } 
    }
    for (int l=0; l < 5; l++){
        for (int k=0; k < 3 ; k++){
            Rag[11*(l+6)  + k]   = Rag_su.element[l][k];
            Rag[11*(l+6)  + 3+k] = Rag_so.element[l][k];
        }
        for (int k=0; k < 5 ; k++){
            Rag[11*(l+6)  + 6+k] = Rag_se.element[l][k];
        }
    }
    DELETE(til_Rfu);
    DELETE(til_Rfo);
    DELETE(til_Rfe);
    DELETE(hat_Rsu);
    DELETE(hat_Rso);
    DELETE(hat_Rse);
    
    return;
}

void System::outFTS_RigidAggregate(){
    cerr << "Fag = ";
    cl_force.cerr();
    cerr << "Tag = ";
    cl_torque.cerr();
    cerr << "Sag = " ;
    for(int l=0; l < 5; l++)
        cerr << cl_stresslet[l] << ' ' ;
    cerr << endl;
    return;
}

void System::calcFTS_RigidAggregate(const vec3d &Uag, 
                                    const vec3d &Oag){
    cl_velocity = Uag;
    cl_omega = Oag;
    
    double Fag[3];
    double Tag[3];
    double Sag[5];
    
    for(int l=0; l<3; l++){
        Fag[l] = 0.0;
        Tag[l] = 0.0;
    }
    for(int l=0; l<5; l++){
        Sag[l]=0.0;
    }
    for (int l=0; l < 3; l++){
        Fag[l] += -Rag_fu.element[l][0]*(cl_velocity.x);
        Fag[l] += -Rag_fu.element[l][1]*(cl_velocity.y);
        Fag[l] += -Rag_fu.element[l][2]*(cl_velocity.z);
        Fag[l] += -Rag_fo.element[l][0]*(cl_omega.x - sd->Oi[0]);
        Fag[l] += -Rag_fo.element[l][1]*(cl_omega.y - sd->Oi[1]);
        Fag[l] += -Rag_fo.element[l][2]*(cl_omega.z - sd->Oi[2]);
        Fag[l] += -Rag_fe.element[l][0]*(-sd->Ei[0]);
        Fag[l] += -Rag_fe.element[l][1]*(-sd->Ei[1]);
        Fag[l] += -Rag_fe.element[l][2]*(-sd->Ei[2]);
        Fag[l] += -Rag_fe.element[l][3]*(-sd->Ei[3]);
        Fag[l] += -Rag_fe.element[l][4]*(-sd->Ei[4]);        
        Tag[l] += -Rag_tu.element[l][0]*(cl_velocity.x);
        Tag[l] += -Rag_tu.element[l][1]*(cl_velocity.y);
        Tag[l] += -Rag_tu.element[l][2]*(cl_velocity.z);
        Tag[l] += -Rag_to.element[l][0]*(cl_omega.x - sd->Oi[0]);
        Tag[l] += -Rag_to.element[l][1]*(cl_omega.y - sd->Oi[1]);
        Tag[l] += -Rag_to.element[l][2]*(cl_omega.z - sd->Oi[2]);
        Tag[l] += -Rag_te.element[l][0]*(-sd->Ei[0]);
        Tag[l] += -Rag_te.element[l][1]*(-sd->Ei[1]);
        Tag[l] += -Rag_te.element[l][2]*(-sd->Ei[2]);
        Tag[l] += -Rag_te.element[l][3]*(-sd->Ei[3]);
        Tag[l] += -Rag_te.element[l][4]*(-sd->Ei[4]);
    }
    
    for (int l=0; l < 5; l++){
        Sag[l] += -Rag_su.element[l][0]*(cl_velocity.x);
        Sag[l] += -Rag_su.element[l][1]*(cl_velocity.y);
        Sag[l] += -Rag_su.element[l][2]*(cl_velocity.z);
        Sag[l] += -Rag_so.element[l][0]*(cl_omega.x - sd->Oi[0]);
        Sag[l] += -Rag_so.element[l][1]*(cl_omega.y - sd->Oi[1]);
        Sag[l] += -Rag_so.element[l][2]*(cl_omega.z - sd->Oi[2]);
        Sag[l] += -Rag_se.element[l][0]*(-sd->Ei[0]);
        Sag[l] += -Rag_se.element[l][1]*(-sd->Ei[1]);
        Sag[l] += -Rag_se.element[l][2]*(-sd->Ei[2]);
        Sag[l] += -Rag_se.element[l][3]*(-sd->Ei[3]);
        Sag[l] += -Rag_se.element[l][4]*(-sd->Ei[4]);
    }
    
    cl_force.set(Fag[0],Fag[1],Fag[2]);
    cl_torque.set(Tag[0],Tag[1],Tag[2]);
    for(int l=0; l<5; l++)
        cl_stresslet[l] = Sag[l];
    
    cl_stresslet_norm = sqrt(0.5*(sq(cl_stresslet[0])+sq(cl_stresslet[4])+sq(cl_stresslet[0]+cl_stresslet[4]))
                             + sq(cl_stresslet[1])+sq(cl_stresslet[2])+sq(cl_stresslet[3]));
    
}

void System::calc_res_3fts_matrix_FDA(){
    double self_force = 1.0 ;
    double self_torque = 4.0/3;
    double self_stresslet = 10.0/9;
    int n11 = np*11;
    for (int i=0; i < n11; i++){
        for (int j=0; j < n11; j++){
            grandResistanceMatrix[ i*n11 + j]=0.0;
        }
    }
    for (int i=0; i < np ; i++){
        int i11 = i*11;
        for (int k=0; k < 3; k++){
            grandResistanceMatrix[(i11+k)*n11   + i11+k] = self_force;
            grandResistanceMatrix[(i11+3+k)*n11 + i11+3+k] = self_torque;
        }
        for (int k=0; k < 3; k++){
            grandResistanceMatrix[(11*i+6+k)*n11 + 11*i+6+k] = self_stresslet;
        }
    }
}

void System::deleteRmatrix(){
    cerr << "np" << np << endl;
    for (int i = 0; i < np; i++){
        DELETE(Rfu[i]);
        DELETE(Rfo[i]);
        DELETE(Rfe[i]);
        DELETE(Rtu[i]);
        DELETE(Rto[i]);
        DELETE(Rte[i]);
        DELETE(Rsu[i]);
        DELETE(Rso[i]);
        DELETE(Rse[i]);
    }   
    DELETE(Rfu);
    DELETE(Rfo);
    DELETE(Rfe);
    DELETE(Rtu);
    DELETE(Rto);
    DELETE(Rte);
    DELETE(Rsu);
    DELETE(Rso);
    DELETE(Rse);
}


/* 
 * output torque-balanced motion 
 * INPUT
 * output: File stream
 */
void System::output_TorqueBalancedMotion(ofstream &fout){
    fout << "# The torque balanced motion. Fag=0 and Tag=0" << endl; //1
    fout << "# The main data is given from line" << endl; //2
    fout << "# Data : pos(1,2,3) f(4,5,6) t(7,8,9) s(10,11,12,13,14)" << endl; //3
	fout << "n " << np << endl; //4
    fout << "Rg " << gyration_radius << endl;//5
    fout << "Uag ";
    fout << setprecision(output_precision);
    fout << cl_velocity.x << ' ' << cl_velocity.y << ' ' << cl_velocity.z << endl; //6
    fout << "Oag ";
    fout << setprecision(output_precision);
    fout << cl_omega.x << ' ' << cl_omega.y << ' ' << cl_omega.z << endl; // 7
    fout << "Sag ";
    fout << setprecision(output_precision);
    fout << cl_stresslet[0] << ' '; //8
    fout << cl_stresslet[1] << ' '; //8
    fout << cl_stresslet[2] << ' '; //8
    fout << cl_stresslet[3] << ' '; //8
    fout << cl_stresslet[4] << ' '; //8
    fout << cl_stresslet_norm << endl; //8
    
    // stress[ xx, xy, xz, yz, yy ];
	for (int i = 0; i < np; i++){
		fout << setprecision(output_precision);
		fout << dr[i].x << ' '; //1 
		fout << dr[i].y << ' ';//2
		fout << dr[i].z << ' ';//3
		fout << -force[3*i] << ' ';//4
		fout << -force[3*i+1] << ' ';//5
		fout << -force[3*i+2] << ' ';//6
		fout << -torque[3*i] << ' ';//7
		fout << -torque[3*i+1] << ' ';//8
		fout << -torque[3*i+2] << ' ';//9
        fout << -stresslet[5*i] << ' '; //10
        fout << -stresslet[5*i+1] << ' '; //11
        fout << -stresslet[5*i+2] << ' '; //12
        fout << -stresslet[5*i+3] << ' '; //13
        fout << -stresslet[5*i+4] << endl; //14
	}
}

void System::output_GrandResistanceMatrix(fstream &output){
    int n11 = np*11;
    int datasize = n11*n11*sizeof(grandResistanceMatrix[0]);
    //  output << np ;
    //  int writepoint = output.tellp();
    //	output.seekp(writepoint, ios::beg);
    //  cerr << sizeof(grandResistanceMatrix) << endl;
	output.write((char *) grandResistanceMatrix,
                 datasize);
    output.flush();
    
    //grandResistanceMatrix
    
}

void System::output_ResistanceMatrixRigidAggregate(ofstream &output){
    output << "# Resistance matrix of rigid aggregate" <<endl;
    for (int i=0; i< 11 ; i++){    
        for (int j=0; j< 11 ; j++){    
            output << Rag[i*11+j] << ' ';
        }
        output << endl;
    }
}

void System::outputForComparison(){
	fout_SD << "# N " << np << endl; 
    fout_SD << "# Re " << 0 << endl;
//    fout_SD << "# Rg " << gyration_radius << endl;
//    fout_SD << "# totalforce " << cl_force.x << ' ' << cl_force.y << ' ' << cl_force.z  << endl;
//    fout_SD << "# totaltorque " << cl_torque.x << ' ' << cl_torque.y << ' ' << cl_torque.z  << endl;
//    fout_SD << "# totalstress " << cl_stresslet[0] << ' ' << cl_stresslet[1] << ' ' << cl_stresslet[2] << ' ';
//    fout_SD << cl_stresslet[3] << ' ' << cl_stresslet[4] << ' ' << cl_stresslet_norm << endl;
	//	fout_sample << "# a " << 1.00 << endl; 
	//	fout_sample << "# shear-rate " << 1.00 << endl; 
	//	fout_sample << "# "<< endl;
	fout_SD << "# x y z vx vy vz ox oy oz fx fy fz tx ty tz"<< endl;
	cerr << sd->pos[0] << ' ' << lx0 << endl;
    
	for (int i = 0; i < np; i ++){


		fout_SD << setw(output_width) << sd->pos[i*3]-lx0;
		fout_SD << setw(output_width) << sd->pos[i*3+1]-ly0;
		fout_SD << setw(output_width) << sd->pos[i*3+2]-lz0;
		fout_SD << setw(output_width) << velocity[i*3];
		fout_SD << setw(output_width) << velocity[i*3+1];
		fout_SD << setw(output_width) << velocity[i*3+2];
		fout_SD << setw(output_width) << omega[i*3];
		fout_SD << setw(output_width) << omega[i*3+1];
		fout_SD << setw(output_width) << omega[i*3+2];
		fout_SD << setw(output_width) << - force[i*3];
		fout_SD << setw(output_width) << - force[i*3+1];
		fout_SD << setw(output_width) << - force[i*3+2];
		fout_SD << setw(output_width) << - torque[i*3];
		fout_SD << setw(output_width) << - torque[i*3+1];
		fout_SD << setw(output_width) << - torque[i*3+2];
		fout_SD << endl;
	}

}




