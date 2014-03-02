/*
 *  system.h
 *  stodyn
 *
 *  Created by seto on 10/04/19.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef system_h
#define system_h 1
#include "stokes.h"
#include "ewald-3fts.h"
#include "ewald-3fts-matrix.h"
#include "ewald-3ft.h"
#include "libiter.h"
#include "vec3d.h"
#include "particle.h"
#include "bond.h"
#include "common.h"
#include "smallmatrix.h"

using namespace std;
class System {
private:    

    ofstream fout_yaplot;
	ofstream fout_SD;
	ofstream fout_sample;	
    vector<vec3d> init_aggregate;
    vec3d U_imposed;
    int lubrication;

	double cl_moi[3];
	double cl_moi_tensor[9];
	double cl_moi_ave;

	int output_precision;
    int output_width;

    // Time-evolutoin of rigid motion
    double cl_mass;
	vec3d cl_mc;
	vec3d x_axis;
	vec3d y_axis;
	vec3d z_axis;
	vec3d p_axes_inertia[3];
    //////////////////////////////////////
    vec3d F_ex;

	/* analise */
	double positive_torque;
	double negative_torque;
	
	/* 0: FD Euler
	 * 1: SD Euler
	 * 2: FD Huen method
	 * 3: SD Huen method
	 */
    
    smallmatrix ** Rfu;
    smallmatrix ** Rfo;
    smallmatrix ** Rfe;
    smallmatrix ** Rtu;
    smallmatrix ** Rto;
    smallmatrix ** Rte;
    smallmatrix ** Rsu;
    smallmatrix ** Rso;
    smallmatrix ** Rse;
    smallmatrix Rag_fu;
    smallmatrix Rag_fo; 
    smallmatrix Rag_fe;
    smallmatrix Rag_tu; 
    smallmatrix Rag_to;
    smallmatrix Rag_te; 
    smallmatrix Rag_su; 
    smallmatrix Rag_so;
    smallmatrix Rag_se;
    
    double *grandResistanceMatrix;
    double Rag[121];
    double Mag[121];

    bool firsttime_yap;
    void shiftCenterOfMass(vector<vec3d> &p);
    void out_xyz(ofstream &out, double xx, double yy, double zz);
	void writeHeader(ofstream &out);
    void memoryForRMatrix();
    void deleteRmatrix();
    void makeSubMatrix();
    void calc_res_3fts_matrix_FDA();

public:
	System();
	~System();

    vec3d clusterOmega(){return cl_omega;}
    int valLubrication(){return lubrication;}
    ////////////////////////////////////////////////////////    
    void set_output_precision(int);
    void setPosition(int, const double&,
                     const double&,const double&);
	void setPosition(int, const vec3d &);
    
    void setMethodHydrodynamicInteraction(int hyd_method);
    void initLibStokes();
    void setSD_IterationMethod();
	void setBox(double lx_, double ly_, double lz_);
    void setFilename(string &);
    void setSimpleShearFlow(double shear_rate_);
	void setImposedFlow_uniform(vec3d &U_);
    void setRigidVelocities();
    void set_dr_from_sdpos();
    void setLubrication(int lub_);
    ////////////////////////////////////////////////////////    
    void importCluster(char*, int skipline);    
    void copyInitialConfiguration( vector<vec3d> &init_positions_);
    void copyPosition(vec3d &p, int i);
    void importTBMfile(char* importfilename,
                       vector<vec3d> &pos_,
                       vector<vec3d> &force_,
                       vector<vec3d> &torque_);
    ////////////////////////////////////////////////////////    
    void openOutputFileStream(); // to be modified    
    void outputYaplot();
	void outputForComparison();

    ////////////////////////////////////////////////////////    
    void calcGyrationRadius();
	void calcMomentOfInertia();

    void calcFTS_RigidAggregate(const vec3d &Uag, const vec3d &Oag);
    void outFTS_RigidAggregate();
    void calcFTS_GrandResistanceMatrix();
	void calcTotalForceTorqueStress();    
    //////////////////////////////////////////////
    /* Torque balanced motion */
    void solveForceFreeRigidMotion(vec3d &Utf, vec3d &Otf);
    void StokesianDynamics();
	void FreeDrainingApproximation();        
    void memAllocateMatrix();
    void calcResistanceMatrixRigidAggregate(); 
    void calcGrandMovMatrix(double* mov);
    
	void output_TorqueBalancedMotion(ofstream &);    
    void output_GrandResistanceMatrix(fstream &);
    void output_ResistanceMatrixRigidAggregate(ofstream &);

    /*************************
	 *   Stokesian Dynamics
	 *************************/
    struct stokes *sd;

    /* method_hydroint
     * 0: Free-draining approximation
     * 1: Stokesian dynamics without lubrication
     * 2: Stokesian dynamics with lubrication
     */
    int method_hydroint;
    int np;
	int nm;
    vec3d *dr;
	double * velocity;
	double * omega;
	double * strain_velocity;
	double * force;
	double * torque;
	double * stresslet;

    vec3d cl_force; // force for the rigid cluster
    vec3d cl_torque; // torque for the rigid cluster
    double cl_stresslet[5]; // stresslet for the rigid cluster
    double cl_stresslet_norm;
	vec3d cl_velocity;
    vec3d cl_omega;
    
	double lx,ly,lz;
	double lx0,ly0,lz0;
	double mass_DLE; // dimension less equation of motion.
	double U0;
	double L0; 
	double F0;
	double T0;	
	double gyration_radius;
	double radius_max;
	double shear_rate0;
	double shear_rate;
	char type_of_flow;
	char type_simu;
	string id_string;
	char filename[64];
};
#endif
