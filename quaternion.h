//
//  quaternion.h
//  stodyn
//
//  Created by SETO Ryohei on 24/08/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef quaternion_h
#define quaternion_h
#include "vec3d.h"

class quaternion {    
public:
    double q0;
    vec3d q; 
    double qq_tmp;
    
    quaternion (void){
        reset();
    }
    
	quaternion (const double &_q0, 
                const double &_q1, 
                const double &_q2, 
                const double &_q3){
        q0 = _q0;
        q.set(_q1, _q2 , _q3);
    }
    
    quaternion (const double &_q0, 
                const vec3d &_u){
        q0 = _q0;
        q  = _u;
    }
        
	~quaternion (void){}
    
    void set(const double &_q0, 
             const double &_q1, 
             const double &_q2, 
             const double &_q3){
        q0 = _q0;
        q.set(_q1, _q2, _q3);
    }

    void set(const double &_q0, const vec3d &_v){
        q0 = _q0;
        q  = _v;
    }
    
    void reset(){
        q0 = 1; 
        q.reset();
    }
    
    
    inline quaternion& operator = (const quaternion& qq){
        q0 = qq.q0;
        q  = qq.q;
		return *this;
	}
    
    inline friend quaternion operator * (const quaternion &qa, const quaternion &qb){
        vec3d prod = qa.q0*qb.q + qb.q0*qa.q + cross(qa.q, qb.q);
        return quaternion(qa.q0*qb.q0 - dot(qa.q, qb.q), prod);
	}

    inline vec3d preparer_rotation_quickcalc(){
        qq_tmp = q0*q0 - q.sq_norm();
    }
    inline void rotation_quickcalc(vec3d & v){
        v = qq_tmp*v - 2*(q0*cross(v,q) - dot(v, q)*q);
    }  
    inline void inv_rotation_quickcalc(vec3d & v){
        v = qq_tmp*v + 2*(q0*cross(v,q) + dot(v, q)*q);
    }   
    inline vec3d rotation_quickcalc_r(const vec3d & v){
        return qq_tmp*v - 2*(q0*cross(v,q) - dot(v, q)*q);
    }  
    inline vec3d inv_rotation_quickcalc_r(const vec3d & v){
        return qq_tmp*v + 2*(q0*cross(v,q) + dot(v, q)*q);
    } 
        
    inline vec3d rotation(const vec3d & v){
        return (q0*q0 - q.sq_norm())*v - 2*(q0*cross(v,q) - dot(v, q)*q);
    }  
    inline vec3d ori_forward(const vec3d & v){
        return (q0*q0 - q.sq_norm())*v - 2*(q0*cross(v,q) - dot(v, q)*q);
    }   
    
    inline vec3d inv_rotation(const vec3d & v){
        return (q0*q0 - q.sq_norm())*v + 2*(q0*cross(v,q) + dot(v, q)*q);
    }   
    inline vec3d ori_backward(const vec3d & v){
        return (q0*q0 - q.sq_norm())*v + 2*(q0*cross(v,q) + dot(v, q)*q);
    }
        
    inline quaternion cong (){  
        return quaternion(q0, - q);
    }
    
    inline void infinitesimalRotation(vec3d & dw){
        q0 -= 0.5*dot(q, dw);
        q += 0.5*(q0*dw + cross(q, dw));
    }
    
    inline double norm(){
        return sqrt(q0*q0 + q.sq_norm());
    }
    
    inline void normalize(){
        double norm = sqrt(q0*q0 + q.sq_norm());
        q0 = q0/norm;
        q.x = q.x/norm;
        q.y = q.y/norm;
        q.z = q.z/norm;
    }
    
    void makeMatrixA_B(double *matA, double *matB){
        double q0q0 = q0*q0;
        double q1q1 = q.x*q.x;
        double q2q2 = q.y*q.y;
        double q3q3 = q.z*q.z;
        double q1q2 = q.x*q.y; 
        double q0q3 = q0 *q.z;
        double q1q3 = q.x*q.z;
        double q0q2 = q0 *q.y;
        double q2q3 = q.y*q.z;
        double q0q1 = q0 *q.x;
        matA[0] = q0q0 + q1q1 - q2q2 - q3q3;
        matA[1] = 2*(q1q2 - q0q3);
        matA[2] = 2*(q1q3 + q0q2);
        matA[3] = 2*(q1q2 + q0q3);
        matA[4] = q0q0  - q1q1 + q2q2 - q3q3;
        matA[5] = 2*(q2q3 - q0q1);
        matA[6] = 2*(q1q3 - q0q2);
        matA[7] = 2*(q2q3 + q0q1);
        matA[8] = q0q0 - q1q1 -q2q2 + q3q3;

        matB[0] = q0q0 + q1q1 - q2q2 - q3q3;
        matB[1] = 2*(q1q2 + q0q3);
        matB[2] = 2*(q1q3 - q0q2);
        matB[3] = 2*(q1q2 - q0q3);
        matB[4] = q0q0  - q1q1 + q2q2 - q3q3;
        matB[5] = 2*(q2q3 + q0q1);
        matB[6] = 2*(q1q3 + q0q2);
        matB[7] = 2*(q2q3 - q0q1);
        matB[8] = q0q0 - q1q1 -q2q2 + q3q3;
    }
    
    
    
    
};

#endif

