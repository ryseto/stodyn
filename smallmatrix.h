//
//  smallmatrix.h
//  stodyn
//
//  Created by seto on 14/07/2011.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef smallmatrix_h
#define smallmatrix_h 1
#include <iostream>
#include <iomanip>
#include <cmath>
#include "my_utilities.h"
#define DELETE(x) delete [] x; x = NULL;

class smallmatrix {
private:
    // matrix[row][column];
    int row; 
    int column; 
public:
	/* variables */
	double **element;

	/* constructor/destructor */
    smallmatrix (void){
        element = NULL;
        return;
    }
    
	~smallmatrix(void){
        for (int i=0; i < row ; i++){
            DELETE(element[i]);
        }
        DELETE(element);
        return;
    }
	inline smallmatrix& operator = (const smallmatrix& mat){
        for (int i=0; i < row; i++){
            for (int j=0; j < column; j++){
                element[i][j] += mat.element[i][j];
            }
        }
		return *this;
	}
    
    void allocate_memory(int row_, int column_){
        row = row_;
        column = column_;
        try{
            element = new double * [row];
            for (int i=0; i < row ; i++){
                element[i] = new double [column];
            }
        } catch (bad_alloc &){
            cerr << "bad_alloc at smallmatrix::allocate_memory()" << endl;
            exit(1);
        }
        

        return;
    }    
    
    inline smallmatrix& operator +=(const smallmatrix &mat){
        for (int i=0; i < row; i++){
            for (int j=0; j < column; j++){
                element[i][j] += mat.element[i][j];
            }
        }

		return *this;
	}
    
    void matprod1 (smallmatrix &prod,
                   const double &x, const double &y, const double &z){
        //        prod.allocate_memory(row,column);
        for(int i=0; i < row; i++){
            prod.element[i][0] = -element[i][1]*z + element[i][2]*y; 
            prod.element[i][1] = -element[i][2]*x + element[i][0]*z; 
            prod.element[i][2] = -element[i][0]*y + element[i][1]*x; 
        }

    }
 
    void matprod2 (smallmatrix &prod,
                   const double &x, const double &y, const double &z){
        for(int i=0; i < row; i++){
            prod.element[i][0] = element[i][0]*x - element[i][2]*z;
            prod.element[i][1] = element[i][0]*y + element[i][1]*x;
            prod.element[i][2] = element[i][0]*z + element[i][2]*x;
            prod.element[i][3] = element[i][1]*z + element[i][2]*y;
            prod.element[i][4] = element[i][1]*y - element[i][2]*z;
        }
    }
    
    void mat_vector_prod (smallmatrix &prod,
                          const double &x, const double &y, const double &z){
        for (int k=0; k < column ; k++){
            prod.element[0][k] =  element[2][k]*y - element[1][k]*z;
            prod.element[1][k] =  element[0][k]*z - element[2][k]*x;             
            prod.element[2][k] =  element[1][k]*x - element[0][k]*y;
        }
    }

    void multiplyVector(const double *vec, double *product){
        for (int j=0; j < row; j++){
            product[j] = 0.0;
            for (int i=0; i < column; i++){
                product[j] += element[j][i]*vec[i];
            }
        }
    }
        
    void set_zero(){
        for (int i=0; i < row; i++){
            for (int j=0; j < column; j++){
                element[i][j] = 0;
            }
        }
    }

};
#endif	


