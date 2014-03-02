/*
 *  my_utilities.h
 *  stodyn
 *
 *  Created by seto on 10/05/03.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef my_utilities_h
#define my_utilities_h 1
#include <fstream>
#include <iomanip>
#include <sys/time.h>
//#define DELETE(x) if(x){delete [] x; x = NULL;}
using namespace std;

/* math tools */
inline double pow_10(int index){
	double val = 1.0;
	if ( index > 0 ){
		while ((index --) > 0){
			val *= 10;
		}
	} else {
		while ((index ++) < 0){
			val /= 10;
		}
	}
	return val ;
}

//
//void samplecode_ImportGRMdata(int argc, char** argv){
//    /*
//     * Grand resistance matrix is huge matrix.
//     * If you want to remain the matrix data in file,
//     * txt file format can be too slwo.
//     * We can keep it as binary format.
//     */
//    fstream import_grm;
//    import_grm.open(argv[2], ios::in | ios::binary);
//    double *tmp;
//    int n = 1024;
//    int n11=n*11;
//    int datasize = n11*n11*sizeof(tmp[0]);
//    tmp = new double [n11*n11];
//    import_grm.read((char *) tmp, datasize);    
//    import_grm.close();
//}


#endif





