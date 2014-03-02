/*
 *  grid.h
 *  CCN_3D
 *
 *  Created by seto on 09/08/18.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef grid_h
#define grid_h 1
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <map>
#include "common.h"
#include "vec3d.h"
#include "particle.h"
using namespace std;
const int max_number_point =20000;
//extern int gy_max, gx_max;
class Particle;

class Grid{
	double h;
	vector<GridPoint> gl;
	vector<int> ***vcell;
	vector<GridPoint> ***neighbor_cell;
	int num_of_particle;
	GridPoint gp_tmp;
	int gx_max;
	int gy_max;
	int gz_max;
    bool allocated;
public:
	Grid();
	~Grid();
    
	void init(int num_of_particle, double L[3] , double grid_size);
	void remake(vector<Particle *> &particle);
	void entry(const vec3d &p, int i);
	GridPoint p_to_grid(const vec3d &p);
	inline vector<int>* particle_in_cell(int x, int y, int z){ return &vcell[x][y][z]; }
	void reset();
	inline int val_gx_max(){ return gx_max; }
	inline int val_gy_max(){ return gy_max; }
	inline int val_gz_max(){ return gz_max; }
	inline int min_gy(int gy){ return ( gy >= 0 ? gy : 0 ) ; }
	inline int max_gy(int gy){ return ( gy <= gy_max ? gy : gy_max ); }
	int gy(int i){ return gl[i].y;}
	int size(int x, int y, int z) { return vcell[x][y][z].size(); }
	void gl_resize(int n){ gl.resize(n); }
    void clear_vcell(){
        for (int i = 0; i < num_of_particle; ++i ){
            vcell[gl[i].x][gl[i].y][gl[i].z].clear();
        }
    }
    inline void remakeOneParticle(int i, vec3d &pos){        
        gl[i] = p_to_grid( pos ); 
		vcell[gl[i].x][gl[i].y][gl[i].z].push_back(i);
    }
	void get_neighbor_list(const vec3d &p, vector<int> &neighbor);
	void get_neighbor_list_pointer(const vec3d &p, vector< vector<int>* > &p_neighbor);
	
};
#endif





