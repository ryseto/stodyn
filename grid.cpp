/*
 *  grid.cpp
 *  CCN_3D
 *
 *  Created by seto on 09/08/18.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "grid.h"
#include <algorithm>

Grid::Grid(){
}

Grid::~Grid(){
    if (allocated){
        for (int gx = 0; gx < gx_max; ++gx){
            for (int gy = 0; gy < gy_max ; ++gy){
                DELETE(vcell[gx][gy]);	
                DELETE(neighbor_cell[gx][gy]);	
            }
            DELETE(vcell[gx]);
            DELETE(neighbor_cell[gx]);	
        }
        DELETE(vcell);
        DELETE(neighbor_cell);	
    }
}

void Grid::init
(int num_of_particle_, double L[3], double grid_size)
{
    allocated= true;
	num_of_particle = num_of_particle_;
	h = grid_size;
	gx_max = (int)( L[0]/h );
	gy_max = (int)( L[1]/h );
	gz_max = (int)( L[2]/h );

    cerr << gx_max << endl;
    
	vcell = new vector<int> ** [gx_max];
	neighbor_cell = new	vector<GridPoint> ** [gx_max];
	for (int gx = 0; gx < gx_max; ++gx){
		vcell[gx] = new vector<int> * [gy_max];
		neighbor_cell[gx] = new	vector<GridPoint> * [gy_max];
		for (int gy = 0; gy < gy_max; ++gy){	
			vcell[gx][gy] = new vector<int> [gz_max];
			neighbor_cell[gx][gy] = new	vector<GridPoint> [gz_max];
		}
	} 
	gl.resize(num_of_particle);
	for (int gx = 0; gx < gx_max; ++gx){
		for (int gy = 0; gy < gy_max; ++gy){	
			for (int gz = 0; gz < gz_max; ++gz){
				////////////////////////////////////////
				GridPoint gp;
				for(gp.x = gx-1; gp.x <= gx+1; ++gp.x){	
					for(gp.y = gy-1; gp.y <= gy+1; ++gp.y){	
						for(gp.z = gz-1; gp.z <= gz+1; ++gp.z){
							if ( gp.z >= 0 && gp.z < gz_max){
								GridPoint gp_tmp = gp;
								if (gp_tmp.x == -1) gp_tmp.x = gx_max - 1;
								else if (gp_tmp.x == gx_max) gp_tmp.x = 0;
								if (gp_tmp.y == -1)	gp_tmp.y = gy_max - 1;
								else if (gp_tmp.y == gy_max) gp_tmp.y = 0;
								neighbor_cell[gx][gy][gz].push_back( gp_tmp );
							}
						}
					}
				}
			}
		}
	}

}

void Grid::remake(vector<Particle *> &particle){
    clear_vcell();
    //	for (int i = 0; i < num_of_particle; ++i ){
    //		vcell[gl[i].x][gl[i].y][gl[i].z].clear();
    //	}	
	for (int i = 0; i < num_of_particle; ++i ){
        gl[i] = p_to_grid( particle[i]->p ); 
		vcell[gl[i].x][gl[i].y][gl[i].z].push_back(i);
	}
}

GridPoint Grid::p_to_grid(const vec3d &p){
    gp_tmp.x = (int)(p.x/h);
	gp_tmp.y = (int)(p.y/h);
	gp_tmp.z = (int)(p.z/h);
    if (gp_tmp.x == gx_max ){
        gp_tmp.x --;
    }
    if (gp_tmp.y == gy_max ){
        gp_tmp.y --;
    }
    if (gp_tmp.z == gz_max ){
        gp_tmp.z --;
    }
	return gp_tmp;  
}

void Grid::entry(const vec3d &p, int i){
	GridPoint gp = p_to_grid(p);
	vcell[gp.x][gp.y][gp.z].push_back(i);
	gl[i] = gp; 
}

void Grid::get_neighbor_list(const vec3d &p, vector<int> &neighbor){
	//GridPoint gp;
    gp_tmp = p_to_grid(p);
	foreach(vector<GridPoint>, neighbor_cell[gp_tmp.x][gp_tmp.y][gp_tmp.z], gp){
		foreach(vector<int>, vcell[(*gp).x][(*gp).y][(*gp).z], iter){
			neighbor.push_back(*iter);
		}
	}
}

void Grid::get_neighbor_list_pointer(const vec3d &p, vector< vector<int>* > &p_neighbor){
	p_neighbor.clear();
    gp_tmp = p_to_grid(p);
    
    foreach(vector<GridPoint>, neighbor_cell[gp_tmp.x][gp_tmp.y][gp_tmp.z], gp){
		p_neighbor.push_back( &(vcell[(*gp).x][(*gp).y][(*gp).z]));
	}
}

