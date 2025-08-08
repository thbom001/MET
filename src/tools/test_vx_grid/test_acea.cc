#include <cstdio>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "vx_math.h"
#include "vx_util.h"
#include "vx_log.h"
#include "acea_grid.h"
#include "laea_grid.h"


using namespace std;

int main(){
	std::cout << "Hello World!\n";
	
	Grid *g = (Grid *) nullptr;
	AlbersData adata;
	LambertData ldata;
	adata.std_parallel_1 = -36.0;
	adata.std_parallel_2 = -18.0;
	adata.lon_orient     = 132.0;
   adata.lat_centre     = 0.0;
	adata.nx					= 1959;
	adata.ny					= 1706;
	adata.ll_x				= -2198800.;
	adata.ll_y				= -4998800.;
	adata.dx_m				= 2400.;
	adata.dy_m				= 2400.;
	adata.eccentricity	= 0.0;

	std::cout << "Hello World!\n";
	g = new Grid ( adata );
	std::cout << "Hello World!\n";

	double x, y, lat, lon;
	x=100.;
	y=-3000.;
	g->xy_to_latlon(x, y, lat, lon);
	printf("(%f, %f) -> (%f, %f)\n", x, y, lat, lon);
	

	return 0;
}
