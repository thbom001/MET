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


using namespace std;

int main(){
   Grid *g = (Grid *) nullptr;
   AlbersData adata;
   adata.name                 = "Albers Equal Area";
   adata.std_parallel_1       = -36.0;
   adata.std_parallel_2       = -18.0;
   adata.lon_orient           = 132.0;
   adata.lat_centre           = 0.0;
   adata.nx                   = 1959;
   adata.ny                   = 1706;
   adata.ll_x                 = -2198800.;
   adata.ll_y                 = -4998800.;
   adata.dx_m                 = 2400.;
   adata.dy_m                 = 2400.;
   adata.semi_major_axis_km   = 6378137.;
   adata.semi_minor_axis_km   = 6356752.31414036; 
   //adata.eccentricity         = 0.0818191910428;

   g = new Grid ( adata );

   double x, y, lat, lon;
   x=-2198800.;
   y=-4998800.;
   g->xy_to_latlon(x, y, lat, lon);
   printf("(%f, %f) -> (%f, %f)\n", x, y, lat, lon);

   return 0;
}
