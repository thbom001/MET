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
   std::cout << "Hello World 0!\n";
   
   Grid *g = (Grid *) nullptr;
   AlbersData adata;
   adata.std_parallel_1    = -36.0;
   adata.std_parallel_2    = -18.0;
   adata.lon_orient        = 132.0;
   adata.lat_centre        = 0.0;
   adata.nx                = 1959;
   adata.ny                = 1706;
   adata.ll_x              = -2198800.;
   adata.ll_y              = -4998800.;
   adata.dx_m              = 2400.;
   adata.dy_m              = 2400.;
   adata.semi_major_axis   = 6378137.;
   adata.semi_minor_axis   = 6356752.31414036; 
   adata.eccentricity      = 0.0818191910428;
   //adata.eccentricity    = 0;

   std::cout << "Hello World 1!\n";
   g = new Grid ( adata );
   std::cout << "Hello World 2!\n";

   double x, y, lat, lon;
   x=-2198800.;
   y=-4998800.;
   std::cout << "Hello World 3!\n";
   g->xy_to_latlon(x, y, lat, lon);
   std::cout << "Hello World 4!\n";
   printf("(%f, %f) -> (%f, %f)\n", x, y, lat, lon);
   std::cout << "Hello World 5!\n";

   return 0;
}
