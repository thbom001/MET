// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2025
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////

static double acea_func(double lat, double Cone, const bool is_north);
static double acea_der_func(double lat, double Cone, const bool is_north);
static double acea_inv_func(double   r, double Cone, const bool is_north);
static void   reduce(double &);
static double albers_segment_area(double u0, double v0, double u1, double v1, double c);
static double albers_beta(double u0, double delta_u, double v0, double delta_v, double c, double t);

////////////////////////////////////////////////////////////////////////


   //
   //  Code for class AlbersGrid
   //


////////////////////////////////////////////////////////////////////////


AlbersGrid::AlbersGrid()

{

clear();

}


////////////////////////////////////////////////////////////////////////


AlbersGrid::~AlbersGrid()

{

clear();

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::clear()

{

IsNorthHemisphere = true;

Lat_LL = 0.0;               // Latitude of lower left corner
Lon_LL = 0.0;               // Longitude of lower left corner

Lon_orient = 0.0;

Alpha = 0.0;

Cone = 0.0;

Bx = 0.0;
By = 0.0;

Nx = 0;
Ny = 0;

Name.clear();

memset(&Data, 0, sizeof(Data));

return;

}


////////////////////////////////////////////////////////////////////////


AlbersGrid::AlbersGrid(const AlbersData & data)

{

clear();

double ratio;

Lat_LL = data.lat_pin;   //  temporarily
Lon_LL = data.lon_pin;   //  temporarily

reduce(Lon_LL);

Lon_orient = data.lon_centre;

reduce(Lon_orient);

Bx = 0.0;
By = 0.0;

Nx = data.nx;
Ny = data.ny;

Name = data.name;

   //
   //  calculate Alpha
   //

ratio = (data.r_km)/(data.d_km);

Alpha = (1.0/acea_der_func(data.std_parallel_1, Cone, IsNorthHemisphere));

Alpha = fabs(Alpha);

Alpha *= ratio;

   //
   //  Calculate Bx, By
   //

double r_pin, theta_pin;

r_pin = acea_func(data.lat_pin, Cone, IsNorthHemisphere);

theta_pin = Cone*(rescale_deg(Lon_orient - data.lon_pin, -180.0, 180.0));

Bx = data.x_pin - Alpha*r_pin*sind(theta_pin);
By = data.y_pin + Alpha*r_pin*cosd(theta_pin);

xy_to_latlon(0.0, 0.0, Lat_LL, Lon_LL);

reduce(Lon_LL);

set_so2(data.so2_angle);

Data = data;

   //
   //  Done
   //

}


////////////////////////////////////////////////////////////////////////


double AlbersGrid::f(double lat) const

{

return acea_func(lat, Cone, IsNorthHemisphere);

}


////////////////////////////////////////////////////////////////////////


double AlbersGrid::df(double lat) const

{

return acea_der_func(lat, Cone, IsNorthHemisphere);

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::latlon_to_xy(double lat, double lon, double & x, double & y) const

{
    // Project (lat, lon) geographical coordinates onto the Albers Conic
    // Equal Area map.
    //
    // Initially we use the spherical formulae given in Snyder
    // (https://pubs.usgs.gov/publication/pp1395), but we will probably need to add
    // support fol the ellipsoidal formulae.
    //
    // Input variables:
    // lat:     latitude in degrees North.
    // lon:     longitude in degrees East.
    //
    // Output variables:
    // x:       Projected cartesian X coordinate (units: metres).
    // y:       Projected cartesion Y coordinate (units: metres).

double r, theta, n, C, rho_0, rho;

if (is_eq(Data.eccentricity, 0.0)) {
    // Spherical Albers conic equal area formulae (Snyder, p. 100).
    n = (sind(Data.std_parallel_1) + sind(Data.std_parallel_2))/2;   // Snyder Eq. 14-6.
    C = cosd(Data.std_parallel_1)*cosd(Data.std_parallel_1) +
        2*n*sind(Data.std_parallel_1);                              	// Snyder Eq. 14-5.
    theta = n*(lon - Data.lon_centre);               		// Snyder Eq. 14-4.
    rho_0 = earth_radius_km * 
        sqrt((C - 2*n*sind(Data.lat_centre)))/n;            		// Snyder Eq. 14-3a.
    rho = earth_radius_km*
            sqrt((C - 2*n*sind(lat)))/n;                    		// Snyder Eq. 14-3.

    x = rho*sind(theta);
    y = rho_0 - rho*cos(theta);
}
else {
    // Ellipsoidal Albers conic equal area formulae.
    // Still to be implemented.
}

reduce(lon);

r = acea_func(lat, Cone, IsNorthHemisphere);

theta = Cone*(Lon_orient - lon);

x = Bx + Alpha*r*sind(theta);

y = By - Alpha*r*cosd(theta);

return;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::xy_to_latlon(double x, double y, double & lat, double & lon) const

{

double r, theta;

x = (x - Bx)/(Alpha);
y = (y - By)/(Alpha);

r = sqrt( x*x + y*y );

lat = acea_inv_func(r, Cone, IsNorthHemisphere);

if ( fabs(r) < 1.0e-5 )  theta = 0.0;
else                     theta = atan2d(x, -y);   //  NOT atan2d(y, x);

lon = Lon_orient - theta/(Cone);

reduce(lon);

return;

}


////////////////////////////////////////////////////////////////////////


double AlbersGrid::calc_area(int x, int y) const

{

double u[4], v[4];
double sum;


// xy_to_uv(x - 0.5, y - 0.5, u[0], v[0]);  //  lower left
// xy_to_uv(x + 0.5, y - 0.5, u[1], v[1]);  //  lower right
// xy_to_uv(x + 0.5, y + 0.5, u[2], v[2]);  //  upper right
// xy_to_uv(x - 0.5, y + 0.5, u[3], v[3]);  //  upper left


xy_to_uv(x      , y      , u[0], v[0]);  //  lower left
xy_to_uv(x + 1.0, y      , u[1], v[1]);  //  lower right
xy_to_uv(x + 1.0, y + 1.0, u[2], v[2]);  //  upper right
xy_to_uv(x      , y + 1.0, u[3], v[3]);  //  upper left


sum = uv_closedpolyline_area(u, v, 4);

sum *= earth_radius_km*earth_radius_km;

return sum;

}


////////////////////////////////////////////////////////////////////////


int AlbersGrid::nx() const

{

return Nx;

}


////////////////////////////////////////////////////////////////////////


int AlbersGrid::ny() const

{

return Ny;

}


////////////////////////////////////////////////////////////////////////


ConcatString AlbersGrid::name() const

{

return Name;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::dump(ostream & out, int depth) const

{

Indent prefix(depth);



out << prefix << "Name       = ";

if ( Name.length() > 0 )  out << '\"' << Name << '\"';
else                      out << "(nul)\n";

out << '\n';

out << prefix << "Projection = Albers Conic Equal Area\n";

out << prefix << "\n";

out << prefix << "Lat_LL     = " << Lat_LL << "\n";
out << prefix << "Lon_LL     = " << Lon_LL << "\n";

out << prefix << "\n";

out << prefix << "Alpha      = " << Alpha << "\n";
out << prefix << "Cone       = " << Cone  << "\n";

out << prefix << "\n";

out << prefix << "Bx         = " << Bx << "\n";
out << prefix << "By         = " << By << "\n";

out << prefix << "\n";

out << prefix << "Nx         = " << Nx << "\n";
out << prefix << "Ny         = " << Ny << "\n";


   //
   //  done
   //

out.flush();

return;

}


////////////////////////////////////////////////////////////////////////


ConcatString AlbersGrid::serialize(const char *sep) const

{

ConcatString a;
char junk[256];

a << "Projection: Albers Conic Equal Area" << sep;

a << "Nx: " << Nx << sep;
a << "Ny: " << Ny << sep;

snprintf(junk, sizeof(junk), "Lat_LL: %.3f", Lat_LL);   a << junk << sep;
snprintf(junk, sizeof(junk), "Lon_LL: %.3f", Lon_LL);   a << junk << sep;

snprintf(junk, sizeof(junk), "Lon_orient: %.3f", Lon_orient);   a << junk << sep;

snprintf(junk, sizeof(junk), "Alpha: %.3f", Alpha);   a << junk << sep;

snprintf(junk, sizeof(junk), "Cone: %.3f", Cone);   a << junk << sep;

snprintf(junk, sizeof(junk), "Bx: %.4f", Bx);   a << junk << sep;
snprintf(junk, sizeof(junk), "By: %.4f", By);   a << junk;

   //
   //  done
   //

return a;

}


////////////////////////////////////////////////////////////////////////


GridInfo AlbersGrid::info() const

{

GridInfo i;

i.set(Data);

return i;

}


////////////////////////////////////////////////////////////////////////


double AlbersGrid::rot_grid_to_earth(int x, int y) const

{

double lat, lon, angle;
double diff, hemi;


xy_to_latlon((double) x, (double) y, lat, lon);

diff = Lon_orient - lon;

// Figure out if the grid is in the northern or southern hemisphere
// by checking whether the first latitude (p1_deg -> Phi1_radians)
// is greater than zero
// NH -> hemi = 1, SH -> hemi = -1
// if(Phi1_radians < 0.0) hemi = -1.0;
// else                   hemi = 1.0;

   //
   //  assume northern hemisphere
   //

hemi = 1.0;

angle = diff*Cone*hemi;

return angle;

}


////////////////////////////////////////////////////////////////////////


bool AlbersGrid::wrap_lon() const

{

return false;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::shift_right(int N)

{

if ( N == 0 )  return;

mlog << Error << "\nAlbersGrid::shift_right(int) -> "
     << "shifting is not allowed for non-global grids\n\n";

exit ( 1 );

}


////////////////////////////////////////////////////////////////////////


GridRep * AlbersGrid::copy() const

{

AlbersGrid * p = new AlbersGrid (Data);

p->Name = Name;

return p;

}


////////////////////////////////////////////////////////////////////////


   //
   //  Code for struct AlbersData
   //


////////////////////////////////////////////////////////////////////////

/*
AlbersData::AlbersData()

{

hemisphere = 'N';

}
*/


////////////////////////////////////////////////////////////////////////


   //
   //  Code for misc functions
   //


////////////////////////////////////////////////////////////////////////


double acea_func(double lat, double Cone, const bool is_north)

{

double r;

r = tand(45.0 - 0.5*lat);

r = pow(r, Cone);

return r;

}


////////////////////////////////////////////////////////////////////////


double acea_inv_func(double r, double Cone, const bool is_north)

{

double lat;

lat = 90.0 - 2.0*atand(pow(r, 1.0/Cone));


return lat;

}


////////////////////////////////////////////////////////////////////////


double acea_der_func(double lat, double Cone, const bool is_north)

{

double a;

a = -(Cone/cosd(lat))*acea_func(lat, Cone, is_north);


return a;

}


////////////////////////////////////////////////////////////////////////


void reduce(double & angle)

{
    // Convert an angle in degrees to lie within the range [-180, 180)

angle -= 360.0*floor( (angle/360.0) + 0.5 );

return;

}


////////////////////////////////////////////////////////////////////////


double albers_segment_area(double u0, double v0, double u1, double v1, double c)

{

int i, j, k, n;
double rom, denom, h, delta_u, delta_v;
double trap, t[15], left, right, sum;
double test = 0.0;
const double a = 0.0, b = 1.0;
const double tol = 1.0e-6;

delta_u = u1 - u0;
delta_v = v1 - v0;

i = 0;
n = 2;

h = (b - a)/n;

sum = albers_beta(u0, delta_u, v0, delta_v, c, a) + albers_beta(u0, delta_u, v0, delta_v, c, b);

t[0] = trap = (h/2.0)*sum + h*albers_beta(u0, delta_u, v0, delta_v, c, a + h);

do {

   ++i;

   n *= 2;

   h = (b - a)/n;

   sum = 0.0;

   for (j=1; j<n; j+=2)   sum += albers_beta(u0, delta_u, v0, delta_v, c, a + j*h);

   trap = 0.5*trap + h*sum;

   left = trap;

   for (k=1; k<=i; ++k)  {

      denom = pow(4.0, (double) k) - 1.0;

      right = left + (left - t[k-1])/denom;

      test = 2.0*(left - t[k-1]);

      t[k-1] = left;

      left = right;

   }

   t[i] = left;

}  while ( (fabs(test) >= tol) && (i <= 14) );

if ( i >= 14 )  {

   mlog << Error << "\nlambert_segment_area() -> "
        << "array bounds error\n\n";

   exit ( 1 );

}

rom = t[i];

rom *= (2.0/c)*(u0*v1 - u1*v0);

return rom;

}


////////////////////////////////////////////////////////////////////////


double albers_beta(double u0, double delta_u, double v0, double delta_v, double c, double t)

{

double answer;
double u, v, r2, e_top, e_bot;

u = u0 + t*delta_u;
v = v0 + t*delta_v;

r2 = u*u + v*v;

e_bot = 1.0/c;

e_top = e_bot - 1.0;

answer = pow(r2, e_top)/(1.0 + pow(r2, e_bot));

return answer;

}


////////////////////////////////////////////////////////////////////////


Grid::Grid(const AlbersData & data)

{

init_from_scratch();

set(data);


}


////////////////////////////////////////////////////////////////////////


void Grid::set(const AlbersData & data)

{

clear();

rep = new AlbersGrid (data);

if ( !rep )  {

   mlog << Error << "\nGrid::set(const AlbersData &) -> "
        << "memory allocation error\n\n";

   exit ( 1 );

}

}


////////////////////////////////////////////////////////////////////////
