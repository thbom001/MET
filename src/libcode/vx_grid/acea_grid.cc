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

static void   reduce(double &);
static double albers_segment_area(double u0, double v0, double u1, double v1);
static double snyder_q_fcn(double lat, double ecc);
static double snyder_m_fcn(double lat, double ecc);
static double snyder_beta_fcn(double q, double ecc);

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

Name.clear();
Std_parallel_1 = 0.0;
Std_parallel_2 = 0.0;
Lon_orient     = 0.0;
Lat_centre     = 0.0;
Nx             = 0;
Ny             = 0;
Ll_x           = 0.0;
Ll_y           = 0.0;
Dx_m           = 0.0;
Dy_m           = 0.0;

memset(&Data, 0, sizeof(Data));

return;

}


////////////////////////////////////////////////////////////////////////


AlbersGrid::AlbersGrid(const AlbersData & data)

{

clear();

Name        = data.name;

Lon_orient  = data.lon_orient;
reduce(Lon_orient);
Lat_centre  = data.lat_centre;
Nx          = data.nx;
Ny          = data.ny;
Ll_x        = data.ll_x;
Ll_y        = data.ll_y;
Dx_m        = data.dx_m;
Dy_m        = data.dy_m;

Data = data;

   //
   //  Done
   //

}


////////////////////////////////////////////////////////////////////////


// double AlbersGrid::f(double lat) const
// 
// {
// 
// return acea_func(lat, Cone, IsNorthHemisphere);
// 
// }


////////////////////////////////////////////////////////////////////////


void AlbersGrid::latlon_to_xy(double lat, double lon, double & x, double & y) const

{
    // Project (lat, lon) geographical coordinates onto the Albers Conic
    // Equal Area map.
    //
    // Formulae are from Snyder (https://pubs.usgs.gov/publication/pp1395).
    //
    // Input variables:
    // lat:     latitude in degrees North.
    // lon:     longitude in degrees East.
    //
    // Output variables:
    // x:       Projected cartesian X coordinate (units: metres).
    // y:       Projected cartesion Y coordinate (units: metres).


if (is_eq(Data.eccentricity, 0.0)) {
   // Spherical Albers conic equal area formulae (Snyder, p. 100).
   double theta, n, C, rho_0, rho;
   reduce(lon);                                                     // Ensure lon is in the
                                                                    // range [-180., 180)
   n = (sind(Data.std_parallel_1) + sind(Data.std_parallel_2))/2;   // Snyder Eq. 14-6.
   C = pow(cosd(Data.std_parallel_1), 2) +
       2*n*sind(Data.std_parallel_1);                               // Snyder Eq. 14-5.
   theta = n*(lon - Data.lon_orient);                               // Snyder Eq. 14-4.
   rho_0 = (earth_radius_km * 1000) * 
       sqrt((C - 2*n*sind(Data.lat_centre)))/n;                     // Snyder Eq. 14-3a.
   rho = (earth_radius_km * 1000) *
           sqrt((C - 2*n*sind(lat)))/n;                             // Snyder Eq. 14-3.

   x = rho*sind(theta);
   y = rho_0 - rho*cos(theta);
}
else {
   // Ellipsoidal Albers conic equal area formulae. From Snyder, p. 101.
   double theta, n, C, rho_0, rho, q, q0, q1, q2, m1, m2;
   q0    = snyder_q_fcn(Data.lat_centre, Data.eccentricity);
   q1    = snyder_q_fcn(Data.std_parallel_1, Data.eccentricity);
   q2    = snyder_q_fcn(Data.std_parallel_2, Data.eccentricity);
   q     = snyder_q_fcn(lat, Data.eccentricity);
   m1    = snyder_m_fcn(Data.std_parallel_1, Data.eccentricity);
   m2    = snyder_m_fcn(Data.std_parallel_2, Data.eccentricity);
   n     = (pow(m1,2)-pow(m2,2))/(q2-q1);
   C     = pow(m1,2)+n*q1;
   rho_0 = (earth_radius_km*1000)*sqrt(C-n*q0)/n;
   theta = n*(lon-Data.lon_orient);
   rho   = (earth_radius_km*1000)*(C-n*q);
}

return;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::xy_to_latlon(double x, double y, double & lat, double & lon) const

{
    // Project (x, y) projected metres on the Albers Conic Equal Area map to
    // geographical (lon, lat) coordinates.
    //
    // Formulae are from Snyder (https://pubs.usgs.gov/publication/pp1395).
    //
    // Input variables:
    // lat:     latitude in degrees North.
    // lon:     longitude in degrees East.
    //
    // Output variables:
    // x:       Projected cartesian X coordinate (units: metres).
    // y:       Projected cartesion Y coordinate (units: metres).

if (is_eq(Data.eccentricity, 0.0)) {
   double theta, n, C, rho_0, rho;
   // Spherical Albers conic equal area inverse formulae (Snyder, p. 101).
   n = (sind(Data.std_parallel_1) + sind(Data.std_parallel_2))/2; // Snyder Eq. 14-6.
   C = cosd(Data.std_parallel_1)*cosd(Data.std_parallel_1) +
      2*n*sind(Data.std_parallel_1);                              // Snyder Eq. 14-5.
   rho_0 = (earth_radius_km * 1000) * 
      sqrt((C - 2*n*sind(Data.lat_centre)))/n;                    // Snyder Eq. 14-3a.
   rho = sqrt(pow(x, 2) + pow((rho_0-y), 2));                     // Snyder Eq. 14-10.
   theta = atand(x/(rho_0 - y));                                  // Snyder Eq. 14-11.

   lat = asind((C-pow((rho*n/(earth_radius_km * 1000)), 2))
      /(2*n));                                                    // Snyder Eq. 14-8.
   lon = Data.lon_orient + theta/n;                               // Snyder Eq. 14-9.
   reduce(lon);
} else {
    // Ellipsoidal Albers conic equal area formulae (p. 102 of Snyder).
   double n, C, rho_0, rho, theta, m1, m2, q, q0, q1, q2;

   m1    = snyder_m_fcn(Data.std_parallel_1, Data.eccentricity);
   m2    = snyder_m_fcn(Data.std_parallel_2, Data.eccentricity);
   q0    = snyder_q_fcn(Data.lat_centre, Data.eccentricity);
   q1    = snyder_q_fcn(Data.std_parallel_1, Data.eccentricity);
   q2    = snyder_q_fcn(Data.std_parallel_2, Data.eccentricity);
   n     = (pow(m1,2)-pow(m2,2))/(q2-q1);                         // Snyder Eq. 14-14.
   C     = pow(m1,2)+n*q1;                                        // Snyder Eq. 14-13.
   rho_0 = (earth_radius_km*1000)*sqrt(C-n*q0)/n;                 // Snyder Eq. 14-12a.

   theta = atand(x/(rho_0-y));                                    // Snyder Eq. 14-11.
   rho   = sqrt(pow(x,2)+pow((rho_0-y),2));                       // Snyder Eq. 14-10.
   q     = (C-pow(rho,2)*pow(n,2)/
            pow((earth_radius_km*1000),2))/n;                     // Snyder Eq. 14-19.

   lat   = snyder_beta_fcn(q, Data.eccentricity);                 // Snyder Eq. 14-18 and 14-21.
   lon   = Data.lon_orient + theta/n;                             // Snyder Eq. 14-9.

}

return;

}


////////////////////////////////////////////////////////////////////////


double AlbersGrid::calc_area(int x, int y) const

{

 double u[4], v[4];
 double sum;
 
 
 xy_to_uv(x - 0.5, y - 0.5, u[0], v[0]);  //  lower left
 xy_to_uv(x + 0.5, y - 0.5, u[1], v[1]);  //  lower right
 xy_to_uv(x + 0.5, y + 0.5, u[2], v[2]);  //  upper right
 xy_to_uv(x - 0.5, y + 0.5, u[3], v[3]);  //  upper left
 
 sum = uv_closedpolyline_area(u, v, 4);
 
 sum *= (earth_radius_km*1000)*(earth_radius_km*1000);

return sum;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::xy_to_uv(double x, double y, double &u, double &v) const

{

aff.reverse(x, y, u, v);

return;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::uv_to_xy(double u, double v, double &x, double &y) const

{

aff.forward(u, v, x, y);

return;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::dump(ostream & out, int depth) const

{

Indent prefix(depth);

out << prefix << "Name         = ";

if ( Name.length() > 0 )  out << '\"' << Name << '\"';
else                      out << "(nul)\n";

//out << prefix << "SpheroidName = ";

//if ( SpheroidName.length() > 0 )  out << '\"' << SpheroidName << '\"';
//else                              out << "(nul)\n";

out << '\n';   //  no prefix

out << prefix << "Projection   = Albers\n";

out << prefix << "Nx           = " << comma_string(Nx) << "\n";

out << prefix << "Ny           = " << comma_string(Ny) << "\n";

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

//a << "SpheroidName: " << SpheroidName << sep;

//snprintf(junk, sizeof(junk), "Lat_LL: %.3f", lat_LL);   a << junk << sep;
//snprintf(junk, sizeof(junk), "Lon_LL: %.3f", lon_LL);   a << junk << sep;

//snprintf(junk, sizeof(junk), "Lat_Pole: %.3f", lat_pole);   a << junk << sep;
//snprintf(junk, sizeof(junk), "Lon_Pole: %.3f", lon_pole);   a << junk << sep;

   //
   //  done
   //

return a;

}


////////////////////////////////////////////////////////////////////////


void AlbersGrid::deserialize(const StringArray &)

{

mlog << Error << "\nAlbersGrid::deserialize(const StringArray &) -> "
     << "not yet implemented\n\n";

exit ( 1 );

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


GridInfo AlbersGrid::info() const

{

GridInfo i;

i.set(Data);

return i;

}


////////////////////////////////////////////////////////////////////////


double AlbersGrid::rot_grid_to_earth(int x, int y) const

{

//
// I don't think the grid to Earth transformation is a simple rotation,
// because the Albers equal area projection is not conformal.
//

return 0.0;

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


   //
   //  Code for misc functions
   //


////////////////////////////////////////////////////////////////////////


void reduce(double & angle)

{
    // Convert an angle in degrees to lie within the range [-180, 180)

angle -= 360.0*floor( (angle/360.0) + 0.5 );

return;

}


////////////////////////////////////////////////////////////////////////


double snyder_q_fcn(double lat, double ecc)

{
   // Compute "q" using Equation 3-12 (p. 101) in Snyder.
   //
   // Input variables:
   // lat:     latitude in degrees North.
   // ecc:     eccentricity of the ellipsoid.
   //
   // Return:
   // q:       Snyder Equation 3-12.
   double q;

   q = (1-pow(ecc,2))*(sind(lat)/(1-pow(ecc,2)*pow(sind(lat), 2))
      - (1/(2*ecc))*log((1-ecc*sind(lat))/(1+ecc*sind(lat))));

   return q;
}


////////////////////////////////////////////////////////////////////////


double snyder_m_fcn(double lat, double ecc)

{
   // Compute "m" using Equation 14-15 (p. 101) in Snyder.
   //
   // Input variables:
   // lat:     latitude in degrees North.
   // ecc:     eccentricity of the ellipsoid.
   //
   // Return:
   // m:       Snyder Equation 14-15.
   double m;

   m = cosd(lat)/sqrt(1-pow(ecc,2)*pow(sind(lat),2));

   return m;
}


////////////////////////////////////////////////////////////////////////


double snyder_beta_fcn(double q, double ecc)

{
   // Compute latitude using Equations 14-18 and 14-21 (p. 102) in Snyder.
   //
   // Input variables:
   // q:			Computed previously (Snyder Equation 14-19).
   // ecc:     eccentricity of the ellipsoid.
   //
   // Return:
   // lat:     Snyder Equation 14-18 and 14-21.
   double beta, lat;

   beta = asind(q/(1-((1-pow(ecc,2))/(2*ecc))*log((1-ecc)/(1+ecc))));
   lat = beta +
      (pow(ecc,2)/3 + 31*pow(ecc,4)/180 + 517*pow(ecc,6)/5040)*sind(2*beta) +
      (23*pow(ecc,4)/360 + 251*pow(ecc,6)/3780)*sind(4*beta) +
      (761*pow(ecc,6)/45360)*sind(6*beta);

   return lat;
}


////////////////////////////////////////////////////////////////////////


double albers_segment_area(double u0, double v0, double u1, double v1)

{
	// This computes the area of the triangular segment defined by two points on
	// the map with coordinates (u0, u1) and (v0, v1) and the origin (0, 0). The
	// computation takes the cross product of these two vectors and then
	// computes half the magnitude of this cross product. See texts on the
	// geometrical interpretation of the cross product for why this works.

double answer;

answer = 0.5*(u0*v1 - u1*v0);

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


double AlbersGrid::uv_closedpolyline_area(const double *u, const double *v, int n) const

{

int k;
double sum;


sum = 0.0;

for (int j=0; j<n; ++j)  {

   k = (j + 1)%n;

   sum += albers_segment_area(u[j], v[j], u[k], v[k]);

}   //  for j

sum = fabs(sum);

return sum;

}


////////////////////////////////////////////////////////////////////////


double AlbersGrid::xy_closedpolyline_area(const double *x, const double *y, int n) const

{

double sum;
double *u = (double *) nullptr;
double *v = (double *) nullptr;

u = new double [n];
v = new double [n];

if ( !u || !v )  {

   mlog << Error << "\nAlbersGrid::xy_closedpolyline_area() -> "
        << "memory allocation error\n\n";

   exit ( 1 );

}

for (int j=0; j<n; ++j)  {

   xy_to_uv(x[j], y[j], u[j], v[j]);

}

sum = uv_closedpolyline_area(u, v, n);

sum *= earth_radius_km*earth_radius_km;

delete [] u;  u = (double *) nullptr;
delete [] v;  v = (double *) nullptr;

return sum;

}


////////////////////////////////////////////////////////////////////////
