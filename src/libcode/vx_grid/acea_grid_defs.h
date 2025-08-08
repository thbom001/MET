// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2025
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


////////////////////////////////////////////////////////////////////////


#ifndef  __ALBERS_GRID_DEFINITIONS_H__
#define  __ALBERS_GRID_DEFINITIONS_H__


////////////////////////////////////////////////////////////////////////


struct AlbersData {

   const char * name;		//  Not allocated.

   double std_parallel_1;  //  First standard parallel.
   double std_parallel_2;  //  Second standard parallel.
   double lon_orient;      //  Central longitude.
   double lat_centre;      //  Central latitude.
   int nx;
   int ny;
	double ll_x;				// Projected X coordinate - lower left corner (units: m)
	double ll_y;				// Projected Y coordinate - lower left corner (units: m)
   double dx_m;				// Projection (not geographic) X grid spacing (units: m)
   double dy_m;				// Projection (not geographic) Y grid spacing (units: m)
   double eccentricity;		// Eccentricity is 0 for the sphere, and non-zero for an ellipsoid.

   void dump() const;

};


////////////////////////////////////////////////////////////////////////


#endif   /*  __ALBERS_GRID_DEFINITIONS_H__  */


////////////////////////////////////////////////////////////////////////
