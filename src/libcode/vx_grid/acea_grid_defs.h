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

   const char * name;

   double std_parallel_1;  //  First standard parallel.
   double std_parallel_2;  //  Secnd standard parallel.

   double lat_pin;
   double lon_pin;

   double x_pin;
   double y_pin;

   double lon_centre;      //  Central longitude.
   double lat_centre;      //  Central latitude.

   double d_km;
   double r_km;

   int nx;
   int ny;

   double so2_angle;       //  rotation about pin point

   // LambertData() { so2_angle = 0.0; };

   void dump() const;

};


////////////////////////////////////////////////////////////////////////


#endif   /*  __ALBERS_GRID_DEFINITIONS_H__  */


////////////////////////////////////////////////////////////////////////
