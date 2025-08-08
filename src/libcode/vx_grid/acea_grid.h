// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2024
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


////////////////////////////////////////////////////////////////////////


#ifndef  __ALBERS_GRID_H__
#define  __ALBERS_GRID_H__


////////////////////////////////////////////////////////////////////////


#include "acea_grid_defs.h"
#include "grid_base.h"


////////////////////////////////////////////////////////////////////////


class AlbersGrid : public GridRep {

      friend class Grid;

   private:

      AlbersGrid();
     ~AlbersGrid();
      AlbersGrid(const AlbersData &);

      void clear();

         //
         //
         //

      ConcatString Name;

      double Std_parallel_1;
		double Std_parallel_2;
		double Lon_orient;
		double Lat_centre;
      int Nx;
      int Ny;
		double Ll_x;
		double Ll_y;
		double Dx_m;
		double Dy_m;

      double  f(double) const;
      double df(double) const;

      AlbersData Data;

         //

      bool is_north() const;
      bool is_south() const;

         //
         //  grid interface
         //

      void xy_to_uv(double x, double y, double & u, double & v) const;
      void uv_to_xy(double u, double v, double & x, double & y) const;

      void latlon_to_xy(double lat, double lon, double & x, double & y) const;
      void xy_to_latlon(double x, double y, double & lat, double & lon) const;

      double calc_area(int x, int y) const;

      int nx() const;
      int ny() const;

      ConcatString name() const;

      void dump(std::ostream &, int = 0) const;

      ConcatString serialize(const char *sep=" ") const;

      GridInfo info () const;

      double rot_grid_to_earth(int x, int y) const;

      bool wrap_lon() const;

      void shift_right(int);

      GridRep * copy() const;

      double scale_km() const;

};


////////////////////////////////////////////////////////////////////////


//inline double AlbersGrid::scale_km() const { return Data.d_km; }


////////////////////////////////////////////////////////////////////////


extern Grid create_oriented_ac(bool is_north_projection,
                               double lat_cen,  double lon_cen,
                               double lat_prev, double lon_prev,
                               double d_km, double r_km,
                               int nx, int ny,
                               double bearing);


////////////////////////////////////////////////////////////////////////


#endif   //  __ALBERS_GRID_H__


////////////////////////////////////////////////////////////////////////
