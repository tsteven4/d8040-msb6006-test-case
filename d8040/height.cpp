/*
	Copyright (C) 2008 Alexander Stapff, a.stapff@gmx.de

	Geoid separation code by Oleg Gusev, from data by Peter Dana.
	This code was published by the gpsd project (http://gpsd.berlios.de/)
	under the BSD license.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#include "height.hpp"
#include <cmath>    // for floor
#include <cstdint>  // for int8_t

// Until c++17 we have to define odr-used constexpr static data members at namespace scope.
#if __cplusplus < 201703L
constexpr double HeightFilter::geoid_grid_deg;
constexpr double HeightFilter::geoid_scale;
constexpr int HeightFilter::geoid_row;
constexpr int HeightFilter::geoid_col;
constexpr int8_t HeightFilter::geoid_delta[];
#endif

double HeightFilter::bilinear(double x1, double y1, double x2, double y2, double x, double y, double z11, double z12, double z21, double z22)
{
  if (y1 == y2 && x1 == x2) {
    return (z11);
  }
  if (y1 == y2 && x1 != x2) {
    return (z22*(x-x1)+z11*(x2-x))/(x2-x1);
  }
  if (x1 == x2 && y1 != y2) {
    return (z22*(y-y1)+z11*(y2-y))/(y2-y1);
  }

  double delta = (y2-y1)*(x2-x1);

  return (z22*(y-y1)*(x-x1)+z12*(y2-y)*(x-x1)+z21*(y-y1)*(x2-x)+z11*(y2-y)*(x2-x))/delta;
}

/* return geoid separation (MSL - WGS84) in meters, given a lat/lot in degrees */
double HeightFilter::wgs84_separation(double lat, double lon)
{
  int ilat = (int)floor((90.0+lat)/geoid_grid_deg);
  int ilon = (int)floor((180.0+lon)/geoid_grid_deg);

  int ilat1 = ilat;
  int ilon1 = ilon;
  int ilat2 = (ilat < geoid_row-1)? ilat+1:ilat;
  int ilon2 = (ilon < geoid_col-1)? ilon+1:ilon;

  return bilinear(
           ilon1*geoid_grid_deg-180.0, ilat1*geoid_grid_deg-90.0,
           ilon2*geoid_grid_deg-180.0, ilat2*geoid_grid_deg-90.0,
           lon, lat,
           (double)geoid_delta[ilon1+ilat1*geoid_col]/geoid_scale,
           (double)geoid_delta[ilon2+ilat1*geoid_col]/geoid_scale,
           (double)geoid_delta[ilon1+ilat2*geoid_col]/geoid_scale,
           (double)geoid_delta[ilon2+ilat2*geoid_col]/geoid_scale
         );
}
