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

#ifndef HEIGHT_H_INCLUDED_
#define HEIGHT_H_INCLUDED_

#include <cstdint>

class HeightFilter
{
public:
  // include static constexpr data member definitions with intializers for grid as private members.
  #include "heightgrid.hpp"

  double bilinear(double x1, double y1, double x2, double y2, double x, double y, double z11, double z12, double z21, double z22);
  double wgs84_separation(double lat, double lon);

};
#endif // HEIGHT_H_INCLUDED_