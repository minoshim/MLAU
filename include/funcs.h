/* Copyright 2020 Takashi Minoshima */

/* This file is part of MLAU. */

/* MLAU is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* MLAU is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with MLAU.  If not, see <https://www.gnu.org/licenses/>. */

#ifndef _FUNCS_H_
#define _FUNCS_H_

#include <math.h>

inline double max(double a, double b)
{
  return( (a >= b)?a:b );
}

inline double min(double a, double b)
{
  return( (a >= b)?b:a );
}

inline double minmod(double a, double b)
{
  return(+min(max(a,b),0)
	 +max(min(a,b),0));
}

#endif

