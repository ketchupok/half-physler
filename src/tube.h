/*
  tube.cpp:

  Copyright (C) 2018 - Sebastian Schmutzhard, Alex Hofmann, Gokberk Erdogan
  and Vasileios Chatziioannou

  This file is part of Csound.

  The Csound Library is free software; you can redistribute it
  and/or modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  Csound is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with Csound; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02110-1301 USA
*/


#ifndef SRC_TUBE_H_
#define SRC_TUBE_H_

void grid_init(MYFLT Len, MYFLT dt, MYFLT *dx, int *M, MYFLT *L,
               MYFLT c = 3.4386e+02);

void update_vp_pointers(int M, const MYFLT& dt, const MYFLT& dx, const MYFLT& c,
                        const MYFLT& rho_user, const MYFLT *S,
                        const MYFLT *pold, const MYFLT *vold,
                        MYFLT *pnew, MYFLT *vnew);

MYFLT cross_area(MYFLT radius, MYFLT x, MYFLT prelength, MYFLT slope);

void interpolation(int M, int Mold, MYFLT Lold, MYFLT dx, MYFLT dxold, \
                   csnd::AuxMem<MYFLT> knew, csnd::AuxMem<MYFLT> kold);

void interpolation_pointers(int M, int Mold, MYFLT Lold, MYFLT dx, MYFLT dxold,
                            MYFLT *knew, MYFLT *kold);

#endif  // SRC_TUBE_H_
