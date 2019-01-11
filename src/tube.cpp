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

#include <plugin.h>
#include <math.h>
#include <atomic>
#include <iostream>
#include <fstream>
#include <random>
#include "const.h"
#include "tube.h"


void grid_init(MYFLT Len, MYFLT dt, MYFLT * dx, int * M, MYFLT * L, MYFLT c) {
    // TODO(AH) - replace weird arrays with [subscript]
  // INPUTS: (total) length
  // calculate time step and grid size, taking stability into account
  // return length and grid size (M and dx)

  // In order to satisfy the stability conditions dx/dt>c and M*dx = L
  L[0] = Len;
  dx[0]  = c*dt;
  dx[0]  = dx[0]/0.9;  // Guarantees the first stability condition
  M[0]   = int(floor(L[0] / dx[0]));  // Number of steps has to be an integer
  M[0]   = std::min(M[0], Mmax);
  dx[0]  = L[0] / (M[0]);  // Guarantees both conditions
  }


void update_vp_pointers(int M, const MYFLT& dt, const MYFLT& dx, \
                const MYFLT& c, const MYFLT& rho_user, \
                const MYFLT *S, \
                const MYFLT *pold, const MYFLT *vold, \
                MYFLT  *pnew, MYFLT *vnew) {
  //  calculate flow up to last grid point
  for (int m = 1; m <= M; m++) {  // Solve the first diff. eq. wrt vnew
    vnew[m] = -dt/rho_user * (pold[m] - pold[m-1]) / dx + vold[m];
  }
  //  calculate pressure up to one grid point before the last
  // (radiation in resonator cpp file)
  for (int m = 0; m <= M-1; m++) {  // Solve the second diff. eq. wrt pnew
    pnew[m] = - dt *rho_user *c*c / S[m] * (S[m+1]*vnew[m+1]
              - S[m]*vnew[m])/dx + pold[m];
  }
  vnew[0]  = 0;  // input flow is zero after the first sample
}


MYFLT cross_area(MYFLT radius, MYFLT x, MYFLT prelength, MYFLT slope) {
  /*x is the current position(m*dx), prelength is the length of the total tube
  up until the current segment (used for concatenation), radius is the
  input radius and slope is a ratio (dr/dx)*/
    MYFLT area;
    area = (radius+(x-prelength)*slope)*(radius+(x-prelength)*slope)*PI;
    return area;
}

void interpolation(int M, int Mold, MYFLT Lold, MYFLT dx, MYFLT dxold, \
                csnd::AuxMem<MYFLT> knew, csnd::AuxMem<MYFLT> kold){
  MYFLT x, xl, xr;
  int m1, m2;
  for (int m = 0; m<=std::min(M, int(floor(Lold/dx))); m++) {
    m1 = floor(m*(dx/dxold));
    m2 = m1 + 1;  // Equivalent to m2 = ceil(m*dx/dxold)

    xl = m1*dxold;
    xr = m2*dxold;
    x  = m*dx;
    kold[m] = knew[m1]*((xr - x)/(xr - xl)) + knew[m2]*((x-xl)/(xr-xl));
  }


  for (int m = std::min(M, int(floor(Lold/dx)))+1; m<= M; m++) {
    //    kold[m] = kold[std::min(M,int(floor(Lold/dx)))];
    kold[m] = knew[Mold];
  }
}

void interpolation_pointers(int M, int Mold, MYFLT Lold, MYFLT dx,
                        MYFLT dxold, MYFLT *knew, MYFLT *kold) {
  MYFLT x, xl, xr;
  int m1, m2;
  for (int m = 0; m<= std::min(M, int(floor(Lold/dx))); m++) {
    m1 = floor(m*(dx/dxold));
    m2 = m1 + 1;  // Equivalent to m2 = ceil(m*dx/dxold)

    xl = m1*dxold;
    xr = m2*dxold;
    x  = m*dx;
    kold[m] = knew[m1]*((xr - x)/(xr - xl)) + knew[m2]*((x-xl)/(xr-xl));
  }

  // when L increases k is set to last old value (e.g. pressure at tube end)
  for (int m = std::min(M, int(floor(Lold/dx)))+1; m<= M; m++ ) {
    //    kold[m] = kold[std::min(M,int(floor(Lold/dx)))];
    kold[m] = knew[Mold];
  }
}
