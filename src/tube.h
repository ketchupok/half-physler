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

//void grid_init_visco(MYFLT Len, MYFLT dt, int Mmax, MYFLT *dx, int *M, MYFLT *L, MYFLT c = 3.4386e+02);
void alloc_visco_memory(csnd::Csound *csound);

void init_loss_state_array(int M);

void update_vp_pointers(int M, const MYFLT& dt, const MYFLT& dx, const MYFLT& c,
                        const MYFLT& rho_user, const MYFLT *S,
                        const MYFLT *pold, const MYFLT *vold,
                        MYFLT *pnew, MYFLT *vnew);

void update_visco_pointers(int M, MYFLT dx, \
			 MYFLT dt, MYFLT rho, \
			 MYFLT* S, MYFLT* vold, MYFLT* pold, \
			 MYFLT* vnew, MYFLT* pnew);

void update_losses(int M, MYFLT* vold, MYFLT* pold, MYFLT* vnew, MYFLT* pnew);

MYFLT cross_area(MYFLT radius, MYFLT x, MYFLT prelength, MYFLT slope);

MYFLT cross_area_concatenation(csnd::AuxMem<MYFLT> cone_lengths, csnd::AuxMem<MYFLT> radii_in, csnd::AuxMem<MYFLT> radii_out, csnd::AuxMem<MYFLT> curve_type, MYFLT x, int size);

MYFLT linear_approx(MYFLT radius_in, MYFLT radius_out, MYFLT x, MYFLT prelength, MYFLT length);

MYFLT parabolic_approx(MYFLT radius_in, MYFLT radius_out, MYFLT x, MYFLT prelength, MYFLT length);

MYFLT exponential_approx(MYFLT radius_in, MYFLT radius_out, MYFLT x, MYFLT prelength, MYFLT length);

void interp_loss(MYFLT rad, MYFLT coeff[4][100], MYFLT * r); // TODO: change to pointer

void interpolation_pointers(int M, int Mold, MYFLT Lold, MYFLT dx, MYFLT dxold,
                            MYFLT *knew, MYFLT *kold);

void interpolation_visco_pointers(int M, int Mold, MYFLT Lold, MYFLT dx, MYFLT dxold, \
            MYFLT* klossnew, MYFLT*  klossold);

void interpolation_visco_arrays(int M, int Mold, MYFLT Lold, MYFLT dx, MYFLT dxold);

MYFLT R0Z(MYFLT r, MYFLT rho, MYFLT eta);

void compute_loss_arrays_pointers(int M, MYFLT* S, MYFLT RsZ[4][100], \
			 MYFLT LsZ[4][100], MYFLT GsY[4][100], MYFLT CsY[4][100],
             MYFLT dt, MYFLT rho, MYFLT c, \
       MYFLT Zmult = 1, MYFLT Ymult = 1);



#endif  // SRC_TUBE_H_
