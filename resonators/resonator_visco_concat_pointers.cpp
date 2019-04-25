/*
  file: resonators/resonator_visco_concat_pointers.cpp
  opcode-name: tube_resonator

  Copyright (C) 2018 - Alex Hofmann, Vasileios Chatziioannou,
                       Sebastian Schmutzhard, Gokberk Erdogan

  C++ Implementation of Wind instrument models within the
  Csound-Plug-in Framwork by Lazzarini (SMC2017 - Paper)

  Resonator models taken from:

  **Schmutzhard, Sebastian; Chatziioannou, Vasileios, and Hofmann, Alex (2017)
    "Parameter Optimisation of a Viscothermal Time-Domain Model for Wind
    Instruments," in Proceedings of the 2017 International Symposium on
    Musical Acoustics (Montreal, CA) p. 27--30.(ISMA2017)**

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
#include <atomic>
#include <iostream>
#include <random>
#include "../src/tube.h"   // Includes the functions
#include "../src/const.h"  // Includes the constants

struct Resonator_Visco_Concat_Pointers : csnd::Plugin<2, 7> {
    /* Resonator with radiation losses and viscothermal losses, driven with
       an initial air velocity. The shape can be defined by three arrays

       aFeedb, aSound tube_resonator aPressure, kLength, kRad, kSlope kRad, kSlope, kEndReflection, kDensity
    */

  csnd::AuxMem<MYFLT> pold;  // Pressure value at (n)th time grid
  csnd::AuxMem<MYFLT> vold;  // Velocity value at (n)th time grid
  csnd::AuxMem<MYFLT> pnew;  // Pressure value at (n+1)th time grid
  csnd::AuxMem<MYFLT> vnew;  // velocity value at (n+1)th time grid
  csnd::AuxMem<MYFLT> S;     // Cross sectional area

  // iterators to be passed to update_vp_pointers()
  csnd::AuxMem<MYFLT>::iterator iter_pold;
  csnd::AuxMem<MYFLT>::iterator iter_vold;
  csnd::AuxMem<MYFLT>::iterator iter_pnew;
  csnd::AuxMem<MYFLT>::iterator iter_vnew;
  csnd::AuxMem<MYFLT>::iterator iter_S;


  int M;             // Number of steps
  int Mold;          // Previous number of steps
  MYFLT L;           // Length of the tube
  MYFLT Lold;        // Previous length of the tube
  MYFLT dt;          // Time grid
  MYFLT fs;          // Sampling rate
  MYFLT dx;          // Step size
  MYFLT dxold;       // Previous step size
  MYFLT rad_alphaS;  // Radiation parameter alpha
  MYFLT rad_betaS;   // Radiation parameter beta

  // -------------------------------
  MYFLT c_user;      //   = 3.4386e+02; // speed of sound
  MYFLT rho_user;    // = 1.2000e+00; // density
  MYFLT eta_user;    // = 1.8200e-05; // air viscosity
  MYFLT Zmult;       // specific impedance multiplier
  MYFLT Ymult;       // shunt admittance multiplier
  MYFLT pickup_pos;  // user defined pickup position in tube

  // user gemoetry settings
  csnd::AuxMem<MYFLT> cone_lengths;
  csnd::AuxMem<MYFLT> radii_in;
  csnd::AuxMem<MYFLT> radii_out;
  csnd::AuxMem<MYFLT> curve_type;

  // viscothermal loss variables

  MYFLT rz_tmp[4];
  MYFLT lz_tmp[4];
  MYFLT gy_tmp[4];
  MYFLT cy_tmp[4];


  int init() {
    L          = inargs[1];  // Length as input
    pickup_pos = inargs[2];  // TODO(AH): relative pickup position M/2 to M
    cone_lengths.allocate(csound, inargs.vector_data<MYFLT>(3).len());
    radii_in.allocate(csound, inargs.vector_data<MYFLT>(4).len());
    radii_out.allocate(csound, inargs.vector_data<MYFLT>(5).len());
    curve_type.allocate(csound, inargs.vector_data<MYFLT>(6).len());

    // copy user provided geometry to MYFLT arrays
    std::copy(inargs.vector_data<MYFLT>(3).begin(), inargs.vector_data<MYFLT>(3).end(), cone_lengths.begin());
    std::copy(inargs.vector_data<MYFLT>(4).begin(), inargs.vector_data<MYFLT>(4).end(), radii_in.begin());
    std::copy(inargs.vector_data<MYFLT>(5).begin(), inargs.vector_data<MYFLT>(5).end(), radii_out.begin());
    std::copy(inargs.vector_data<MYFLT>(6).begin(), inargs.vector_data<MYFLT>(6).end(), curve_type.begin());

    c_user            = 3.4386e+02;  // u-m speed of sound
    rho_user          = 1.2000e+00;  // u-m density
    eta_user          = 1.8200e-05;  // u-m air viscosity
    Zmult             = 1;
    Ymult             = 1;

    // ------Compute values for spatial and temporal grid-----------------------
    fs = csound->sr();
    dt = 1.0/fs;
    grid_init(L, dt, &dx, &M, &L);  // setup grid for finite difference

    // -----Allocate memory for grid state arrays------------------------------
    // using Mmax instead of M, allows 'Length' changes w/o new alloc in aperf()
    pold.allocate(csound, Mmax+1);
    pnew.allocate(csound, Mmax+1);
    vold.allocate(csound, Mmax+1);
    vnew.allocate(csound, Mmax+1);
    S.allocate(csound, Mmax+1);  // Cross-sectional area

    iter_pold  = pold.begin();
    iter_pnew  = pnew.begin();
    iter_vold  = vold.begin();
    iter_vnew  = vnew.begin();
    iter_S     = S.begin();

    alloc_visco_memory(csound);

    // Init grid points with 0
    for (int m = 0; m<= M; m++) {
      pold[m] = 0;
      vold[m] = 0;
      pnew[m] = 0;
      vnew[m] = 0;
      MYFLT mdx = m*dx;
      printf("point=%d, pos = %f\n", m, mdx);
      // Setting Cross sectional area for each grid point (see ./src/tube.cpp)
      S[m]    = cross_area_concatenation(cone_lengths, radii_in, radii_out, \
                                         curve_type, m*dx, cone_lengths.len() );
    }

    rad_alphaS = rad_alpha / sqrt(S[M]); // normalize radiation parameters
    rad_betaS  = rad_beta / c_user;

    // --------Compute loss arrays---------------------------------------------

    compute_loss_arrays_pointers(M, iter_S, RsZ, LsZ, GsY, CsY, rz_tmp, lz_tmp,\
                        gy_tmp, cy_tmp, dt, rho_user, \
                        c_user, Zmult, Ymult);

    init_loss_state_array(M);
    Lold  =  L;
    Mold  =  M;
    dxold  =  dx;

    return OK;
  }

  int aperf() { // calculate one audio block
    csnd::AudioSig out_feedback(this, outargs(0));   // feedback output (M = 0)
    csnd::AudioSig out_sound(this, outargs(1));   // sound out (at vari M != 0)
    csnd::AudioSig in(this, inargs(0));              // Csound opcode in
    L = inargs[1];   // Length as input
    pickup_pos = inargs[2];  // TODO(AH): relative pickup position M/2 to M
    std::copy(inargs.vector_data<MYFLT>(3).begin(),inargs.vector_data<MYFLT>(3).end(),cone_lengths.begin());
    std::copy(inargs.vector_data<MYFLT>(4).begin(),inargs.vector_data<MYFLT>(4).end(),radii_in.begin());
    std::copy(inargs.vector_data<MYFLT>(5).begin(),inargs.vector_data<MYFLT>(5).end(),radii_out.begin());
    std::copy(inargs.vector_data<MYFLT>(6).begin(),inargs.vector_data<MYFLT>(6).end(),curve_type.begin());

    // ------------ Re-calculate the grid ---------------------------
    if(inargs[1]!=Lold) {  // new geometry calculations only when length changed
        grid_init(L, dt, &dx, &M, &L);  // setup grid for finite difference
        for (int m = 0; m<= M; m++){  // new cross sectional area
          S[m] = cross_area_concatenation(cone_lengths, radii_in, radii_out,
                                       curve_type, m*dx, cone_lengths.len());
        }

        rad_alphaS = rad_alpha / sqrt(S[M]);  // normalize radiation parameters

        // -------- interpolate old grid status to new grid for each point--------
        interpolation_pointers(M, Mold, Lold, dx, dxold, iter_pnew, iter_pold);
        interpolation_pointers(M, Mold, Lold, dx, dxold, iter_vnew, iter_vold);

        // AH Problem, as same function called with diff varis.. qloss and wloss
        //AHref interpolation_visco_pointers(M, Mold, Lold, dx, dxold, iter_wloss, iter_wlossold);
        //AHref interpolation_visco_pointers(M, Mold, Lold, dx, dxold, iter_qloss, iter_qlossold);
        compute_loss_arrays_pointers(M, iter_S, RsZ, LsZ, GsY, CsY, rz_tmp, lz_tmp,\
                            gy_tmp, cy_tmp, dt, rho_user, \
                            c_user, Zmult, Ymult);
    } //Ending bracket for changed length

    int i  = 0;
    for (auto &o_sound : out_sound) {  // For each sample ..
        out_feedback[i] = pnew[0];  // Output pressure at tube begin for coupling
        vnew[0]  = in[i];  // Input external velocity at beginning of tube
        update_visco_pointers(M,
                dx, dt, rho_user, iter_S, iter_vold, iter_pold, iter_vnew, iter_pnew);

        // Boundary condition at tube end has radiation losses, damps traveling wave
        pnew[M]  = (pold[M]*rad_betaS/rho + vnew[M]-vold[M]) / (rad_betaS/rho + rad_alphaS/rho*dt);
        i++;

        // sound output at variable grid point
        int pickup_idx = std::min(int(ceil(pickup_pos * L/dx)),M-1);
        o_sound = pnew[pickup_idx];

        update_losses(M, iter_vold, iter_pold, iter_vnew, iter_pnew);
        // Copying p(n+1) to p(n) and v(n+1) to v(n),
        // i.e. new becomes old grid for next call
        std::copy(pnew.begin(), pnew.end(), pold.begin());
        std::copy(vnew.begin(), vnew.end(), vold.begin());
        //AHref std::copy(wloss.begin(), wloss.end(), wlossold.begin());
        //AHref std::copy(qloss.begin(), qloss.end(), qlossold.begin());

    }
    Lold = L;
    Mold = M;
    dxold = dx;

    return OK;
  }
};
