/*
  file: resonators/cone_radiation_losses_vp_fixedM_Bela.cpp
  opcode-name: halfphysler_bela

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

struct Cone_Radiation_Losses_fixedM_Bela : csnd::Plugin<2, 7> {
  /* Resonator with radiation losses, driven within
     an initial air velocity. The shape can be defined as a cone by passing the slope
     as an input argument. The number of grid points is set to fixed M=32 to run
     efficiently on Bela Mini. This restricts kLength to a range of 0.25 - 0.75 m.

     Here Latency is compensated for running with audioFrames=32 (32 sample delay)

     aFeedb, aSound cone_radiation_losses aPressure, kLength, kRad, kSlope, kEndReflection, kDensity
     AH (was cone_radiation_losses.cpp: )
  */

  csnd::AuxMem<MYFLT> pold;     // Pressure value at (n)th time grid
  csnd::AuxMem<MYFLT> vold;     // Velocity value at (n)th time grid
  csnd::AuxMem<MYFLT> pnew;     // Pressure value at (n+1)th time grid
  csnd::AuxMem<MYFLT> vnew;     // velocity value at (n+1)th time grid
  csnd::AuxMem<MYFLT> S;        // Cross sectional area

  // iterators to be passed to update_vp()
  csnd::AuxMem<MYFLT>::iterator iter_pold;
  csnd::AuxMem<MYFLT>::iterator iter_vold;
  csnd::AuxMem<MYFLT>::iterator iter_pnew;
  csnd::AuxMem<MYFLT>::iterator iter_vnew;
  csnd::AuxMem<MYFLT>::iterator iter_S;

  int M;                        // Number of steps
  int Mold;                     // Previous number of steps
  MYFLT L;                      // Length of the tube
  MYFLT Lold;                   // Previous length of the tube
  MYFLT dt;                     // Time grid
  MYFLT fs;                     // Sampling rate
  MYFLT dx;                     // Step size
  MYFLT dxold;                  // Previous step size
  MYFLT rad_alphaS;             // Radiation parameter alpha
  MYFLT rad_betaS;              // Radiation parameter beta
  MYFLT mult_alpha;             // user_multiplier for radiation losses
  MYFLT r;                      // Tube radius
  MYFLT r_old;                  // Prev. tube radius
  MYFLT slope;                  // Slope of the tube
  MYFLT slope_old;              // Prev. slope of the tube
  MYFLT mult_rho;               // user_multiplier for density
  MYFLT pickup_pos;             // user defined pickup position in tube

  int init() {
    L          = inargs[1];  // Length as input
    r          = inargs[2];  // Radius as input
    slope      = inargs[3];  // Slope as input
    mult_alpha = inargs[4];  // tube end reflection coefficient (radiation)
    mult_rho   = inargs[5];  // density coefficient
    pickup_pos = inargs[6];  // TODO(AH): relative pickup position M/2 to M.

    fs = csound->sr();
    dt = 1.0/fs;
    M = 32;  // max M for a length of 0.75 (max Bela-supported length)
    dx = L/M;

    // Memory Allocation
    pold.allocate(csound, M+1);
    pnew.allocate(csound, M+1);
    vold.allocate(csound, M+1);
    vnew.allocate(csound, M+1);
    S.allocate(csound, M+1);

    iter_pnew = pnew.begin();
    iter_vnew = vnew.begin();
    iter_pold = pold.begin();
    iter_vold = vold.begin();
    iter_S = S.begin();

    for (int m = 0; m<= M; m++) {    // Set all array elements to 0
      pold[m] = 0;
      vold[m] = 0;
      pnew[m] = 0;
      vnew[m] = 0;

      // Cross sectional area of tube  see ./src/tube.cpp
      S[m] = cross_area(r, m*dx, 0, slope);  // Input radius of cone
    }
    rad_alphaS = (rad_alpha * mult_alpha) / sqrt(S[M]);  // normalization
    rad_betaS = rad_beta / c;
    Lold      = L;
    r_old     = r;
    slope_old = slope;
    Mold      = M;
    dxold     = dx;

    return OK;
  }

  int aperf() {     // Calculate one audio block (ksmps samples)
    csnd::AudioSig out_feedback(this, outargs(0));   // feedback output (M = 0)
    csnd::AudioSig out_sound(this, outargs(1));   // sound out (at vari M != 0)
    csnd::AudioSig in(this, inargs(0));     // Csound opcode in
    L          = inargs[1];  // Length as input
    r          = inargs[2];  // Radius as input
    slope      = inargs[3];  // Slope as input
    mult_alpha = inargs[4];  // tube end reflection coefficient (radiation)
    mult_rho   = inargs[5];  // density coefficient
    pickup_pos = inargs[6];  //  TODO: relative pickup position M/2 to M

    rad_alphaS = (rad_alpha * mult_alpha) / sqrt(S[M]);  // normalization
    // re-calculate the grid
    if (L != Lold || r != r_old || slope != slope_old) {
    // Ensures that new calculations are made only when L,r or slope changed
      dx = L/M;

      // new cross-sectional area
      for (int m = 0; m<= M; m++) {
          S[m] = cross_area(r, m*dx, 0, slope);
      }

      // interpolate old grid status to new grid for each point
      interpolation(M, Mold, Lold, dx, dxold, iter_pnew, iter_pold);
      interpolation(M, Mold, Lold, dx, dxold, iter_vnew, iter_vold);
      rad_alphaS = (rad_alpha * mult_alpha) / sqrt(S[M]);  // normalization
  }  // Ending bracket of changed geometry

    int i = 0;
    for (auto &o_sound : out_sound) {  // For each sample ..
    /* ---  Coupling reed and resonator ---
       Output feedback read from an earlier position of the pressure wave to
       compensate for I/O latency.
       would be "out_feedback[i] = pnew[0];" if we had no latency at all */
     int del_smps = 32;
     int feedback_point = floor((del_smps * c) / (fs * dx));
     // printf("New Feedbp M equals to= %i\n", (feedback_point % M));
     out_feedback[i] = pnew[feedback_point % M];  // signal to shaker
     vnew[0] = in[i];   // reed signal (feedback)

     update_vp(M, dt, dx, c, (mult_rho * rho), iter_S, iter_pold, iter_vold,
                                                   iter_pnew, iter_vnew);

     // Boundary condition at tubeend has radiation losses, damps traveling wave
     pnew[M]  = (pold[M]*rad_betaS/(mult_rho * rho) + vnew[M]-vold[M]) /
                (rad_betaS/(mult_rho * rho) + rad_alphaS/(mult_rho * rho)*dt);

     i = i+1;
      // Copying p(n+1) to p(n) and v(n+1) to v(n),
      // i.e. new becomes old grid for next call
      std::copy(pnew.begin(), pnew.end(), pold.begin());
      std::copy(vnew.begin(), vnew.end(), vold.begin());

      /* ---  Sound output --- */
      //o_sound = pnew[M];  // Output the damped ending of the tube to speaker

      // sound is obtained at given position
      int pickup_idx = std::min(int(ceil(pickup_pos * L/dx)),M-1);
      o_sound = pnew[pickup_idx];  // Output the damped ending of the tube to csound
    }

    Lold      = L;
    r_old     = r;
    slope_old = slope;
    Mold      = M;
    dxold     = dx;

    return OK;
  }
};
