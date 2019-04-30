/*
  file: resonators/resonator_visco_concat.cpp
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

struct Resonator_Visco_Concat : csnd::Plugin<2, 10> {
    /* Resonator with radiation losses and viscothermal losses, driven with
       an initial air velocity. The shape can be defined by three arrays

       aFeedb, aSound tube_resonator aPressure, kLength, kRad, kSlope kRad, kSlope, kEndReflection, kDensity
    */

  // grid arrays for pressure and velocity
  csnd::AuxMem<MYFLT> pold;  // Pressure value at (n)th time grid
  csnd::AuxMem<MYFLT> vold;  // Velocity value at (n)th time grid
  csnd::AuxMem<MYFLT> pnew;  // Pressure value at (n+1)th time grid
  csnd::AuxMem<MYFLT> vnew;  // velocity value at (n+1)th time grid
  csnd::AuxMem<MYFLT> S;     // Cross sectional area

  // iterators to be passed to update_vp()
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
  MYFLT mult_alpha;  // user_multiplier for radiation losses

  // -------------------------------
  //MYFLT c_user;      //   = 3.4386e+02; // speed of sound
  //MYFLT eta_user;    // = 1.8200e-05; // air viscosity
  //MYFLT rho_user;    // = 1.2000e+00; // density

  MYFLT mult_rho;     // user_multiplier for density
  MYFLT Zmult;       // specific impedance multiplier
  MYFLT Ymult;       // shunt admittance multiplier
  MYFLT pickup_pos;  // user defined pickup position in tube

  // user gemoetry settings
  csnd::AuxMem<MYFLT> cone_lengths;
  csnd::AuxMem<MYFLT> radii_in;
  csnd::AuxMem<MYFLT> radii_out;
  csnd::AuxMem<MYFLT> curve_type;

  int maxGeoSegments;
  csnd::AuxMem<MYFLT> allGeoSettings_new;
  csnd::AuxMem<MYFLT> allGeoSettings_old;
  bool geometryChanged;
  int computeVisco;


  int init() {
    // --------------  user inputs --------------------------------
    L          = 0.1;  // inargs[1] Length as input
    pickup_pos = 0.5;  // inargs[2] TODO(AH): relative pickup position M/2 to M
    mult_alpha = 1.0;  // inargs[3] tube end reflection coefficient (radiation)
    mult_rho   = 1.0;  // inargs[4] density coefficient

    // it is currently (Apr. 2019) discussed in the Csnd-dev list,
    // if k-time array size changes are allowed or not
    // if not, we could set maxGeoSegments = (1 + 4*cone_lengths.len())
    maxGeoSegments = 25;      // Gemoetry can have up to 25 segments
    cone_lengths.allocate(csound, maxGeoSegments);
    radii_in.allocate(csound, maxGeoSegments);
    radii_out.allocate(csound, maxGeoSegments);
    curve_type.allocate(csound, maxGeoSegments);

    // inargs[5..8] copy user provided geometry to MYFLT arrays
    std::copy(inargs.vector_data<MYFLT>(5).begin(),
              inargs.vector_data<MYFLT>(5).end(), cone_lengths.begin());
    std::copy(inargs.vector_data<MYFLT>(6).begin(),
              inargs.vector_data<MYFLT>(6).end(), radii_in.begin());
    std::copy(inargs.vector_data<MYFLT>(7).begin(),
              inargs.vector_data<MYFLT>(7).end(), radii_out.begin());
    std::copy(inargs.vector_data<MYFLT>(8).begin(),
              inargs.vector_data<MYFLT>(8).end(), curve_type.begin());
    // TODO(AH): How to make this optional?
    // computeVisco = inargs[9];  // save CPU
     computeVisco = true;


    // copy all geometry into one long array, for comparison if changed
    allGeoSettings_new.allocate(csound, 1 + (maxGeoSegments*4));
    allGeoSettings_old.allocate(csound, 1 + (maxGeoSegments*4));
    geometryChanged = false;

    allGeoSettings_new[0] = L;
    std::copy(cone_lengths.begin(), cone_lengths.end(),
                allGeoSettings_new.begin()+1);
    std::copy(radii_in.begin(), radii_in.end(),
                allGeoSettings_new.begin() + 1 + maxGeoSegments);
    std::copy(radii_out.begin(), radii_out.end(),
                allGeoSettings_new.begin() + (1 + (maxGeoSegments*2)));
    std::copy(curve_type.begin(), curve_type.end(),
                allGeoSettings_new.begin() + (1 + (maxGeoSegments*3)));

    // TODO(Seb): What are these?? Do we want them as user inputs?
    Zmult             = 1;
    Ymult             = 1;

    // ------Compute values for spatial and temporal grid-----------------------
    fs = csound->sr();
    dt = 1.0/fs;
    grid_init(0.5, dt, &dx, &M, &L);  // setup grid for finite difference

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

    if (computeVisco)
        alloc_visco_memory(csound);

    // Init grid points with 0
    for (int m = 0; m<= M; m++) {
      pold[m] = 0;
      vold[m] = 0;
      pnew[m] = 0;
      vnew[m] = 0;
      MYFLT mdx = m*dx;
      //printf("point=%d, pos = %f\n", m, mdx);
      // Setting Cross sectional area for each grid point (see ./src/tube.cpp)
      S[m]    = cross_area_concatenation(cone_lengths, radii_in, radii_out, \
                                         curve_type, m*dx, cone_lengths.len() );
    }

    // normalize radiation parameters
    rad_alphaS = (rad_alpha * mult_alpha) / sqrt(S[M]);
    rad_betaS  = rad_beta / c;

    // --------Compute loss arrays---------------------------------------------
    if (computeVisco) {
        compute_loss_arrays(M, iter_S, RsZ, LsZ, GsY, CsY, dt, (mult_rho * rho), \
                            c, Zmult, Ymult);

        init_loss_state_array(M);
    }

    Lold  =  L;
    Mold  =  M;
    dxold  =  dx;
    std::copy(allGeoSettings_new.begin(),
                allGeoSettings_new.end(),
                allGeoSettings_old.begin());
    printf("init done\n");
    return OK;
  }



  int aperf() {  // calculate one audio block (at k-time)
    csnd::AudioSig out_feedback(this, outargs(0));   // feedback output (M = 0)
    csnd::AudioSig out_sound(this, outargs(1));   // sound out (at vari M != 0)
    csnd::AudioSig in(this, inargs(0));              // Csound opcode in

    // --------------  user inputs --------------------------------
    L          = inargs[1];  // Length as input
    pickup_pos = inargs[2];  // TODO(AH): relative pickup position M/2 to M
    mult_alpha = inargs[3];  // tube end reflection coefficient (radiation)
    mult_rho   = inargs[4];  // density coefficient

    cone_lengths.allocate(csound, maxGeoSegments);
    radii_in.allocate(csound, maxGeoSegments);
    radii_out.allocate(csound, maxGeoSegments);
    curve_type.allocate(csound, maxGeoSegments);

    // copy user provided geometry to MYFLT arrays
    std::copy(inargs.vector_data<MYFLT>(5).begin(),
              inargs.vector_data<MYFLT>(5).end(), cone_lengths.begin());
    std::copy(inargs.vector_data<MYFLT>(6).begin(),
              inargs.vector_data<MYFLT>(6).end(), radii_in.begin());
    std::copy(inargs.vector_data<MYFLT>(7).begin(),
              inargs.vector_data<MYFLT>(7).end(), radii_out.begin());
    std::copy(inargs.vector_data<MYFLT>(8).begin(),
              inargs.vector_data<MYFLT>(8).end(), curve_type.begin());
    // TODO(AH): How to make this optional?
    computeVisco = inargs[9];  // save CPU



    // ------------ check if geometry has changed in one array -----------------
    allGeoSettings_new[0] = L;
    // this is a bit redundant, but ...nicer for the comparison
    std::copy(cone_lengths.begin(), cone_lengths.end(),
                allGeoSettings_new.begin()+1);
    std::copy(radii_in.begin(), radii_in.end(),
                allGeoSettings_new.begin() + 1 + maxGeoSegments);
    std::copy(radii_out.begin(), radii_out.end(),
                allGeoSettings_new.begin() + (1 + (maxGeoSegments*2)));
    std::copy(curve_type.begin(), curve_type.end(),
                allGeoSettings_new.begin() + (1 + (maxGeoSegments*3)));

    // compare until the first difference, in best case its kLength at elem 0
    for (int x = 0; x <= (1 + (maxGeoSegments*4)); x++) {
        if (allGeoSettings_new[x] != allGeoSettings_old[x]) {
            geometryChanged = true;
            break;
        }
    }

    // ------------ Re-calculate the grid ---------------------------
    if (geometryChanged) {  // new geometry calculations only when length change
        grid_init(L, dt, &dx, &M, &L);  // setup grid for finite difference
        for (int m = 0; m<= M; m++) {  // new cross sectional area
          S[m] = cross_area_concatenation(cone_lengths, radii_in, radii_out,
                                       curve_type, m*dx, cone_lengths.len());
        }
        rad_alphaS = (rad_alpha * mult_alpha) / sqrt(S[M]);  // normalization

        // -------- interpolate old grid status to new grid for each point------
        interpolation(M, Mold, Lold, dx, dxold, iter_pnew, iter_pold);
        interpolation(M, Mold, Lold, dx, dxold, iter_vnew, iter_vold);
        if (computeVisco) {
            interpolation_visco_arrays(M, Mold, Lold, dx, dxold);
            compute_loss_arrays(M, iter_S, RsZ, LsZ, GsY, CsY, dt, (mult_rho * rho), \
                                c, Zmult, Ymult);
            }
    }  // Ending bracket for geometryChanged

    int i  = 0;
    for (auto &o_sound : out_sound) {  // For each sample ..
        out_feedback[i] = pnew[0];  // Output pressure, tube begin for coupling
        vnew[0]  = in[i];  // Input external velocity at beginning of tube
        if (computeVisco) {
            update_visco(M, dx, dt, (mult_rho * rho), iter_S,
                        iter_vold, iter_pold, iter_vnew, iter_pnew);
        } else {
            int mult_rho = 1;
            update_vp(M, dt, dx, c, (mult_rho * rho), iter_S, iter_pold,
                    iter_vold, iter_pnew, iter_vnew);
            }
        // Boundary condition at tube end has radiation losses, damps traveling wave
        pnew[M]  = (pold[M]*rad_betaS/rho + vnew[M]-vold[M]) / (rad_betaS/rho + rad_alphaS/rho*dt);
        i++;

        // sound output at variable grid point
        int pickup_idx = std::min(int(ceil(pickup_pos * L/dx)), M-1);
        o_sound = pnew[pickup_idx];

        if (computeVisco) {
            update_losses(M, iter_vold, iter_pold, iter_vnew, iter_pnew);
        }
        // Copying p(n+1) to p(n) and v(n+1) to v(n),
        // i.e. new becomes old grid for next call
        std::copy(pnew.begin(), pnew.end(), pold.begin());
        std::copy(vnew.begin(), vnew.end(), vold.begin());
    }

    Lold            = L;
    Mold            = M;
    dxold           = dx;
    std::copy(allGeoSettings_new.begin(), allGeoSettings_new.end(),
              allGeoSettings_old.begin());
    geometryChanged = false;

    return OK;
  }
};
