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
  int Mmax;          // Maximal number of steps
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

  // viscothermal loss variables
  MYFLT rz_tmp[4];
  MYFLT lz_tmp[4];
  MYFLT gy_tmp[4];
  MYFLT cy_tmp[4];

  csnd::AuxMem<MYFLT> eLZ;
  csnd::AuxMem<MYFLT> eCY;
  csnd::AuxMem<MYFLT> MATZ;
  csnd::AuxMem<MYFLT> MATSY;
  csnd::AuxMem<MYFLT> qloss;
  csnd::AuxMem<MYFLT> wloss;
  csnd::AuxMem<MYFLT> qlossold;
  csnd::AuxMem<MYFLT> wlossold;

  // iterators
  csnd::AuxMem<MYFLT>::iterator iter_eLZ;
  csnd::AuxMem<MYFLT>::iterator iter_eCY;
  csnd::AuxMem<MYFLT>::iterator iter_MATZ;
  csnd::AuxMem<MYFLT>::iterator iter_MATSY;
  csnd::AuxMem<MYFLT>::iterator iter_qloss;
  csnd::AuxMem<MYFLT>::iterator iter_wloss;
  csnd::AuxMem<MYFLT>::iterator iter_qlossold;
  csnd::AuxMem<MYFLT>::iterator iter_wlossold;

  csnd::AuxMem<MYFLT> sumZ1;
  csnd::AuxMem<MYFLT> factors_v;
  csnd::AuxMem<MYFLT> sumY1;
  csnd::AuxMem<MYFLT> factors_p;

  csnd::AuxMem<MYFLT> sumZ2;
  csnd::AuxMem<MYFLT> sumY2;

  // iterators
  csnd::AuxMem<MYFLT>::iterator iter_sumZ1;
  csnd::AuxMem<MYFLT>::iterator iter_factors_v;
  csnd::AuxMem<MYFLT>::iterator iter_sumY1;
  csnd::AuxMem<MYFLT>::iterator iter_factors_p;
  csnd::AuxMem<MYFLT>::iterator iter_sumZ2;
  csnd::AuxMem<MYFLT>::iterator iter_sumY2;

  // user gemoetry settings
  csnd::AuxMem<MYFLT> cone_lengths;
  csnd::AuxMem<MYFLT> radii_in;
  csnd::AuxMem<MYFLT> radii_out;
  csnd::AuxMem<MYFLT> curve_type;

  MYFLT x;            //  Position for m*dx
  MYFLT xl;           //  Position for m1*dxold
  MYFLT xr;           //  Position for m2*dxold
  int m1;
  int m2;


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
    Mmax = 400;
    fs = csound->sr();
    dt = 1.0/fs;
    grid_init_visco(.2, dt, Mmax, &dx, &M, &L, c_user); // setup the grid for finite differences

    // -----Allocate memory for grid state arrays------------------------------
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

    // -----Memory allocation for viscothermal loss calculations---------------
    // Allocate loss arrays
    // Note: To avoid matrices in aperf, eLZ[m][k]
    // for k=4 (precision of Sebastian??)
    // Arrays hold concatenated viscothermal loss factors instead.
    eLZ.allocate(csound, 4*(Mmax+1));
    iter_eLZ  =  eLZ.begin();
    eCY.allocate(csound, 4*(Mmax+1));
    iter_eCY  =  eCY.begin();
    MATZ.allocate(csound, 4*(Mmax+1));
    iter_MATZ  =  MATZ.begin();
    MATSY.allocate(csound, 4*(Mmax+1));
    iter_MATSY  =  MATSY.begin();
    qloss.allocate(csound, 4*(Mmax+1));
    iter_qloss  =  qloss.begin();
    wloss.allocate(csound, 4*(Mmax+1));
    iter_wloss  =  wloss.begin();
    qlossold.allocate(csound, 4*(Mmax+1));
    iter_qlossold  =  qlossold.begin();
    wlossold.allocate(csound, 4*(Mmax+1));
    iter_wlossold  =  wlossold.begin();

    sumZ1.allocate(csound, Mmax+1);
    iter_sumZ1  =  sumZ1.begin();
    factors_v.allocate(csound, Mmax+1);
    iter_factors_v  =  factors_v.begin();
    sumY1.allocate(csound, Mmax+1);
    iter_sumY1  =  sumY1.begin();
    factors_p.allocate(csound, Mmax+1);;
    iter_factors_p  =  factors_p.begin();
    sumZ2.allocate(csound, Mmax+1);;
    iter_sumZ2  =  sumZ2.begin();
    sumY2.allocate(csound, Mmax+1);;
    iter_sumY2  =  sumY2.begin();

    // Init grid
    for (int m = 0; m<= M; m++){
      pold[m] = 0;
      vold[m] = 0;
      pnew[m] = 0;
      vnew[m] = 0;
      // Setting Cross sectional area of tube (see ./src/tube.cpp)
      S[m]    = cross_area_concatenation(cone_lengths, radii_in, radii_out, \
                                         curve_type, m*dx, cone_lengths.len() );
    }

    rad_alphaS = rad_alpha / sqrt(S[M]); // normalize radiation parameters
    rad_betaS  = rad_beta / c_user;

    // --------Compute loss arrays---------------------------------------------
    compute_loss_arrays_pointers(M, iter_S, RsZ, LsZ, GsY, CsY, rz_tmp, lz_tmp,\
                        gy_tmp, cy_tmp, iter_sumZ1, iter_sumY1, dt, rho_user, \
                        c_user, iter_eLZ, iter_eCY, iter_MATZ, iter_MATSY, \
                        iter_factors_v, iter_factors_p, Zmult, Ymult);

    // -------initialize loss state arrays w and q = 0-----------------------
    for (int m = 0; m <= M; m++) {
      for (int k = 0; k < 4; k++) {
        wloss[m*4+k] = 0;
        qloss[m*4+k] = 0;
        wlossold[m*4+k] = 0;
        qlossold[m*4+k] = 0;
      }
    }
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


    if(inargs[1]!=Lold)
      { //Ensures that new calculations are made only when L is changed

        fs = csound->sr();
        dt = 1.0/fs;
        grid_init_visco(inargs[1], dt, Mmax, &dx, &M, &L, c_user);  // reset grid due to length changes

	for (int m = 0; m<= M; m++){
	  S[m] = cross_area_concatenation(cone_lengths, radii_in, radii_out, curve_type, m*dx, cone_lengths.len() );
  }

	rad_alphaS = rad_alpha / sqrt(S[M]);

	//  interpolation(M, Mold, Lold, dx, dxold, pnew, pold);
	//  interpolation(M, Mold, Lold, dx, dxold, vnew, vold);
	//  interpolation_visco(M, Mold, Lold, dx, dxold, wloss, wlossold);
	//  interpolation_visco(M, Mold, Lold, dx, dxold, qloss, qlossold);

	//  compute_loss_arrays(M, S, RsZ, LsZ, GsY, CsY, rz_tmp, lz_tmp, gy_tmp, \
  cy_tmp, sumZ1, sumY1, dt, rho, c, eLZ, eCY, MATZ, MATSY, factors_v, factors_p);

// Why are both called here??
    interpolation_pointers(M, Mold, Lold, dx, dxold, iter_pnew, iter_pold);
    interpolation_pointers(M, Mold, Lold, dx, dxold, iter_vnew, iter_vold);
    interpolation_visco_pointers(M, Mold, Lold, dx, dxold, iter_wloss, iter_wlossold);
    interpolation_visco_pointers(M, Mold, Lold, dx, dxold, iter_qloss, iter_qlossold);

    compute_loss_arrays_pointers(M, iter_S, RsZ, LsZ, GsY, CsY, rz_tmp, lz_tmp, gy_tmp, \
                        cy_tmp, iter_sumZ1, iter_sumY1, dt, rho_user, c_user, iter_eLZ, iter_eCY,
                        iter_MATZ, iter_MATSY, iter_factors_v, iter_factors_p, eta_user, Zmult, Ymult);


      } //Ending bracket of if(L!=Lold)



    int i  = 0;
    for (auto &o_sound : out_sound) { // for each sample ..

      //      update_visco(M, sumZ1, sumY1, sumZ2, sumY2, eLZ, eCY, wlossold, qlossold, \
      // dx, dt, rho, factors_v, factors_p, S, vold, pold, vnew, pnew);

        out_feedback[i] = pnew[0];
        vnew[0]  = in[i]; //Put comment here
        update_visco_pointers(M, iter_sumZ1, iter_sumY1, iter_sumZ2, iter_sumY2, iter_eLZ, iter_eCY, iter_wlossold, iter_qlossold,
                dx, dt, rho_user, iter_factors_v, iter_factors_p, iter_S, iter_vold, iter_pold, iter_vnew, iter_pnew);

        pnew[M]  = (pold[M]*rad_betaS/rho + vnew[M]-vold[M]) / (rad_betaS/rho + rad_alphaS/rho*dt);
        i++;
        int pickup_idx = std::min(int(ceil(pickup_pos * L/dx)),M-1);
        //printf("%f\n", pickup_pos);
        o_sound = pnew[pickup_idx];


      // Updating arrays
      for (int m = 0; m <= M; m++) {
        for (int k = 0; k < 4; k++) {
            wloss[m*4+k]  =  wlossold[m*4+k]*eLZ[m*4+k]
                            + (vnew[m]-vold[m])*MATZ[m*4+k];
            qloss[m*4+k]  =  qlossold[m*4+k]*eCY[m*4+k]
                            + (pnew[m]-pold[m])*MATSY[m*4+k];
        }
      }
      std::copy(pnew.begin(), pnew.end(), pold.begin());
      std::copy(vnew.begin(), vnew.end(), vold.begin());
      std::copy(wloss.begin(), wloss.end(), wlossold.begin());
      std::copy(qloss.begin(), qloss.end(), qlossold.begin());



    }
    Lold = L;
    Mold = M;
    dxold = dx;

    return OK;
  }
};
