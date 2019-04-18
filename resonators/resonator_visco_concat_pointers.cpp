/*
  file: resonators/resonator_visco_concat_pointers.cpp
  opcode-name: tube_resonator opcode

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

#include <atomic>
#include <iostream>
#include <plugin.h>
#include <random>
#include "../src/tube.h"   // Includes the functions
#include "../src/const.h"  // Includes the constants

struct Resonator_Visco_Concat_Pointers : csnd::Plugin<1, 6> {
  // TODO: add infos here!
  csnd::AuxMem<MYFLT> pold; // Pressure value at (n)th time grid
  csnd::AuxMem<MYFLT> vold; // Velocity value at (n)th time grid
  csnd::AuxMem<MYFLT> pnew; // Pressure value at (n+1)th time grid
  csnd::AuxMem<MYFLT> vnew; // velocity value at (n+1)th time grid
  csnd::AuxMem<MYFLT> S;    //Cross sectional area

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
  MYFLT c_user;    // = 3.4386e+02;  // speed of sound
  MYFLT rho_user;  // = 1.2000e+00;  // density
  MYFLT eta_user;  // = 1.8200e-05;  // air viscosity
  MYFLT Zmult;  // specific impedance multiplier
  MYFLT Ymult;  // shunt admittance multiplier



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

  MYFLT x;            //Position for m*dx
  MYFLT xl;           //Position for m1*dxold
  MYFLT xr;           //Position for m2*dxold
  int m1;
  int m2;



  int init() {

    c_user            = 3.4386e+02;  //u-m speed of sound
    rho_user          = 1.2000e+00;  //u-m density
    eta_user          = 1.8200e-05;  //u-m air viscosity
    Zmult             = 1;
    Ymult             = 1;

    cone_lengths.allocate(csound,inargs.vector_data<MYFLT>(2).len());
    radii_in.allocate(csound,inargs.vector_data<MYFLT>(3).len());
    radii_out.allocate(csound,inargs.vector_data<MYFLT>(4).len());
    curve_type.allocate(csound,inargs.vector_data<MYFLT>(5).len());

    std::copy(inargs.vector_data<MYFLT>(2).begin(),inargs.vector_data<MYFLT>(2).end(),cone_lengths.begin());
    std::copy(inargs.vector_data<MYFLT>(3).begin(),inargs.vector_data<MYFLT>(3).end(),radii_in.begin());
    std::copy(inargs.vector_data<MYFLT>(4).begin(),inargs.vector_data<MYFLT>(4).end(),radii_out.begin());
    std::copy(inargs.vector_data<MYFLT>(5).begin(),inargs.vector_data<MYFLT>(5).end(),curve_type.begin());

    // --------------------------------------------------------------------------
    // Compute values for spatial and temporal grid
    Mmax = 400;
    fs = csound->sr();
    dt = 1.0/fs;
    grid_init_visco(.2, dt, Mmax, &dx, &M, &L, c_user); // setup the grid for finite differences



    // ---------------------------------------------------------------------------
    // Allocate state arrays
    pold.allocate(csound, Mmax+1); //Memory Allocation
    pnew.allocate(csound, Mmax+1);
    vold.allocate(csound, Mmax+1);
    vnew.allocate(csound, Mmax+1);
    iter_pold  = pold.begin();
    iter_pnew  = pnew.begin();
    iter_vold  = vold.begin();
    iter_vnew  = vnew.begin();

    // Cross-sectional area
    S.allocate(csound, Mmax+1);
    iter_S     = S.begin();

    // ---------------------------------------------------------------------------
    // Allocate loss arrays
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

    for (int m = 0; m<= M; m++){ //Set all array elements to 0
      pold[m] = 0;
      vold[m] = 0;
      pnew[m] = 0;
      vnew[m] = 0;
      S[m]    = cross_area_concatenation(cone_lengths, radii_in, radii_out, curve_type, m*dx, cone_lengths.len() ); //Specify cross sectional area
    }

    rad_alphaS = rad_alpha / sqrt(S[M]); // radiation parameters are normalized
    rad_betaS  = rad_beta/c_user;


    //    compute_loss_arrays(M, S, RsZ, LsZ, GsY, CsY, rz_tmp, lz_tmp, gy_tmp, \
    // cy_tmp, sumZ1, sumY1, dt, rho, c, eLZ, eCY, MATZ, MATSY, factors_v, factors_p);

    // ---------------------------------------------------------------------------
    // compute loss arrays
    compute_loss_arrays_pointers(M, iter_S, RsZ, LsZ, GsY, CsY, rz_tmp, lz_tmp, gy_tmp, \
                        cy_tmp, iter_sumZ1, iter_sumY1, dt, rho_user, c_user, iter_eLZ, iter_eCY, iter_MATZ,
                        iter_MATSY, iter_factors_v, iter_factors_p, Zmult, Ymult);


    // ---------------------------------------------------------------------------
    // initialize loss state arrays w and q to 0
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
    csnd::AudioSig out(this, outargs(0));
    csnd::AudioSig in(this, inargs(0)); // Csound opcode in
    //    L = inargs[1]; // Change in length according to input value from Csound
    std::copy(inargs.vector_data<MYFLT>(2).begin(),inargs.vector_data<MYFLT>(2).end(),cone_lengths.begin());
    std::copy(inargs.vector_data<MYFLT>(3).begin(),inargs.vector_data<MYFLT>(3).end(),radii_in.begin());
    std::copy(inargs.vector_data<MYFLT>(4).begin(),inargs.vector_data<MYFLT>(4).end(),radii_out.begin());
    std::copy(inargs.vector_data<MYFLT>(5).begin(),inargs.vector_data<MYFLT>(5).end(),curve_type.begin());

    L = inargs[1];

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
    for (auto &o : out) { // for each sample ..

      //      update_visco(M, sumZ1, sumY1, sumZ2, sumY2, eLZ, eCY, wlossold, qlossold, \
     // dx, dt, rho, factors_v, factors_p, S, vold, pold, vnew, pnew);

      update_visco_pointers(M, iter_sumZ1, iter_sumY1, iter_sumZ2, iter_sumY2, iter_eLZ, iter_eCY, iter_wlossold, iter_qlossold,
                dx, dt, rho_user, iter_factors_v, iter_factors_p, iter_S, iter_vold, iter_pold, iter_vnew, iter_pnew);

      vnew[0]  = in[i]; //Put comment here
      i++;
      pnew[M]  = (pold[M]*rad_betaS/rho + vnew[M]-vold[M]) / (rad_betaS/rho + rad_alphaS/rho*dt);
      o = pnew[M];



      //Updating arrays
//      for (int m = 0; m <= M; m++) {
//	for (int k = 0; k<4; k++){
//	  wloss[m][k]  =  wlossold[m][k]*eLZ[m][k] + (vnew[m]-vold[m])*MATZ[m][k];
//	  qloss[m][k]  =  qlossold[m][k]*eCY[m][k] + (pnew[m]-pold[m])*MATSY[m][k];
//	  wlossold[m][k] = wloss[m][k];
//	  qlossold[m][k] = qloss[m][k];
//	}
//	pold[m] = pnew[m];
//	vold[m] = vnew[m];
//      }


      // Updating arrays
      for (int m = 0; m <= M; m++) {
        for (int k = 0; k < 4; k++) {
            wloss[m*4+k]  =  wlossold[m*4+k]*eLZ[m*4+k]
                            + (vnew[m]-vold[m])*MATZ[m*4+k];
            qloss[m*4+k]  =  qlossold[m*4+k]*eCY[m*4+k]
                            + (pnew[m]-pold[m])*MATSY[m*4+k];
	    //            wlossold[m*4+k] = wloss[m*4+k];
	    //            qlossold[m*4+k] = qloss[m*4+k];
        }
	//        pold[m] = pnew[m];
	//        vold[m] = vnew[m];
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
