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

  MYFLT rz_tmp[4];
  MYFLT lz_tmp[4];
  MYFLT gy_tmp[4];
  MYFLT cy_tmp[4];

  // iterators
  csnd::AuxMem<MYFLT>::iterator iter_sumZ1;
  csnd::AuxMem<MYFLT>::iterator iter_factors_v;
  csnd::AuxMem<MYFLT>::iterator iter_sumY1;
  csnd::AuxMem<MYFLT>::iterator iter_factors_p;
  csnd::AuxMem<MYFLT>::iterator iter_sumZ2;
  csnd::AuxMem<MYFLT>::iterator iter_sumY2;


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

void alloc_visco_memory(csnd::Csound* csound) {
    // -----Memory allocation for viscothermal loss calculations---------------
    // Allocate loss arrays
    // Note: To avoid matrices in aperf, eLZ[m][k]
    // for k=4 (Bilbao Harrison, 2016)
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
}

void init_loss_state_array(int M) {
    // -------initialize loss state arrays w and q = 0-----------------------
// was after compute_loss_arrays()
    for (int m = 0; m <= M; m++) {
      for (int k = 0; k < 4; k++) {
        wloss[m*4+k] = 0;
        qloss[m*4+k] = 0;
        wlossold[m*4+k] = 0;
        qlossold[m*4+k] = 0;
      }
    }
}

void update_vp(int M, const MYFLT& dt, const MYFLT& dx, \
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
  //vnew[0]  = 0;  // input flow is zero after the first sample
}

void update_visco(int M, MYFLT dx, \
                        MYFLT dt, MYFLT rho_user, \
                        MYFLT* S, MYFLT* vold, MYFLT* pold, \
                        MYFLT* vnew, MYFLT* pnew) {
      for (int m = 0; m <=  M; m++) {  // Solve the first diff. eq. wrt vnew
        sumZ2[m]  =  0;
        sumY2[m]  =  0;
        for (int k = 0; k < 4; k++) {
          sumZ2[m] += eLZ[m*4+k]*wlossold[m*4+k];
          sumY2[m] += eCY[m*4+k]*qlossold[m*4+k];
        }
      }

      for (int m = 1; m <=  M; m++) {  // Solve the first diff. eq. wrt vnew
        vnew[m] = vold[m]*(rho_user/dt + sumZ1[m]) - 1/dx *(pold[m]-pold[m-1]) - sumZ2[m];
        vnew[m] = vnew[m]/factors_v[m];
      }

      for (int m = 0; m <= M-1; m++) {  // Solve the second diff. eq. wrt pnew
        pnew[m] = pold[m] - (1/dx * (S[m+1]*vnew[m+1] - S[m]*vnew[m]) + sumY2[m])/factors_p[m];
      }
      // TODO(Vasi/Seb) why is not vnew[0] = 0 in this case??
      // vnew[0]  = 0;  // input flow is zero after the first sample
}

void update_losses(int M, MYFLT* vold, MYFLT* pold, MYFLT* vnew, MYFLT* pnew) {
        // Updating losses at each grid point
        for (int m = 0; m <= M; m++) {
            for (int k = 0; k < 4; k++) {
                wloss[m*4+k]  =  wlossold[m*4+k]*eLZ[m*4+k]
                                + (vnew[m]-vold[m])*MATZ[m*4+k];
                qloss[m*4+k]  =  qlossold[m*4+k]*eCY[m*4+k]
                                + (pnew[m]-pold[m])*MATSY[m*4+k];
            }
        }
        std::copy(wloss.begin(), wloss.end(), wlossold.begin());
        std::copy(qloss.begin(), qloss.end(), qlossold.begin());
}

MYFLT cross_area(MYFLT radius, MYFLT x, MYFLT prelength, MYFLT slope) {
  /*x is the current position(m*dx), prelength is the length of the total tube
  up until the current segment (used for concatenation), radius is the
  input radius and slope is a ratio (dr/dx)*/
    MYFLT area;
    area = (radius+(x-prelength)*slope)*(radius+(x-prelength)*slope)*PI;
    return area;
}

MYFLT cross_area_concatenation(csnd::AuxMem<MYFLT> cone_lengths, \
                                csnd::AuxMem<MYFLT> radii_in,    \
                                csnd::AuxMem<MYFLT> radii_out,   \
                                csnd::AuxMem<MYFLT> curve_type,  \
                                MYFLT grid_pos, int num_segs) {
    /* complex geometries are given as segments by the user
        cross sectional areas for each are calculated according to these geometries
        jumps between segment radii can be interpolated using three different algorithms
        a) linear, b) exponential, c) parabolic
    */

  MYFLT prelength = 0.0;  // tube length until the current segment
  int seg = 0;
  MYFLT area;

  // determine in which segment we are in, set prelength to end of last segment
  for (int i = 0; i < num_segs; i++) {
    prelength += cone_lengths[i];  // sum up total length to determine m and dx
    if ((grid_pos <=  prelength) || (i == num_segs-1) ){ // due to round off error, x might exceed total length in the last segment. || (i == size-1) fixes that.
      seg = i;      // Location of current m*dx, 1st cone, 2nd segment.
      prelength -= cone_lengths[i];
      break;
    }
  }

//printf("prelength before switch: %f\n", prelength);
  switch (int(curve_type[seg])) {  // converting MYFLT to int for switch()
    case 1:  // 1 for linear approximation
      area = linear_approx(radii_in[seg], radii_out[seg], grid_pos,  \
                            prelength, cone_lengths[seg]);
      break;
    case 2:  // 2 for parabolic approximation
      area = parabolic_approx(radii_in[seg], radii_out[seg], grid_pos,  \
                                prelength, cone_lengths[seg]);
      break;
    case 3:  // 3 for exponential approximation
      area = exponential_approx(radii_in[seg], radii_out[seg], grid_pos,  \
                                prelength, cone_lengths[seg]);
      break;
  }
  //printf("area: %f\n", area);
  return area;
}

MYFLT linear_approx(MYFLT radius_in, MYFLT radius_out, MYFLT grid_pos,  \
                    MYFLT prelength, MYFLT seg_length) {
  /*grid_pos is the current position(m*dx), prelength is the length of the total tube
  up until the current segment, seg_length is the length of the current segment*/
    MYFLT area;
    MYFLT r;
    if (seg_length == 0) {  // TODO(AH): how to deal with negative numbers?? <=
        r = radius_in;
    } else {
        MYFLT distance_to_seg_begin = grid_pos - prelength;
        MYFLT radius_diff_per_meter = (radius_out - radius_in) / seg_length;
        MYFLT radius_change = distance_to_seg_begin * radius_diff_per_meter;
        r = radius_in + radius_change;
    }
    area = r * r * PI;
    return area;
}

MYFLT parabolic_approx(MYFLT radius_in, MYFLT radius_out, MYFLT x,  \
                       MYFLT prelength, MYFLT length) {
  /*x is the current position(m*dx), prelength is the length of the total tube
  up until the current segment, length is the length of the current segment*/
  /*Formula for the parabola is taken as r=a*x^2+b, and variables a and b
  are determined using radius_in and radius_out*/
    MYFLT area;
    MYFLT a, b, r;
    if (length == 0) {
        r = radius_in;
    } else {
        b = radius_in;/*At x=0: a*0^2+b=b, b must be equal to initial radius */
        a = (radius_out-b)/(length*length);/*At x=length: a*length^2+b=radius_out. Solve wrt a */
        r = a*(x-prelength)*(x-prelength)+b;
    }
    area = r*r*PI;
    return area;
}

MYFLT exponential_approx(MYFLT radius_in, MYFLT radius_out, MYFLT x, MYFLT prelength, MYFLT length){
  /*x is the current position(m*dx), prelength is the length of the total tube
  up until the current segment, length is the length of the current segment*/
  /*Formula for the parabola is taken as r=a*e^(k*x), and variables a and k
  are determined using radius_in and radius_out*/
    MYFLT area;
    MYFLT a, k, r;
    if (length == 0) {
        r = radius_in;
    } else {
        a = radius_in;/*At x=0: a*e^0=a, so a must equal to initial radius */
        k = (log(radius_out/a)/length); /*At x=length: a*e^(k*length)=radius_out. Solve wrt k */
        r = a*exp(k*(x-prelength));
    }
    area = r*r*PI;
    return area;
}

void interp_loss(MYFLT rad, MYFLT coeff[4][100], MYFLT * r) {
    // TODO(AH): change to pointers??
  if (rad < 0.0035) {
    r[0] = coeff[0][0];
    r[1] = coeff[1][0];
    r[2] = coeff[2][0];
    r[3] = coeff[3][0];
    //    return r;
  }

  if (rad > 0.05) {
    r[0] = coeff[0][99];
    r[1] = coeff[1][99];
    r[2] = coeff[2][99];
    r[3] = coeff[3][99];
    //    return r;
  }

  int i = 1;
  while (radii[i] < rad) {
      i++;
    }

  r[0] = coeff[0][i-1]*(rad - radii[i])/(radii[i-1]-radii[i]) + coeff[0][i]*(radii[i-1]-rad)/(radii[i-1]-radii[i]);
  r[1] = coeff[1][i-1]*(rad - radii[i])/(radii[i-1]-radii[i]) + coeff[1][i]*(radii[i-1]-rad)/(radii[i-1]-radii[i]);
  r[2] = coeff[2][i-1]*(rad - radii[i])/(radii[i-1]-radii[i]) + coeff[2][i]*(radii[i-1]-rad)/(radii[i-1]-radii[i]);
  r[3] = coeff[3][i-1]*(rad - radii[i])/(radii[i-1]-radii[i]) + coeff[3][i]*(radii[i-1]-rad)/(radii[i-1]-radii[i]);
  //  return r;

}





void interpolation(int M, int Mold, MYFLT Lold, MYFLT dx,
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

void interpolation_visco(int M, int Mold, MYFLT Lold, MYFLT dx,  \
                            MYFLT dxold, MYFLT* klossnew, MYFLT* klossold) {
  MYFLT x, xl, xr;
  int m1, m2;
  for (int m = 0; m<= std::min(M,int(floor(Lold/dx))); m++){
    m1 = floor(m*(dx/dxold));
    m2 = m1 + 1; // Equivalent to m2 = ceil(m*dx/dxold)
    xl = m1*dxold;
    xr = m2*dxold;
    x  = m*dx;
    for (int k = 0; k < 4; k++){
      klossold[m*4+k] = klossnew[m1*4+k]*(x - xr)/(xl - xr) + klossnew[m2*4+k]*(x-xl)/(xr-xl);
    }

    for (int m = std::min(M,int(floor(Lold/dx)))+1; m<= M; m++){
      for (int k = 0; k < 4; k++){
        klossold[m*4+k] = klossnew[Mold*4+k];
      }
    }
  }
}


void interpolation_visco_arrays(int M, int Mold, MYFLT Lold, MYFLT dx, MYFLT dxold) {
    interpolation_visco(M, Mold, Lold, dx, dxold, iter_wloss, iter_wlossold);
    interpolation_visco(M, Mold, Lold, dx, dxold, iter_qloss, iter_qlossold);
}


MYFLT R0Z(MYFLT r, MYFLT rho_user, MYFLT eta){
    // constant R0 term in formular 6 (ISMA)
  MYFLT K = r*sqrt(rho_user/eta);
  MYFLT out  =  8*rho_user/(K*K);

  return out;
};

void compute_loss_arrays(int M, MYFLT* S, MYFLT RsZ[4][100], \
                                MYFLT LsZ[4][100], MYFLT GsY[4][100], \
                                MYFLT CsY[4][100], MYFLT dt, MYFLT rho_user, MYFLT c, \
                                MYFLT Zmult, MYFLT Ymult) {
  for (int m = 0; m <= M; m++) {
    MYFLT r  = sqrt(S[m]/PI);
    // read out constants or interpolate if radius inbetween
    interp_loss(r, RsZ, rz_tmp);
    interp_loss(r, LsZ, lz_tmp);
    interp_loss(r, GsY, gy_tmp);
    interp_loss(r, CsY, cy_tmp);
    sumZ1[m] = 0;  // needs to be interpolated from old value
    sumY1[m] = 0;

    for (int k = 0; k < 4; k++) {
      eLZ[m*4+k]  = exp(-lz_tmp[k]*dt);
      eCY[m*4+k]  = exp(-cy_tmp[k]*dt);

      sumZ1[m] += rz_tmp[k]*exp(-lz_tmp[k]*dt*.5);
      sumY1[m] += S[m]*gy_tmp[k] * exp(-cy_tmp[k]*dt*.5);

      MATZ[m*4+k]  =      Zmult*rz_tmp[k] * exp(-lz_tmp[k]*dt*.5);
      MATSY[m*4+k] = Ymult*S[m]*gy_tmp[k] * exp(-cy_tmp[k]*dt*.5);
    }

    factors_v[m] = rho_user/dt + sumZ1[m] + R0Z(r, rho_user, eta);  // formular 16 in ISMA
    factors_p[m] = S[m]/(rho_user*c*c*dt) + sumY1[m];          // formular 17
  }
}
