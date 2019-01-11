#ifndef CONST_H
#define CONST_H
#include <plugin.h>


extern MYFLT radii[100]; // radii for the interpolation of viscothermal losses coefficients
extern MYFLT GsY[4][100];   // precalculated G-coeffs of viscothermal losses of Y at radii (K = 4)
extern MYFLT CsY[4][100];   // precalculated C-coeffs of viscothermal losses of Y at radii (K = 4)
extern MYFLT RsZ[4][100];   // precalculated R-coeffs of viscothermal losses of Z at radii (K = 4)
extern MYFLT LsZ[4][100];   // precalculated L-coeffs of viscothermal losses of Z at radii (K = 4)


extern MYFLT c;
extern MYFLT rho;
extern MYFLT eta;

extern int   Mmax;

extern MYFLT rad_alpha;
extern MYFLT rad_beta;


#endif
