#include "../resonators/cone_radiation_losses_vp_PC.cpp"
#include "../resonators/cone_radiation_losses_vp_fixedM_Bela.cpp"
#include "../resonators/resonator_visco_concat.cpp"



/** Library loading */
// see Csound Book (2016) p.140 for Argument type overview
#include <modload.h>
void csnd::on_load(Csound *csound) {
  csnd::plugin<Cone_Radiation_Losses>(csound, "halfphysler", "aa", "akkkkkk", csnd::thread::ia);
  csnd::plugin<Cone_Radiation_Losses_fixedM_Bela>(csound, "halfphysler_bela", "aa", "akkkkkk", csnd::thread::ia);
  csnd::plugin<Resonator_Visco_Concat>(csound, "tube_resonator", "aa", "akk[]k[]k[]k[]PPPPo", csnd::thread::ia);
}
