#include "../resonators/cone_radiation_losses_vp_pointers_PC.cpp"
#include "../resonators/cone_radiation_losses_vp_pointers_fixedM_Bela.cpp"
#include "../resonators/resonator_visco_concat_pointers.cpp"



/** Library loading */
#include <modload.h>
void csnd::on_load(Csound *csound) {
  csnd::plugin<Cone_Radiation_Losses>(csound, "halfphysler", "aa", "akkkkkk", csnd::thread::ia);
  csnd::plugin<Cone_Radiation_Losses_fixedM_Bela>(csound, "halfphysler_bela", "aa", "akkkkkk", csnd::thread::ia);
  csnd::plugin<Resonator_Visco_Concat_Pointers>(csound, "tube_resonator", "aa", "akki[]i[]i[]i[]", csnd::thread::ia);
}
