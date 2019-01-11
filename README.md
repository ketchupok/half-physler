# Overview

The _Half-Phylser_ as presented in NIME 2019

This repository contains material to build a virtual single-reed instrument called the _Half-Physler_. In this design the excitation mechanism of the instrument is based on actual instrument parts (mouthpiece, clarinet reed, ligature) but the tube is simulated on the [Bela](http://www.bela.io) board.

In this repository you find the C++ code to compile the tube model as an opcode for [Csound](https://csound.com). The opcode is implemented following the guidlines of the _Csound Plugin Opcode Framwork_ by Lazzarini (2017).



## Install
 git clone https://github.com/ketchupok/half-physler.git

 ```
 cd half-physler

 make
```
## Usage

### Mac

```
csound Demo_PC_Mac/demo_halfphysler.csd
```

### Linux / Ubuntu

```
csound Demo_PC_Mac/demo_halfphysler.csd
```


## Bela

```
make installBela
belacsound --csd=/Bela/projects/Bela_HalfPhysler_DemoProject/_main.csd --period=32
```

### Windows
- currently not supported


# Csound Opcodes:

- halfphysler

Resonator with radiation losses, driven by an initial air velocity.

Csound code:
```
   aFeedb, aSound halfphysler aVelocity, kLength, kRad, kSlope, kEndReflection, kDensity
```

aVelocity      = input signal to drive the tube model
kLength        = length of resonator in meters
kRad           = radius of beginning section in meters
kEndReflection = multiplier for end reflection
coefficient
kDensity       = multiplier for air density


- halfphysler_bela

The halfphysler_bela version uses a fixed number of grid points for the model (M=32). This ensures that the no memory will be allocated in the real-time thread, and a number of gridpoints is used that does not cause dropouts on _Bela Mini_. Furthermore, the _feedback_ output is set to compensate for the latency of the _Bela Mini Board_ (--period=32).

## Further Information / References:

- S. Schmutzhard, V. Chatziioannou, and A. Hofmann. _Parameter optimisation of a viscothermal time-domain model for wind instruments._ In Proceedings of the 2017 International Symposium on Musical Acoustics, pages 27–30, Montreal, CA, 2017.

- V. Lazzarini. _The csound plugin opcode framework._ In Proceedings of the 14th Sound and Music Computing
Conference, pages 267–274, 2017.
