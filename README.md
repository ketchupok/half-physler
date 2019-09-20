# Overview

The _Half-Physler_ as presented in NIME 2019 + the _resontube_ opcode as presented at ICSC 2019

This repository contains material to build a virtual single-reed instrument called the _Half-Physler_. In this design the excitation mechanism of the instrument is based on actual instrument parts (mouthpiece, clarinet reed, ligature) but the tube is simulated on the [Bela](http://www.bela.io) board.

In this repository you find the C++ code to compile the tube model as an opcode for [Csound](https://csound.com). The opcode is implemented following the guidelines of the _Csound Plugin Opcode Framework_ by Lazzarini (2017).


## Requirements
This is a Plug-In for [Csound](https://csound.com). So make sure you have the latest Csound (Version 6.10 and higher) installed.

The best way to install Csound is to follow the instructions on the [Csound Website](https://csound.com).



## Install the Halfphysler (general info)

### Mac/Linux
 ```
git clone https://github.com/ketchupok/half-physler.git

cd half-physler

make
```

### on Bela
It can be a bit difficult to access the Internet on the Bela Mini, as it has neither Wifi nor Ethernet on bord. On the Mac, you can share Internet with the Bela Mini plugged to your Mac via the USB port.

In short:
In Mac setting allow under 'Sharing' that internet is shared with Bela and Beaglbone.
(then the Bela IDE in the browser freezes)


In your terminal log into Bela:
```
myMac$> ssh root@bela.local

bela$> dhclient usb1

bela$> git clone https://github.com/ketchupok/half-physler.git

bela$> cd half-physler

bela$> make
```
[Find more details here.](https://www.hackster.io/hologram/sharing-internet-with-the-pocketbeagle-on-osx-cd62b2).

## Using the Half-physler Opcode

First run a simple demo.

### Mac

```
csound Demo_PC_Mac/demo_halfphysler.csd
```

### Linux / Ubuntu

```
csound Demo_PC_Mac/demo_halfphysler.csd
```


## Bela

This example needs an external signal! So connect a sound source to the audioIn0 of your Bela board. And a speaker to your right output, to hear some sound.

```
make installBela
belacsound --csd=/Bela/projects/Bela_HalfPhysler_DemoProject/_main.csd --period=32
```

In some cases we observed that building on Bela fails with error "plugin.h" not found. If you observe this, please report us.

### Windows
- currently not supported


# Csound Opcodes:

The opcodes implements a tube resonator with one closed and one
open end, similar to the resonator of a clarinet or a saxophone. An overview
of all user parameters of the opcode including a description of the underlying
physical parameters are given below.

## halfphysler

Resonator with radiation losses, driven by an initial air velocity.

Csound code:
```
   aFeedb, aSound halfphysler aVelocity, kLength, kRad, kSlope, kEndReflection, kDensity, kPos
```

- aVelocity      = input air velocity (signal) to drive the tube model
- kLength        = length of resonator in meters
- kRad           = radius of beginning section in meters
- kEndReflection = multiplier for end reflection
coefficient
- kDensity       = multiplier for air density
- kPos          = pickup position along the tube (0-1) relative to length (this only affects the aSound output, not the aFeedb!)

## halfphysler_bela

The halfphysler_bela version uses a fixed number of grid points for the model (M=32). This ensures that the no memory will be allocated in the real-time thread, and a number of gridpoints is used that does not cause dropouts on _Bela Mini_. Furthermore, the _feedback_ output is set to compensate for the latency of the _Bela Mini Board_ (--period=32).

## resontube

The _resontube_ opcode is the most flexible tube resonator opcode in this collection. It is not optimised to run on Bela and better suited to be used on Mac/Linux as an opcode in Csound. In comparison to the _halfphysler_ opcode it allows for more complex geometry settings and to switch on and off the calculation of viscothermal losses.

The resonator is driven by an input particle velocity (\texttt{aVelocity}), a parameter that describes the speed of the air entering the tube, for example via a single-reed instrument mouthpiece. The second input parameter \texttt{kLen} was introduced, to allow to change the resonance frequency of the resonator in a woodwind instrument-like style. A given \texttt{kLen} in meters cuts the resonator at this point, similar to the function of opening toneholes at acoustic woodwind instruments.

The initial geometry of the entire resonator is given in segments via the Csound arrays \texttt{kSegLengths[], kRadiiIn[], kRadiiOut[], kCurveType[]}. Each segment is defined by its length, input radius, output radius and an interpolation curve type. We allow up to 25 segments. The example below gives a good approximation of a Bb-flat clarinet geometry using only 4 segments.



```
kSegLengths[] fillarray 0.0316, 0.051, .3, 0.02
kRadiiIn[] fillarray 0.0055, 0.00635, 0.0075, 0.0075
kRadiiOut[] fillarray 0.0055, 0.0075, 0.0075, 0.0275
kCurveType[] fillarray 1, 1, 1, 2

aFeedb, aSnd resontube aVelocity, kLen, kSegLengths[], kRadiiIn[], kRadiiOut[], kCurveType[], [kEndReflect, kDensity, kPickPos, kComputeVisco]
```

Find a working example in:
```
csound Demo_PC_Mac/demo_resontube.csd
```

Additional parameters to shape the sound by modifying the end reflection, the air density and the pick-up position along the tube are provided as k-rate inputs. An option to switch between the computation of viscothermal losses or not was added, which allows real-time playback also for complex geometries and long resonators on consumer PCs. However, this will have a direct effect on the sound.



| Variable | type | Parameter | Range | Functionality |
|----------|:----:|----------:|-------|---------------|
|aVelocity | a-rate | v |  | Input particle velocity (m/s) |
|kLen | k-rate | L | 0.01–3 m (with sr = 44100) | length of the tube |
|kSegLengths[] | k-rate | L_part | sum =< 3 m | length of each segment, given as Csound array |
|kRadiiIn[] | k-rate | r_in  | 0.0075–0.0095 m (= 7.5--9.5 mm) | radius at the beginning of each segment, given as Csound array |
|kRadiiOut[] | k-rate | r_out|  0.0075–0.0095 m (= 7.5--9.5 mm) | end-radius of segment, given as Csound array |
|kCurveType[] | k-rate |  | 1)=linear, 2)=parabolic, 3)=exponential | interpolation mode for computation of S, given as Csound array |
|kEndReflection | k-rate | alpha | 0.1–4.0 | multiplier for end reflection coefficient (radiation resistance) |
|kDensity | k-rate | beta | 0.1–30.0 | multiplier for air density |
|kPickPos | k-rate |  | 0.0–1.0 | scaled pickup position along the tube |
|kComputeVisco | k-rate | | 0 / 1 | turn On/Off the computation of viscothermal losses |




## Further Information / References:

- S. Schmutzhard, V. Chatziioannou, and A. Hofmann. _Parameter optimisation of a viscothermal time-domain model for wind instruments._ In Proceedings of the 2017 International Symposium on Musical Acoustics, pages 27–30, Montreal, CA, 2017.

- V. Lazzarini. _The csound plugin opcode framework._ In Proceedings of the 14th Sound and Music Computing
Conference, pages 267–274, 2017.

- Hofmann, A., Chatziioannou, V., Schmutzhard, S., Erdoğan, G., & Mayer, A. (2019). _The Half-Physler: An oscillating real-time interface to a tube resonator model._ In Proceedings of the International Conference on New Interfaces for Musical Expression (NIME) 2019 (pp. 130–133). Porto Allegre, BR: Federal University of Rio Grande do Sul.
