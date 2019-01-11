<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so -b128
-odac -iadc -m0d
</CsOptions>
<CsInstruments>

ksmps = 16
sr = 44100

0dbfs = 1
nchnls  = 2

opcode AtanLimit, a, a
  ain xin
  aout = 2 * taninv(ain) / 3.1415927
  xout aout
endop

instr 6
  a0, a1 ins  ; we will only use the left input (0)
  a1 atone a1, 400
  kLength_m init 0.75
    aIn0 chnget "analogIn0"
    kLength_m = (k(aIn0)*0.5) + 0.35

  kCylinder_Radius_m init 0.0075    ; initial radius of cone
    ; Analog in 1 controls radius in range 0.0075, 0.0095
    aIn1 chnget "analogIn1"
    kCylinder_Radius_m = (k(aIn1)*0.002) + 0.0075

  kSlope init 0.00 ; we can/should allow negative slope, up to -r/0.75
    ; Analog in 2
    aIn2 chnget "analogIn2"  ; 0 - 0.87
    kSlope = (k(aIn2 - 0.3)*0.03)

  kEndReflection init 1.0           ; 0.1 - 2 - default (1)
    ; Analog in 3
    aIn3 chnget "analogIn3"
    kEndReflection = (k(aIn3)*2) + 0.1

  kDensity init 1.0                 ; 0.1 - 30.0 - default (1)
    ; Analog in 4
    aIn4 chnget "analogIn4"
    kDensity = (k(aIn4)*30) + 0.1
    ;printk2 kDensity

  kFeedback init 1.0                 ; 0.1 - 1.0 - default (1)
    ; Analog in 5
    aIn5 chnget "analogIn5"
    kFeedback = (k(aIn5)*1.3) + 0.01

  kPos init 0.0                     ; 0.0 = tube end - 1.0 = half tube
  aFeed, aSound halfphysler_bela a0*kFeedback, kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPos

  prints "PLAYING Half-Physler: Cone_Radiation_Losses\n"
  aFeed = aFeed * 0.5
  aFeed AtanLimit aFeed
  aIm mpulse 1, .1
 aSound = aSound *0.1
 aSound AtanLimit aSound
  outs aFeed, aSound    ; left output goes to actuator, right output to speaker
endin

</CsInstruments>
<CsScore>

i6 0.5 100000

</CsScore>
</CsoundSynthesizer>
