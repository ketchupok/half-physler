<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so
-odac -b256 -B2048
</CsOptions>

<CsInstruments>

ksmps = 32
sr = 44100
0dbfs = 1

instr 6
aIn mpulse 1, .1
kLength_m init .3
kLength_m linseg .3,  1, .6
kCylinder_Radius_m init 0.0075    ; initial radius of cone
kSlope init 0.00                  ; 0.0 - 0.01
kEndReflection init 1.0           ; 0.1 - 2 - default (1)
kDensity init 1.0                 ; 0.1 - 30.0 - default (1)
kPos init 0.0                     ; 0.0 = tube end - 1.0 = half tube
aFeed, aSound halfphysler aIn, kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPos
prints "PLAYING Half-Physler: Simple Cone with radiation losses with increasing the length.\n"
out aFeed*0.002
endin

</CsInstruments>
<CsScore>

i6 0.5 2

</CsScore>
</CsoundSynthesizer>
