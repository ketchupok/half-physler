<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so
-odac -iadc -+rtmidi=alsa -M hw:1  --midi-velocity-amp=4 --midi-key-cps=5 -b256 -B2048
</CsOptions>

<CsInstruments>

ksmps = 32
sr = 44100
0dbfs = 1
gklength init  0.7
massign 1,1

instr 1
gklength = 20/p5
endin

instr 6
;aIn mpulse 1, .1
aIn inch 1
kCylinder_Radius_m init 0.0075    ; initial radius of cone 
kSlope init 0.00                 ; -0.01 - 0.02  - default (0)
kEndReflection init 1           ;  0.1 - 30      - default (1)
kDensity init 1.0                 ;  0.1 - 30.0    - default (1)
kPos init 0.0                     ;  0.0 = left tube end - 1.0 = right tube end - default (0)
aFeed, aSound halfphysler aIn, gklength, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPos
prints "PLAYING Half-Physler: Simple Cone with radiation losses with increasing the length.\n"
printks "%f", 0.1, gklength

;out  0.1*aSound

out aFeed*0.002
endin

</CsInstruments>
<CsScore>

i6 0.5 3600

</CsScore>
</CsoundSynthesizer>
