<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so
-odac1 -b256 -B2048 -Ma
;-otest.wav
</CsOptions>

<CsInstruments>

ksmps = 32
sr = 44100
0dbfs = 1

opcode FeedbackHalfphysler,a,akikkkkkk
setksmps 1
asig, kFbk, idel,kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPick xin
 kpos init 0
 isize = idel > 1/sr ? round(idel*sr) : 1
 adelay[] init isize
 kFbk = abs(kFbk) < 1 ? kFbk : 0
     kCylinder_Radius_m init 0.0075    ; initial radius of cone
     kSlope init 0.000                 ; -0.01 - 0.02  - default (0)
     kEndReflection init 1.0           ;  0.1 - 30      - default (1)
     kDensity init 1.0                 ;  0.1 - 30.0    - default (1)
     kPick init 0.0                    ;  0.0 = left tube end - 1.0 = right tube end - default (0)
     kLength_m init .7 ; =122Hz plays correctly
     ;kLength_m linseg .3,  1, .6
     ay, aSound halfphysler 0.03*asig+adelay[kpos]*kFbk, kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPick
 ;xout adelay[kpos]
 xout aSound

 ay = 2 * taninv(ay) / 3.1415927 ; limiter
 adelay[kpos] = ay
 kpos = kpos == isize-1 ? 0 : kpos+1
endop

opcode AtanLimit, a, a
  ain xin
  aout = 2 * taninv(ain) / 3.1415927
  xout aout
endop

opcode Saw, a, kk
    kfreq, kamp xin
    asig     vco2 kamp, kfreq
    xout asig
endop


instr 1

    aImpulse mpulse .5, 1000
    kLength_m           ctrl7 1, 41, 0.3, 0.9
    kCylinder_Radius_m  ctrl7 1, 42, 0.0075, 0.0095
    kSlope              ctrl7 1, 43, -0.007, 0.007
    kEndReflection      ctrl7 1, 44, 0.1, 4.0
    kDensity            ctrl7 1, 45, 0.1, 30.0
    kPick_Pos           ctrl7 1, 46, 0.0, 1.0
    ;printk 0.5, kLength_m
    kFeedback           ctrl7 1, 47, 0.001, 0.005
    iDelay init 0.3
    aL FeedbackHalfphysler aImpulse, kFeedback, 0.0003, kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPick_Pos
    aL AtanLimit aL
    out aL
endin


</CsInstruments>
<CsScore>
i1 0.5 200
;i6 0.5 200

</CsScore>
</CsoundSynthesizer>
