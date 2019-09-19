<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so
;-odac1 -b256 -B2048 -Ma
;-odac -iadc -+rtmidi=alsa -Ma  -b512 -B2048
-odac -iadc -b512 -B2048
;-otest.wav
</CsOptions>

<CsInstruments>

ksmps = 32
sr = 22050
0dbfs = 1

opcode FeedbackTube,a,akikkkkkkk
setksmps 1
asig, kFbk, idel,kLength_m, kCylinder_Radius_m, kR_out, kEndReflection, kDensity, kPick_Pos, kComputeVisco xin
 kpos init 0
 isize = idel > 1/sr ? round(idel*sr) : 1
 adelay[] init isize
 kFbk = abs(kFbk) < 1 ? kFbk : 0
     kCylinder_Radius_m init 0.0075    ; initial radius of cone
     kR_out init 0.0075                ; -0.01 - 0.02  - default (0)
     kEndReflection init 1.0           ;  0.1 - 4      - default (1)
     kDensity init 1.0                 ;  0.1 - 30.0    - default (1)
     kPick_Pos init 0.0                    ;  0.0 = left tube end - 1.0 = right tube end - default (0)
     kLength_m init .7 ; =122Hz plays correctly
     ;kLength_m linseg .3,  1, .6

     kcone_lengths[] fillarray .5, .5       ; lengths of cone segments [m]
     kradii_in[] fillarray 0.0075  ; radii of cone segments [m]
     kradii_out[]  fillarray 0.0075   ; radius of cone segment end
     kcurve_type[] fillarray 1
     ;kComputeVisco = 1
     kradii_out[0] = kR_out
     kradii_in[0] = kCylinder_Radius_m

     ;printk 1, kradii_out[0]
     ;ay, aSound halfphysler 0.03*asig+adelay[kpos]*kFbk, kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPick
     ;kDensity = 1.0 ; keep 1.0 when visco ON
     aFeedback, aSound tube_resonator 0.03*asig+adelay[kpos]*kFbk, kLength_m, kcone_lengths, kradii_in, kradii_out, kcurve_type, kEndReflection, kDensity, kPick_Pos, kComputeVisco
 ;xout adelay[kpos]
 xout aSound

 aFeedback = 2 * taninv(aFeedback) / 3.1415927 ; limiter
 adelay[kpos] = aFeedback
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
    kLength_m           ctrl7 1, 21, 0.3, 0.9
    kCylinder_Radius_m  ctrl7 1, 22, 0.0075, 0.0095
    kSlope              ctrl7 1, 23, 0.0035, 0.0135
    kEndReflection      ctrl7 1, 24, 0.1, 4.0
    kDensity            ctrl7 1, 25, 0.5, 30.0
    kPick_Pos           ctrl7 1, 26, 0.0, 1.0
    ;printk 0.5, kLength_m
    kFeedback           ctrl7 1, 27, 0.00001, 0.005
    kComputeVisco       ctrl7 1, 28, 0, 1.0


    ;kLength_m port kLength_m, 0.01
    ;kLength_m linseg 0.1, 20, 0.9
    kCylinder_Radius_m port kCylinder_Radius_m, 0.1
    ;kSlope port kSlope, 0.1
    kEndReflection port kEndReflection, 0.1
    ;kDensity port kDensity, 0.1
    kPick_Pos port kPick_Pos, 0.1
    kFeedback port kFeedback, 0.1

    aL FeedbackTube aImpulse, kFeedback, 0.0003, kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    aL AtanLimit aL
    out aL
endin

instr 2

    aImpulse mpulse .1, 1000
    kLength_m           = 0.5;
    kCylinder_Radius_m  = 0.0075 ;ctrl7 1, 22, 0.0075, 0.0095
    kR_out              = 0.0035 ; ctrl7 1, 23, 0.0035, 0.0135
    kR_out linseg 0.0075, 2, 0.0075, 3, 0.001, 1, 0.003, 3, 0.009
    kEndReflection      = 0.3 ; ctrl7 1, 24, 0.1, 4.0
    kDensity            = 1.0; ctrl7 1, 25, 0.5, 30.0
    kPick_Pos           = 1.0; ctrl7 1, 26, 0.0, 1.0
    ;printk 0.5, kLength_m
    kFeedback           = 0.001; ctrl7 1, 27, 0.00001, 0.005
    kComputeVisco       = 1; ctrl7 1, 28, 0, 1.0


    ;kLength_m port kLength_m, 0.01
    ;kLength_m linseg 0.1, 20, 0.9
    kCylinder_Radius_m port kCylinder_Radius_m, 0.1
    kR_out port kR_out, 0.01
    kEndReflection port kEndReflection, 0.1
    ;kDensity port kDensity, 0.1
    kPick_Pos port kPick_Pos, 0.1
    kFeedback port kFeedback, 0.1

    aL FeedbackTube aImpulse, kFeedback, 0.0003, kLength_m, kCylinder_Radius_m, kR_out, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    aL AtanLimit aL
    out aL
endin

opcode Sine, a, kk
    kfreq, kamp xin
    asig     poscil kamp, kfreq
    xout asig
endop


instr 3
; slow changes
    aImpulse mpulse .1, 1000
    kLength_m           = 0.5;
    kLength_m  linseg 0.5, 20, 0.9

    kCylinder_Radius_m  = 0.0075 ;ctrl7 1, 22, 0.0075, 0.0095
    aCyi Sine 0.04, 0.004
    kCylinder_Radius_m = aCyi + 0.0075

    kR_out              = 0.0035 ; ctrl7 1, 23, 0.0035, 0.0135
    ;kR_out linseg 0.0075, 2, 0.0075, 3, 0.001, 1, 0.003, 3, 0.009
    aR_out Sine 0.03, 0.004
    kR_out = (aR_out) + 0.0075

    kEndReflection      = 0.3 ; ctrl7 1, 24, 0.1, 4.0
    kDensity            = 1.0; ctrl7 1, 25, 0.5, 30.0
    kPick_Pos           = 1.0; ctrl7 1, 26, 0.0, 1.0
    ;apick Sine 0.07, 0.4
    ;kPick_Pos = (apick) + 0.5
    ;printk 0.5, kLength_m
    kFeedback           = 0.001; ctrl7 1, 27, 0.00001, 0.005
    kComputeVisco       = 1; ctrl7 1, 28, 0, 1.0


    ;kLength_m port kLength_m, 0.01
    ;kLength_m linseg 0.1, 20, 0.9
    kCylinder_Radius_m port kCylinder_Radius_m, 0.1
    kR_out port kR_out, 0.01
    kEndReflection port kEndReflection, 0.1
    ;kDensity port kDensity, 0.1
    kPick_Pos port kPick_Pos, 0.1
    kFeedback port kFeedback, 0.1

    aL FeedbackTube aImpulse, kFeedback, 0.0003, kLength_m, kCylinder_Radius_m, kR_out, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    aL AtanLimit aL
    out aL
endin


</CsInstruments>
<CsScore>
;i1 0.5 200
;i2 0.5 200
i3 0.5 200
;i6 0.5 200

</CsScore>
</CsoundSynthesizer>
