<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so
-odac -b256 -B2048
;-otest.wav
</CsOptions>

<CsInstruments>

ksmps = 1
sr = 44100
0dbfs = 1

opcode FDelay,a,aki
setksmps 1
asig,kg,idel xin
 kpos init 0
 isize = idel > 1/sr ? round(idel*sr) : 1
 print isize
 adelay[] init isize
 kg = abs(kg) < 1 ? kg : 0
 printk 1, kg
 ay = asig + adelay[kpos]*kg
 xout adelay[kpos]
 adelay[kpos] = ay
 kpos = kpos == isize-1 ? 0 : kpos+1
endop

opcode allpass,a,aki
setksmps 1
asig,kg,idel xin
 isize = idel > 1/sr ? int(idel*sr) : 1
 adelay[] init isize
 kpos init 0
 ay = asig + adelay[kpos]*kg
 xout adelay[kpos] - kg*ay
 adelay[kpos] = ay
 kpos = kpos == isize-1 ? 0 : kpos+1
endop

opcode tonea,a,aii
 setksmps  1
 ain,icps,iDecay xin
 kcfenv    transeg icps*4,iDecay,-8,1,1,0,1
 aout      tone ain, kcfenv
 xout      aout
endop

opcode TubeDelay,a,aki
setksmps 1
asig,kg,idel xin
 kpos init 0
 isize = idel > 1/sr ? round(idel*sr) : 1
 print isize
 adelay[] init isize
 kg = abs(kg) < 1 ? kg : 0
 printk 1, kg
    /*
     icone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
     iradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
     iradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; slopes of cone segments
     icurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation
     */
     icone_lengths[] fillarray 1       ; lengths of cone segments [m]
     iradii_in[] fillarray 0.0075  ; radii of cone segments [m]
; is this slope or radius??
     iradii_out[]  fillarray 0.0075  ; radius of cone segment end
     icurve_type[] fillarray 1                     ; 1 = linear, 2 = parabolic; 3 = exponential approximation
     kPick init .9    ;  0.0 = left tube end - 1.0 = right tube end - default (0)
     kPick = 0.0
     kPick linseg 0.0, 1, 1.0
     ;kLength randomh .7, .5, 4
     kLength = 0.7
; resonator creates wrong sound
     ay_f, ay tube_resonator 0.03*asig+adelay[kpos]*kg, kLength, kPick, icone_lengths, iradii_in, iradii_out, icurve_type

; this resonator works
     kCylinder_Radius_m init 0.0075    ; initial radius of cone
     kSlope init 0.000                 ; -0.01 - 0.02  - default (0)
     kEndReflection init 1.0           ;  0.1 - 30      - default (1)
     kDensity init 1.0                 ;  0.1 - 30.0    - default (1)
     kPick init 0.0                     ;  0.0 = left tube end - 1.0 = right tube end - default (0)
     kLength_m init .7 ; =122Hz plays correctly
     ;kLength_m linseg .3,  1, .6
     ;ay, aysound halfphysler 0.03*asig+adelay[kpos]*kg, kLength_m, kCylinder_Radius_m, kSlope, kEndReflection, kDensity, kPick

 ;ay = asig + adelay[kpos]*kg
 xout adelay[kpos]
 ;ay tone ay, 200
 ay = 2 * taninv(ay) / 3.1415927
 adelay[kpos] = ay
 kpos = kpos == isize-1 ? 0 : kpos+1
endop

opcode AtanLimit, a, a
  ain xin
  aout = 2 * taninv(ain) / 3.1415927
  xout aout
endop


instr 1
    asig mpulse .5, 10
    ;asig  noise .3, 0
    ;asig     vco2 .3, 100
    ;kup linseg 0.0003, 4, 0.02
    ;adel FDelay asig, 0.45, 0.0024
    adel TubeDelay asig, 0.1, 0.0001
    aSound AtanLimit adel
    ;adel tonea asig, 400, 0.1
    out aSound

endin

instr 6
    icone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
    iradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
    iradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; slopes of cone segments
    icurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation
    aImpulse mpulse 1, 10
    kLength randomh .7, .5, 4
    ;a2 tube_resonator 0.03*aImpulse, kLength, icone_lengths, iradii_in, iradii_out, icurve_type
    a2 init 0
    prints "PLAYING Resonator_Visco_Concat_Pointers\n"
    out a2
endin

</CsInstruments>
<CsScore>
i1 0.5 2
;i6 0.5 200

</CsScore>
</CsoundSynthesizer>
