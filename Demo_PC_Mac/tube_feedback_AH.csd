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

opcode TubeFeedbackClarinet,a,aki
setksmps 1
asig,kg,idel xin
 kpos init 0
 isize = idel > 1/sr ? round(idel*sr) : 1
 adelay[] init isize
 kg = abs(kg) < 1 ? kg : 0

     icone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
     iradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
     iradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; slopes of cone segments
     icurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation
     /*
     icone_lengths[] fillarray 1       ; lengths of cone segments [m]
     iradii_in[] fillarray 0.0075  ; radii of cone segments [m]
     iradii_out[]  fillarray 0.0075  ; radius of cone segment end
     icurve_type[] fillarray 1                     ; 1 = linear, 2 = parabolic; 3 = exponential approximation
     */
     kPick init .9    ;  0.0 = left tube end - 1.0 = right tube end - default (0)
     kPick = 0.0
     kLength init 0.35
     ;kLength randomh .7, .5, 8
     ;kLength = 0.7
     ay_f, ay tube_resonator 0.03*asig+adelay[kpos]*kg, kLength, kPick, icone_lengths, iradii_in, iradii_out, icurve_type
  xout adelay[kpos]
 ;ay tone ay, 200
 ;ay = 2 * taninv(ay) / 3.1415927
     adelay[kpos] = ay
     kpos = kpos == isize-1 ? 0 : kpos+1
endop

opcode TubeFeedbackCylinder,a,aki
setksmps 1
asig,kg,idel xin
 kpos init 0
 isize = idel > 1/sr ? round(idel*sr) : 1
 adelay[] init isize
 kg = abs(kg) < 1 ? kg : 0

     icone_lengths[] fillarray 1       ; lengths of cone segments [m]
     iradii_in[] fillarray 0.0075  ; radii of cone segments [m]
     iradii_out[]  fillarray 0.0075  ; radius of cone segment end
     icurve_type[] fillarray 1                     ; 1 = linear, 2 = parabolic; 3 = exponential approximation

     kPick init .9    ;  0.0 = left tube end - 1.0 = right tube end - default (0)
     kPick = 0.0
     kLength init 0.35 ;235Hz
     kLength init 0.7  ;118Hz
     ;kLength randomh .7, .5, 8
     ;kLength = 0.7
     ay_f, ay tube_resonator 0.03*asig+adelay[kpos]*kg, kLength, kPick, icone_lengths, iradii_in, iradii_out, icurve_type
  xout adelay[kpos]
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
    ;adel TubeFeedback asig, 0.00012, 0.0003
    adel TubeFeedbackCylinder asig, 0.00007, 0.0003
    aSound AtanLimit adel
    out aSound

endin


</CsInstruments>
<CsScore>
i1 0.5 10
;i6 0.5 200

</CsScore>
</CsoundSynthesizer>
