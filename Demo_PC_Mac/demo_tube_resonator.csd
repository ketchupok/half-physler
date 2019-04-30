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
    /*
    icone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
    iradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
    iradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; end-radius of cone segments
    icurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation
    */

    kcone_lengths[] fillarray .1, .1, .1, .1       ; lengths of cone segments [m]
    kradii_in[] fillarray 0.0055, 0.0055, 0.0075, 0.0075   ; radii of cone segments [m]
    kradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; end-radius of cone segments
    kcurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation

    aImpulse mpulse .5, .1
    ;kLength randomh .7, .5, 4
    ;kLength init .3
    kLength linseg 0.2, 3, 0.3

    ;kPick_Pos init 0.0
    kPick_Pos = 0.4
    kEndReflection init 1.0           ;  0.1 - 30      - default (1)
    kEndReflection = 1.0

    ;kDensity init 1.0                 ;  0.1 - 30.0    - default (1)
    kDensity = 1.0

    ;kComputeVisco init 1
    kComputeVisco = 1

    aFeedback, aSound tube_resonator 0.005*aImpulse, kLength, kPick_Pos, kEndReflection, kDensity, kcone_lengths, kradii_in, kradii_out, kcurve_type, kComputeVisco
    prints "PLAYING Resonator_Visco_Concat\n"
    out aSound
endin

</CsInstruments>
<CsScore>

i6 0.5 4

</CsScore>
</CsoundSynthesizer>
