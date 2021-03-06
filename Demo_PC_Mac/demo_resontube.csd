<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so
-odac -b256 -B2048
</CsOptions>

<CsInstruments>

ksmps = 32
sr = 44100
0dbfs = 1

instr 1
    ; clarinet geometry
    kcone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
    kradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
    kradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; end-radius of cone segments
    kcurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation

    kLength linseg 0.2, 2, 0.3
    kPick_Pos = 1.0
    kEndReflection init 1.0           ;  0.1 - 30      - default (1)
    kEndReflection = 1.0
    kDensity = 1.0
    kComputeVisco = 1

    aImpulse mpulse .5, .1
    aFeedback, aSound resontube 0.005*aImpulse, kLength, kcone_lengths, kradii_in, kradii_out, kcurve_type, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    prints "PLAYING: Clarinet geometry & viscothermal losses\n"
    out aSound
endin


instr 2
    ; clarinet geometry no visco thermal losses
    kcone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
    kradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
    kradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; end-radius of cone segments
    kcurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation

    kLength linseg 0.2, 2, 0.3
    kPick_Pos = 1.0
    kEndReflection init 1.0           ;  0.1 - 30      - default (1)
    kEndReflection = 1.0
    kDensity = 1.0
    kComputeVisco = 0

    aImpulse mpulse .5, .1
    aFeedback, aSound resontube 0.005*aImpulse, kLength, kcone_lengths, kradii_in, kradii_out, kcurve_type, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    prints "PLAYING: Clarinet geometry without viscothermal losses\n"
    out aSound
endin

instr 3
    ; simple cylinder
    kcone_lengths[] fillarray 1.0       ; lengths of cone segments [m]
    kradii_in[] fillarray 0.0075        ; radii of cone segments [m]
    kradii_out[]  fillarray 0.0075      ; end-radius of cone segments
    kcurve_type[] fillarray 1           ; 1 = linear, 2 = parabolic; 3 = exponential approximation

    kLength linseg 0.2, 2, 0.3
    kPick_Pos = 1.0
    kEndReflection init 1.0           ;  0.1 - 30      - default (1)
    kEndReflection = 1.0
    kDensity = 1.0
    kComputeVisco = 1

    aImpulse mpulse .5, .1
    aFeedback, aSound resontube 0.005*aImpulse, kLength, kcone_lengths, kradii_in, kradii_out, kcurve_type, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    prints "PLAYING simple cylinder with viscothermal losses\n"
    out aSound
endin

instr 4
    ; simple cylinder no viscothermal losses
    kcone_lengths[] fillarray 1.0       ; lengths of cone segments [m]
    kradii_in[] fillarray 0.0075        ; radii of cone segments [m]
    kradii_out[]  fillarray 0.0075      ; end-radius of cone segments
    kcurve_type[] fillarray 1           ; 1 = linear, 2 = parabolic; 3 = exponential approximation

    kLength linseg 0.2, 2, 0.3
    kPick_Pos = 1.0
    kEndReflection init 1.0           ;  0.1 - 30      - default (1)
    kEndReflection = 1.0
    kDensity = 1.0
    kComputeVisco = 0

    aImpulse mpulse .5, .1
    aFeedback, aSound resontube 0.005*aImpulse, kLength, kcone_lengths, kradii_in, kradii_out, kcurve_type, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    prints "PLAYING simple cylinder without viscothermal losses\n"
    out aSound
endin

opcode WNoise, a, k
    kamp xin
    awhite unirand 2.0
    ; Normalize to +/-1.0
    awhite = awhite - 1.0
    xout awhite
endop

instr 5
    ; clarinet geometry
    kcone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
    kradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
    kradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; end-radius of cone segments
    kcurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation

    kLength linseg 0.2, 2, 0.3
    kPick_Pos = 1.0

    kEndReflection init 1.0           ;  0.1 - 30      - default (1)
    kEndReflection = 0.3
    kDensity = 1.0
    kComputeVisco = 1

    aWNoise WNoise 0.1
    aFeedback, aSound resontube 0.0007*aWNoise, kLength, kcone_lengths, kradii_in, kradii_out, kcurve_type, kEndReflection, kDensity, kPick_Pos, kComputeVisco
    prints "PLAYING: Clarinet geometry with viscothermal losses\n"
    out aSound
endin

</CsInstruments>
<CsScore>

i4 1 2 ; tube no visco
i3 4 2 ; tube + visco
i2 7 2 ; clari no visco
i1 10 2 ; clari + visco

i5 13 4

</CsScore>
</CsoundSynthesizer>
