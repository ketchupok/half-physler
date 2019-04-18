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
    icone_lengths[] fillarray 0.0316, 0.051, .3, 0.2       ; lengths of cone segments [m]
    iradii_in[] fillarray 0.0055, 0.00635, 0.0075, 0.0075   ; radii of cone segments [m]
    iradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; ??slopes of cone segments?? or radius??
    icurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation
    aImpulse mpulse 1, .1
    k1 randomh .7, .5, 4
    a2 tube_resonator 0.005*aImpulse, k1, icone_lengths, iradii_in, iradii_out, icurve_type
    prints "PLAYING Resonator_Visco_Concat_Pointers\n"
    out a2
endin

</CsInstruments>
<CsScore>

i6 0.5 2

</CsScore>
</CsoundSynthesizer>
