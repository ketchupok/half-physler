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
    iradii_out[]  fillarray 0.0055, 0.0075, 0.0075, 0.0275  ; end-radius of cone segments
    icurve_type[] fillarray 1, 1, 1, 2                      ; 1 = linear, 2 = parabolic; 3 = exponential approximation
    aImpulse mpulse 1, .1
    kLength randomh .7, .5, 4
    kPick_Pos init 0.0
    aFeedback, aSound tube_resonator 0.005*aImpulse, kLength, kPick_Pos, icone_lengths, iradii_in, iradii_out, icurve_type
    prints "PLAYING Resonator_Visco_Concat_Pointers\n"
    out aSound
endin

</CsInstruments>
<CsScore>

i6 0.5 2

</CsScore>
</CsoundSynthesizer>
