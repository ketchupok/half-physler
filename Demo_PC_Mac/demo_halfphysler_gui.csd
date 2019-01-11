<CsoundSynthesizer>
<CsOptions>
--opcode-lib=./opcodes_iwk.so
-odac -b256 -B2048
</CsOptions>

<CsInstruments>

ksmps = 32
sr = 44100
0dbfs = 1


ginit_Length  = 0.3
ginit_Radius  = 0.0075
ginit_Slope   = 0
ginit_EndReflection  = 30
ginit_Density = 1.0
ginit_Pos     = 1.0

FLpanel "The Half-Physler Parameter Panel", 660, 700, 50, 50 ;iwidth, iheight, ix, iy

    ;Display boxes for slider values
    ;                       iwidth, iheight,  ix,  iy
    idisp1 FLvalue " ",       100,      30,  540,  50
    idisp2 FLvalue " ",       100,      30,  540, 130
    idisp3 FLvalue "",        100,      30,  540, 210
    idisp4 FLvalue "",        100,      30,  540, 290
    idisp5 FLvalue "",        100,      30,  540, 370
    idisp6 FLvalue "",        100,      30,  540, 450
 
    ;FLslider puts a slider into the corresponding container.
    
    ;                                                                             imin,    imax, iexp, itype,      idisp, iwidth, iheight, ix,   iy
    gk_Length,         gihandle1  FLslider "Length",                               0.3,    1.20,    0,      5,    idisp1,   500,       30, 20,   50
    gk_Radius,         gihandle2  FLslider "Initial Radius",                       0.0035, 0.01,    0,      5,    idisp2,   500,       30, 20,  130
    gk_Slope,          gihandle3  FLslider "Slope",                               -0.01,   0.02,    0,      5,    idisp3,   500,       30, 20,  210
    gk_EndReflection,  gihandle4  FLslider "Reflection Parameter Multiplier",      0.1,   30.00,    0,      5,    idisp4,   500,       30, 20,  290
    gk_Density,        gihandle5  FLslider "Density Multiplier",                   0.1,   30.00,    0,      5,    idisp5,   500,       30, 20,  370
    gk_Pos,            gihandle6  FLslider "Relative Position",                    0.0,    1.00,    0,      5,    idisp6,   500,       30, 20,  450

    gkt1, iht1 FLbutton "Pling!",           0, 0, 1, 200, 100,  20, 550, 0, 1, 0, 100
    gkt2, iht2 FLbutton "Reset Parameters", 0, 0, 1, 200, 100, 230, 550, 0, 2, 0, 3
    gkt3, iht3 FLbutton "Stop",             0, 0, 1, 200, 100, 440, 550, 0, 3, 0, 0.001	
    
FLpanelEnd

FLrun


FLsetVal_i  ginit_Length,        gihandle1
FLsetVal_i  ginit_Radius,        gihandle2
FLsetVal_i  ginit_Slope,         gihandle3
FLsetVal_i  ginit_EndReflection, gihandle4
FLsetVal_i  ginit_Density,       gihandle5
FLsetVal_i  ginit_Pos,           gihandle6


instr 1
aIn mpulse 1, .5
aFeed, aSound halfphysler aIn, gk_Length, gk_Radius, gk_Slope, gk_EndReflection, gk_Density, gk_Pos
prints "PLAYING Half-Physler: Simple Cone with radiation losses with increasing the length.\n"
out aFeed*0.002
endin


instr 2


; Set the widget's initial value
FLsetVal_i     ginit_Length,        gihandle1
FLsetVal_i     ginit_Radius,        gihandle2
FLsetVal_i     ginit_Slope,         gihandle3
FLsetVal_i     ginit_EndReflection, gihandle4
FLsetVal_i     ginit_Density,       gihandle5
FLsetVal_i     ginit_Pos,           gihandle6
endin


instr 3
turnoff2 1,0,0
endin



</CsInstruments>
<CsScore>

i2 0 3600

</CsScore>
</CsoundSynthesizer>
