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
ginit_EndRadius  = 0.0075
ginit_EndReflection  = 30
ginit_Density = 1.0
ginit_Pos     = 1.0
gk_bool  =  1
gk_Slope  = (ginit_EndRadius - ginit_Radius)/ginit_Length
FLpanel "The Half-Physler Parameter Panel", 780, 700, 50, 50 ;iwidth, iheight, ix, iy

    ;Display boxes for slider values
    ;                       iwidth, iheight,  ix,  iy
    idisp1 FLvalue " ",       100,      30,  540,  50
    idisp2 FLvalue " ",       100,      30,  540, 130
    idisp3 FLvalue "",        100,      30,  540, 210
    gidisp7 FLvalue "Slope ",        100,      30,  660, 210
    idisp4 FLvalue "",        100,      30,  540, 290
    idisp5 FLvalue "",        100,      30,  540, 370
    idisp6 FLvalue "",        100,      30,  540, 450

    ;FLslider puts a slider into the corresponding container.

    ;                                                                             imin,    imax, iexp, itype,      idisp, iwidth, iheight, ix,   iy
    gk_Length,         gihandle1  FLslider "Length",                               0.3,    1.20,    0,      5,    idisp1,   500,       30, 20,   50
    gk_Radius,         gihandle2  FLslider "Initial Radius",                       0.0035, 0.01,    0,      5,    idisp2,   500,       30, 20,  130
    gk_EndRadius,      gihandle3  FLslider "End Radius",                           0.001,  0.01,    0,      5,    idisp3,   500,       30, 20,  210
    gk_EndReflection,  gihandle4  FLslider "Reflection Parameter Multiplier",      0.1,   30.00,    0,      5,    idisp4,   500,       30, 20,  290
    gk_Density,        gihandle5  FLslider "Density Multiplier",                   0.1,   30.00,    0,      5,    idisp5,   500,       30, 20,  370
    gk_Pos,            gihandle6  FLslider "Relative Position",                    0.0,    1.00,    0,      5,    idisp6,   500,       30, 20,  450

    gkt1, iht1  FLbutton "Pling!",           1, 0, 2, 200, 100,  20, 550, 0, 3, 0, 100
    gkt2, iht2  FLbutton "Reset Parameters", 0, 0, 1, 200, 100, 230, 550, 0, 2, 0, 3
    gkt3, iht3 FLbutton "Fix Slope \n (End Radius inactive)",        1, 0, 2, 200, 100, 440, 550, -1

FLpanelEnd

FLrun

FLsetVal_i  ginit_Length,        gihandle1
FLsetVal_i  ginit_Radius,        gihandle2
FLsetVal_i  ginit_EndRadius,     gihandle3
FLsetVal_i  ginit_EndReflection, gihandle4
FLsetVal_i  ginit_Density,       gihandle5
FLsetVal_i  ginit_Pos,           gihandle6


instr 1



aIn mpulse 1, .5
aFeed, aSound halfphysler aIn, gk_Length, gk_Radius, gk_Slope, gk_EndReflection, gk_Density, gk_Pos
prints "PLAYING Half-Physler: Simple Cone with radiation losses with increasing the length.\n"
out aSound*0.002


endin


instr 2


; Set the widget's initial value
FLsetVal_i     ginit_Length,        gihandle1
FLsetVal_i     ginit_Radius,        gihandle2
FLsetVal_i     ginit_EndRadius,     gihandle3
FLsetVal_i     ginit_EndReflection, gihandle4
FLsetVal_i     ginit_Density,       gihandle5
FLsetVal_i     ginit_Pos,           gihandle6

; Capture when Slope should be kept fix
if (gkt3 == 1) then

   ; compute slope only once
   if (gk_bool  == 1) then
      gk_bool = 0
      gk_Slope  =  (gk_EndRadius - gk_Radius)/gk_Length
   endif

   ; Set gk_EndRadius accordingly
   FLsetVal  1, gk_Slope*gk_Length + gk_Radius, gihandle3

   ; Capture when EndRadius is too small
   if (gk_Slope*gk_Length + gk_Radius < 0.001) then
      FLsetVal  1, 0.001, gihandle3

      ; Set Slope accordingly
      gk_Slope  =  (gk_EndRadius - gk_Radius)/gk_Length

   endif
endif

; Display Slope
FLprintk 0.01, gk_Slope, gidisp7

; Capture when Slope is variable
if (gkt3 == 0) then
   gk_bool  = 1
   gk_Slope  =  (gk_EndRadius - gk_Radius)/gk_Length

endif


endin


instr 3
if (gkt1 == 0) then
  turnoff2 1,0,0
  turnoff2 3,0,0
elseif (gkt1 == 1) then
  event_i "i", 1, 0, 100
endif
endin



</CsInstruments>
<CsScore>

i2 0 3600

</CsScore>
</CsoundSynthesizer>
