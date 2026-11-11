* defs.h
* Header file
* This file is part of LoopInts

***[for debugging
c#define DEBUG1
c#define DEBUG2
c#define DEBUG5

c#define NO_CACHE
***]for debugging

***[possible option for pseudo-coll singular cases
* only for 3-point and 4-point scalar integrals
c#define SCALAR_COLL_QUAD
***]possible option for pseudo-coll singular cases

#define CVAR double complex
#define RVAR double precision

#define CVAR32 complex*32
#define RVAR16 real*16

#ifdef QUADPREC

#define QPREC(x) QEXT(x)

#define YPDefault DefaultQP
#define YPSetPrec SetPrecQP
#define YPffexi ffexiQP
#define YPPREGUL PREGULQP
#define YPCayleyE CayleyEQP
#define YPCayleyF CayleyFQP
#define YPparanames paranamesQP
#define YPcoeffnames coeffnamesQP
#ifdef NO_CACHE
#define YPltiget ltigetQP
#endif

#define YPclearcache clearcacheQP
#define YPmarkcache markcacheQP
#define YPrestorecache restorecacheQP
#define YPsetmudim setmudimQP
#define YPgetmudim getmudimQP
#define YPsetdelta setdeltaQP
#define YPgetdelta getdeltaQP
#define YPsetlambda setlambdaQP
#define YPgetlambda getlambdaQP
#define YPsetstab0 setstab0QP
#define YPsetstab2 setstab2QP
#define YPsetstab3 setstab3QP
#define YPsetstab4 setstab4QP
#define YPsetstab5 setstab5QP
#define YPsetstab6 setstab6QP
#define YPsetdebugkey setdebugkeyQP
#define YPsetversionkey setversionkeyQP
#define YPsetmaxdev setmaxdevQP
#define YPsetrunmode setrunmodeQP
#define YPNameData NameDataQP

#define CYPLOG LOGCQP
#define CYPeta etaCQP
#define CYPetax etaxCQP
#define CYPRoot2 Root2CQP
#define CYPSp SpCQP
#define YPLi2 Li2QP
#define CYPSp1 Sp1CQP
#define CYPLi2n Li2nCQP
* Ninh working
#define CYPLi2nx Li2nxCQP
#define CYPLambda LambdaCQP
#define CYPKam KamCQP
#define YPirk irkQP
#define YPDltB DltBQP
#define YPErrId ErrIdQP

#define CYPDecomp DecompCQP
#define CYPSolve SolveCQP
#define CYPDet DetCQP
#define CYPInverse InverseCQP

#define CYPA0s A0sCQP
#define CYPA00s A00sCQP
#define CYPBcoeff BcoeffCQP
#define CYPCcoeff CcoeffCQP
#define CYPDcoeff DcoeffCQP
#define CYPEcoeff EcoeffCQP
#define CYPFcoeff FcoeffCQP
#define CYPDumpPara DumpParaCQP
#define CYPDumpCoeff DumpCoeffCQP
#define CYPBscalar BscalarCQP
#define CYPBtensor BtensorCQP
#define CYPBnx BnxCQP
#define CYPFnx FnxCQP
#define CYPdBnx dBnxCQP
#define CYPGnx GnxCQP
#define CYPHnx HnxCQP
#define CYPKnx KnxCQP
#define CYPLnx LnxCQP
#define CYPCscalar CscalarCQP
#define CYPC0ir C0irCQP
#define CYPCtensor CtensorCQP
#define CYPRotate RotateCQP
#define CYPRotateP5 RotateP5CQP
#define CYPSetM0 SetM0CQP
#define CYPRotateM0 RotateM0CQP
#define CYPRotate24 Rotate24CQP
#define CYPD0_abc D0_abcCQP
#define CYPDscalar DscalarCQP
#define CYPD0per D0perCQP
#define CYPD0ir_3m0q_2p0a D0ir_3m0q_2p0aCQP
#define CYPD0ir_3m0q_2p0b D0ir_3m0q_2p0bCQP
#define CYPD0ir_2m0q_0p0a D0ir_2m0q_0p0aCQP
#define CYPD0ir_2m0q_0p0b D0ir_2m0q_0p0bCQP
#define CYPD0ir_1m0q_2p0a D0ir_1m0q_2p0aCQP
#define CYPD0ir_0m0q_4p0a D0ir_0m0q_4p0aCQP
#define CYPD0ir_1m0q_0p0a D0ir_1m0q_0p0aCQP
#define CYPD0ir_1m0q_0p0b D0ir_1m0q_0p0bCQP
#define CYPD0ir_1m0q_0p0c D0ir_1m0q_0p0cCQP

#define CYPD0ir_4m0_4p0a D0ir_4m0_4p0aCQP
#define CYPD0ir_4m0_3p0a D0ir_4m0_3p0aCQP
#define CYPD0ir_4m0_2p0a D0ir_4m0_2p0aCQP
#define CYPD0ir_3m0_2p0a D0ir_3m0_2p0aCQP
#define CYPD0ir_3m0_2p0b D0ir_3m0_2p0bCQP
#define CYPD0ir_3m0_2p0c D0ir_3m0_2p0cCQP
#define CYPD0ir_2m0_1p0a D0ir_2m0_1p0aCQP
#define CYPD0ir_2m0_1p0b D0ir_2m0_1p0bCQP
#define CYPD0ir_2m0_1p0c D0ir_2m0_1p0cCQP
#define CYPD0ir_2m0_0p0a D0ir_2m0_0p0aCQP
#define CYPD0ir_1m0_0p0a D0ir_1m0_0p0aCQP
#define CYPD0M0 D0M0CQP
#define CYPD0DD D0DDCQP
#define CYPDtensor DtensorCQP
#define CYPEscalar EscalarCQP
#define CYPEtensor EtensorCQP
* Ninh checking
#define CYPEscalar1 Escalar1CQP
#define CYPEtensor1 Etensor1CQP
#define CYPEscalar2 Escalar2CQP
#define CYPEtensor2 Etensor2CQP

#define CYPFscalar FscalarCQP
#define CYPFtensor FtensorCQP

#define CYPfT fTCQP
#define CYPtriangle triangleCQP
#define CYPfT_lin fT_linCQP
#define CYPfT13 fT13CQP
#define CYPfS3n fS3nCQP
#define CYPfS3n_lin fS3n_linCQP
#define CYPfS2 fS2CQP
#define CYPfS2_lin fS2_linCQP
#define CYPfS2_linr fS2_linrCQP
#define CYPfRn fRnCQP
#define CYPfR fRCQP

#else
#define QPREC(x) x

#define YPDefault DefaultDP
#define YPSetPrec SetPrecDP
#define YPffexi ffexiDP
#define YPPREGUL PREGULDP
#define YPCayleyE CayleyEDP
#define YPCayleyF CayleyFDP
#define YPparanames paranamesDP
#define YPcoeffnames coeffnamesDP
#ifdef NO_CACHE
#define YPltiget ltigetDP
#endif

#define YPclearcache clearcacheDP
#define YPmarkcache markcacheDP
#define YPrestorecache restorecacheDP
#define YPsetmudim setmudimDP
#define YPgetmudim getmudim
#define YPsetdelta setdeltaDP
#define YPgetdelta getdelta
#define YPsetlambda setlambdaDP
#define YPgetlambda getlambda
#define YPsetstab0 setstab0DP
#define YPsetstab2 setstab2DP
#define YPsetstab3 setstab3DP
#define YPsetstab4 setstab4DP
#define YPsetstab5 setstab5DP
#define YPsetstab6 setstab6DP
#define YPsetdebugkey setdebugkeyDP
#define YPsetversionkey setversionkeyDP
#define YPsetmaxdev setmaxdevDP
#define YPsetrunmode setrunmodeDP
#define YPNameData NameDataDP

#define CYPLOG LOGC
#define CYPeta etaC
#define CYPetax etaxC
#define CYPRoot2 Root2C
#define CYPSp SpC
#define YPLi2 Li2
#define CYPSp1 Sp1C
#define CYPLi2n Li2nC
* Ninh working
#define CYPLi2nx Li2nxC
#define CYPLambda LambdaC
#define CYPKam KamC
#define YPirk irk
#define YPDltB DltB
#define YPErrId ErrId

#define CYPDecomp DecompC
#define CYPSolve SolveC
#define CYPDet DetC
#define CYPInverse InverseC

#define CYPA0s A0sC
#define CYPA00s A00sC
#define CYPBcoeff BcoeffC
#define CYPCcoeff CcoeffC
#define CYPDcoeff DcoeffC
#define CYPEcoeff EcoeffC
#define CYPFcoeff FcoeffC
#define CYPDumpPara DumpParaC
#define CYPDumpCoeff DumpCoeffC
#define CYPBscalar BscalarC
#define CYPBtensor BtensorC
#define CYPBnx BnxC
#define CYPFnx FnxC
#define CYPdBnx dBnxC
#define CYPGnx GnxC
#define CYPHnx HnxC
#define CYPKnx KnxC
#define CYPLnx LnxC
#define CYPCscalar CscalarC
#define CYPC0ir C0irC
#define CYPCtensor CtensorC
#define CYPRotate RotateC
#define CYPRotateP5 RotateP5C
#define CYPSetM0 SetM0C
#define CYPRotateM0 RotateM0C
#define CYPRotate24 Rotate24C
#define CYPD0_abc D0_abcC
#define CYPDscalar DscalarC
#define CYPD0per D0perC
#define CYPD0ir_3m0q_2p0a D0ir_3m0q_2p0aC
#define CYPD0ir_3m0q_2p0b D0ir_3m0q_2p0bC
#define CYPD0ir_2m0q_0p0a D0ir_2m0q_0p0aC
#define CYPD0ir_2m0q_0p0b D0ir_2m0q_0p0bC
#define CYPD0ir_1m0q_2p0a D0ir_1m0q_2p0aC
#define CYPD0ir_0m0q_4p0a D0ir_0m0q_4p0aC
#define CYPD0ir_1m0q_0p0a D0ir_1m0q_0p0aC
#define CYPD0ir_1m0q_0p0b D0ir_1m0q_0p0bC
#define CYPD0ir_1m0q_0p0c D0ir_1m0q_0p0cC

#define CYPD0ir_4m0_4p0a D0ir_4m0_4p0aC
#define CYPD0ir_4m0_3p0a D0ir_4m0_3p0aC
#define CYPD0ir_4m0_2p0a D0ir_4m0_2p0aC
#define CYPD0ir_3m0_2p0a D0ir_3m0_2p0aC
#define CYPD0ir_3m0_2p0b D0ir_3m0_2p0bC
#define CYPD0ir_3m0_2p0c D0ir_3m0_2p0cC
#define CYPD0ir_2m0_1p0a D0ir_2m0_1p0aC
#define CYPD0ir_2m0_1p0b D0ir_2m0_1p0bC
#define CYPD0ir_2m0_1p0c D0ir_2m0_1p0cC
#define CYPD0ir_2m0_0p0a D0ir_2m0_0p0aC
#define CYPD0ir_1m0_0p0a D0ir_1m0_0p0aC
#define CYPD0M0 D0M0C
#define CYPD0DD D0DDC
#define CYPDtensor DtensorC
#define CYPEscalar EscalarC
#define CYPEtensor EtensorC
* Ninh checking
#define CYPEscalar1 Escalar1C
#define CYPEtensor1 Etensor1C
#define CYPEscalar2 Escalar2C
#define CYPEtensor2 Etensor2C

#define CYPFscalar FscalarC
#define CYPFtensor FtensorC

#define CYPfT fTC
#define CYPtriangle triangleC
#define CYPfT_lin fT_linC
#define CYPfT13 fT13C
#define CYPfS3n fS3nC
#define CYPfS3n_lin fS3n_linC
#define CYPfS2 fS2C
#define CYPfS2_lin fS2_linC
#define CYPfS2_linr fS2_linrC
#define CYPfRn fRnC
#define CYPfR fRC

#endif

#define RC 2
#ifdef COMPLEXPARA
#define XVAR CVAR
#else
#define XVAR RVAR
#endif


* for double precision
#ifndef QUADPREC

* for complex masses
#ifdef COMPLEXPARA

#define XYPA0 A0C
#define XYPA00 A00C

#define XYPBget BgetC
#define XYPB0i B0iC
#define XYPCget CgetC
#define XYPDget DgetC
#define XYPEget EgetC
#define XYPFget FgetC

* for real masses
#else

#define XYPA0 A0
#define XYPA00 A00

#define XYPBget Bget
#define XYPB0i B0i
#define XYPCget Cget
#define XYPDget Dget
#define XYPEget Eget
#define XYPFget Fget

#endif

* for quadruple precision
#else

* for complex masses
#ifdef COMPLEXPARA

#define XYPA0 A0CQP
#define XYPA00 A00CQP

#define XYPBget BgetCQP
#define XYPB0i B0iCQP
#define XYPCget CgetCQP
#define XYPDget DgetCQP
#define XYPEget EgetCQP
#define XYPFget FgetCQP

* for real masses
#else

#define XYPA0 A0QP
#define XYPA00 A00QP

#define XYPBget BgetQP
#define XYPB0i B0iQP
#define XYPCget CgetQP
#define XYPDget DgetQP
#define XYPEget EgetQP
#define XYPFget FgetQP

* for complex masses
#endif

* for double precision
#endif

***[versions
#define KeyC0 4
#define KeyD0 6
#define KeyE0 8
***]versions

***[DEBUG
#define DebugB 0
#define DebugC 1
#define DebugD 2
#define DebugE 3
#define DebugF 4
***]DEBUG

***[cache
#ifndef KIND
#define KIND 1
#endif

#ifndef QUADPREC
#define BvalYP(id,p) cache(p+id,RC)
#define CvalYP(id,p) cache(p+id,RC+2)
#define DvalYP(id,p) cache(p+id,RC+4)
#define EvalYP(id,p) cache(p+id,RC+6)
#define FvalYP(id,p) cache(p+id,RC+8)

#else
#define BvalYP(id,p) cacheQP(p+id,RC)
#define CvalYP(id,p) cacheQP(p+id,RC+2)
#define DvalYP(id,p) cacheQP(p+id,RC+4)
#define EvalYP(id,p) cacheQP(p+id,RC+6)
#define FvalYP(id,p) cacheQP(p+id,RC+8)
#endif
***]cache

***[number of tensor ranks and input parameters
#define bb0 1
#define bb1 2
#define bb00 3
#define bb11 4
#define bb001 5
#define bb111 6
#define dbb0 7
#define dbb1 8
#define dbb00 9
#define dbb11 10
#define bb0000 11
#define bb0011 12
#define bb1111 13
#define bb00001 14
#define bb00111 15
#define bb11111 16
#define Pbb 3
#define Nbb 16

#define cc0 1
#define cc1 2
#define cc2 3
#define cc00 4
#define cc11 5
#define cc12 6
#define cc22 7
#define cc001 8
#define cc002 9
#define cc111 10
#define cc112 11
#define cc122 12
#define cc222 13
#define cc0000 14
#define cc0011 15
#define cc0012 16
#define cc0022 17
#define cc1111 18
#define cc1112 19
#define cc1122 20
#define cc1222 21
#define cc2222 22
#define cc00001 23
#define cc00002 24
#define cc00111 25
#define cc00112 26
#define cc00122 27
#define cc00222 28
#define cc11111 29
#define cc11112 30
#define cc11122 31
#define cc11222 32
#define cc12222 33
#define cc22222 34
#define Pcc 6
#define Ncc 34

#define dd0 1
#define dd1 2
#define dd2 3
#define dd3 4
#define dd00 5
#define dd11 6
#define dd12 7
#define dd13 8
#define dd22 9
#define dd23 10
#define dd33 11
#define dd001 12
#define dd002 13
#define dd003 14
#define dd111 15
#define dd112 16
#define dd113 17
#define dd122 18
#define dd123 19
#define dd133 20
#define dd222 21
#define dd223 22
#define dd233 23
#define dd333 24
#define dd0000 25
#define dd0011 26
#define dd0012 27
#define dd0013 28
#define dd0022 29
#define dd0023 30
#define dd0033 31
#define dd1111 32
#define dd1112 33
#define dd1113 34
#define dd1122 35
#define dd1123 36
#define dd1133 37
#define dd1222 38
#define dd1223 39
#define dd1233 40
#define dd1333 41
#define dd2222 42
#define dd2223 43
#define dd2233 44
#define dd2333 45
#define dd3333 46
#define dd00001 47
#define dd00002 48
#define dd00003 49
#define dd00111 50
#define dd00112 51
#define dd00113 52
#define dd00122 53
#define dd00123 54
#define dd00133 55
#define dd00222 56
#define dd00223 57
#define dd00233 58
#define dd00333 59
#define dd11111 60
#define dd11112 61
#define dd11113 62
#define dd11122 63
#define dd11123 64
#define dd11133 65
#define dd11222 66
#define dd11223 67
#define dd11233 68
#define dd11333 69
#define dd12222 70
#define dd12223 71
#define dd12233 72
#define dd12333 73
#define dd13333 74
#define dd22222 75
#define dd22223 76
#define dd22233 77
#define dd22333 78
#define dd23333 79
#define dd33333 80
#define Pdd 10
#define Ndd 80

#define ee0 1
#define ee1 2
#define ee2 3
#define ee3 4
#define ee4 5
#define ee00 6
#define ee11 7
#define ee12 8
#define ee13 9
#define ee14 10
#define ee22 11
#define ee23 12
#define ee24 13
#define ee33 14
#define ee34 15
#define ee44 16
#define ee001 17
#define ee002 18
#define ee003 19
#define ee004 20
#define ee111 21
#define ee112 22
#define ee113 23
#define ee114 24
#define ee122 25
#define ee123 26
#define ee124 27
#define ee133 28
#define ee134 29
#define ee144 30
#define ee222 31
#define ee223 32
#define ee224 33
#define ee233 34
#define ee234 35
#define ee244 36
#define ee333 37
#define ee334 38
#define ee344 39
#define ee444 40
#define ee0000 41
#define ee0011 42
#define ee0012 43
#define ee0013 44
#define ee0014 45
#define ee0022 46
#define ee0023 47
#define ee0024 48
#define ee0033 49
#define ee0034 50
#define ee0044 51
#define ee1111 52
#define ee1112 53
#define ee1113 54
#define ee1114 55
#define ee1122 56
#define ee1123 57
#define ee1124 58
#define ee1133 59
#define ee1134 60
#define ee1144 61
#define ee1222 62
#define ee1223 63
#define ee1224 64
#define ee1233 65
#define ee1234 66
#define ee1244 67
#define ee1333 68
#define ee1334 69
#define ee1344 70
#define ee1444 71
#define ee2222 72
#define ee2223 73
#define ee2224 74
#define ee2233 75
#define ee2234 76
#define ee2244 77
#define ee2333 78
#define ee2334 79
#define ee2344 80
#define ee2444 81
#define ee3333 82
#define ee3334 83
#define ee3344 84
#define ee3444 85
#define ee4444 86
#define Pee 15
#define Nee 86

#define ff0 1
#define ff1 2
#define ff2 3
#define ff3 4
#define ff4 5
#define ff5 6
#define ff00 7
#define ff11 8
#define ff12 9
#define ff13 10
#define ff14 11
#define ff15 12
#define ff22 13
#define ff23 14
#define ff24 15
#define ff25 16
#define ff33 17
#define ff34 18
#define ff35 19
#define ff44 20
#define ff45 21
#define ff55 22
#define ff001 23
#define ff002 24
#define ff003 25
#define ff004 26
#define ff005 27
#define ff111 28
#define ff112 29
#define ff113 30
#define ff114 31
#define ff115 32
#define ff122 33
#define ff123 34
#define ff124 35
#define ff125 36
#define ff133 37
#define ff134 38
#define ff135 39
#define ff144 40
#define ff145 41
#define ff155 42
#define ff222 43
#define ff223 44
#define ff224 45
#define ff225 46
#define ff233 47
#define ff234 48
#define ff235 49
#define ff244 50
#define ff245 51
#define ff255 52
#define ff333 53
#define ff334 54
#define ff335 55
#define ff344 56
#define ff345 57
#define ff355 58
#define ff444 59
#define ff445 60
#define ff455 61
#define ff555 62
#define ff0000 63
#define ff0011 64
#define ff0012 65
#define ff0013 66
#define ff0014 67
#define ff0015 68
#define ff0022 69
#define ff0023 70
#define ff0024 71
#define ff0025 72
#define ff0033 73
#define ff0034 74
#define ff0035 75
#define ff0044 76
#define ff0045 77
#define ff0055 78
#define ff1111 79
#define ff1112 80
#define ff1113 81
#define ff1114 82
#define ff1115 83
#define ff1122 84
#define ff1123 85
#define ff1124 86
#define ff1125 87
#define ff1133 88
#define ff1134 89
#define ff1135 90
#define ff1144 91
#define ff1145 92
#define ff1155 93
#define ff1222 94
#define ff1223 95
#define ff1224 96
#define ff1225 97
#define ff1233 98
#define ff1234 99
#define ff1235 100
#define ff1244 101
#define ff1245 102
#define ff1255 103
#define ff1333 104
#define ff1334 105
#define ff1335 106
#define ff1344 107
#define ff1345 108
#define ff1355 109
#define ff1444 110
#define ff1445 111
#define ff1455 112
#define ff1555 113
#define ff2222 114
#define ff2223 115
#define ff2224 116
#define ff2225 117
#define ff2233 118
#define ff2234 119
#define ff2235 120
#define ff2244 121
#define ff2245 122
#define ff2255 123
#define ff2333 124
#define ff2334 125
#define ff2335 126
#define ff2344 127
#define ff2345 128
#define ff2355 129
#define ff2444 130
#define ff2445 131
#define ff2455 132
#define ff2555 133
#define ff3333 134
#define ff3334 135
#define ff3335 136
#define ff3344 137
#define ff3345 138
#define ff3355 139
#define ff3444 140
#define ff3445 141
#define ff3455 142
#define ff3555 143
#define ff4444 144
#define ff4445 145
#define ff4455 146
#define ff4555 147
#define ff5555 148
#define ff00001 149
#define ff00002 150
#define ff00003 151
#define ff00004 152
#define ff00005 153
#define ff00111 154
#define ff00112 155
#define ff00122 156
#define ff00222 157
#define ff00113 158
#define ff00123 159
#define ff00223 160
#define ff00133 161
#define ff00233 162
#define ff00333 163
#define ff00114 164
#define ff00124 165
#define ff00224 166
#define ff00134 167
#define ff00234 168
#define ff00334 169
#define ff00144 170
#define ff00244 171
#define ff00344 172
#define ff00444 173
#define ff00115 174
#define ff00125 175
#define ff00225 176
#define ff00135 177
#define ff00235 178
#define ff00335 179
#define ff00145 180
#define ff00245 181
#define ff00345 182
#define ff00445 183
#define ff00155 184
#define ff00255 185
#define ff00355 186
#define ff00455 187
#define ff00555 188
#define ff11111 189
#define ff11112 190
#define ff11122 191
#define ff11222 192
#define ff12222 193
#define ff22222 194
#define ff11113 195
#define ff11123 196
#define ff11223 197
#define ff12223 198
#define ff22223 199
#define ff11133 200
#define ff11233 201
#define ff12233 202
#define ff22233 203
#define ff11333 204
#define ff12333 205
#define ff22333 206
#define ff13333 207
#define ff23333 208
#define ff33333 209
#define ff11114 210
#define ff11124 211
#define ff11224 212
#define ff12224 213
#define ff22224 214
#define ff11134 215
#define ff11234 216
#define ff12234 217
#define ff22234 218
#define ff11334 219
#define ff12334 220
#define ff22334 221
#define ff13334 222
#define ff23334 223
#define ff33334 224
#define ff11144 225
#define ff11244 226
#define ff12244 227
#define ff22244 228
#define ff11344 229
#define ff12344 230
#define ff22344 231
#define ff13344 232
#define ff23344 233
#define ff33344 234
#define ff11444 235
#define ff12444 236
#define ff22444 237
#define ff13444 238
#define ff23444 239
#define ff33444 240
#define ff14444 241
#define ff24444 242
#define ff34444 243
#define ff44444 244
#define ff11115 245
#define ff11125 246
#define ff11225 247
#define ff12225 248
#define ff22225 249
#define ff11135 250
#define ff11235 251
#define ff12235 252
#define ff22235 253
#define ff11335 254
#define ff12335 255
#define ff22335 256
#define ff13335 257
#define ff23335 258
#define ff33335 259
#define ff11145 260
#define ff11245 261
#define ff12245 262
#define ff22245 263
#define ff11345 264
#define ff12345 265
#define ff22345 266
#define ff13345 267
#define ff23345 268
#define ff33345 269
#define ff11445 270
#define ff12445 271
#define ff22445 272
#define ff13445 273
#define ff23445 274
#define ff33445 275
#define ff14445 276
#define ff24445 277
#define ff34445 278
#define ff44445 279
#define ff11155 280
#define ff11255 281
#define ff12255 282
#define ff22255 283
#define ff11355 284
#define ff12355 285
#define ff22355 286
#define ff13355 287
#define ff23355 288
#define ff33355 289
#define ff11455 290
#define ff12455 291
#define ff22455 292
#define ff13455 293
#define ff23455 294
#define ff33455 295
#define ff14455 296
#define ff24455 297
#define ff34455 298
#define ff44455 299
#define ff11555 300
#define ff12555 301
#define ff22555 302
#define ff13555 303
#define ff23555 304
#define ff33555 305
#define ff14555 306
#define ff24555 307
#define ff34555 308
#define ff44555 309
#define ff15555 310
#define ff25555 311
#define ff35555 312
#define ff45555 313
#define ff55555 314
#define Pff 21
#define Nff 314
***]number of tensor ranks and input parameters
