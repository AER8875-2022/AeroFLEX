// Gen mesh with no gui using command:
// gmsh naca0012.geo -2 -o naca0012_coarse.msh

// Choose mesh size // number of elements
lc1=1.0; lc2=0.075; // 4 000
// lc1=0.5; lc2=0.025; // 20 000
// lc1=0.275; lc2=0.02; // 50 000
// lc1=0.2; lc2=0.0125; // 100 000
// lc1=0.12; lc2=0.01; // 250 000

// Radius of farfield. For element numbers listed above,
//  keep wc equal to 5.
wc = 5;
aoa = 0.;

Point(1)={0+0.5, 0, 0, lc1};
Point(2)={wc+0.5, 0, 0, lc1};
Point(3)={0+0.5, wc, 0, lc1};
Point(4)={-wc+0.5, 0, 0, lc1};
Point(5)={0+0.5, -wc, 0, lc1};

// NACA 0012
Point(101)={+1.0000000000000000e+000,+0.0000000000000000e+000,+0.0000000000000000e+000,lc2};
Point(102)={+9.9199000000000004e-001,-1.3209999999999999e-003,+0.0000000000000000e+000,lc2};
Point(103)={+9.8340000000000005e-001,-2.7246000000000002e-003,+0.0000000000000000e+000,lc2};
Point(104)={+9.7423999999999999e-001,-4.2063999999999999e-003,+0.0000000000000000e+000,lc2};
Point(105)={+9.6453999999999995e-001,-5.7619000000000004e-003,+0.0000000000000000e+000,lc2};
Point(106)={+9.5430999999999999e-001,-7.3866000000000001e-003,+0.0000000000000000e+000,lc2};
Point(107)={+9.4355999999999995e-001,-8.9133000000000007e-003,+0.0000000000000000e+000,lc2};
Point(108)={+9.3230000000000002e-001,-1.0378999999999999e-002,+0.0000000000000000e+000,lc2};
Point(109)={+9.2056000000000004e-001,-1.1886000000000001e-002,+0.0000000000000000e+000,lc2};
Point(110)={+9.0834999999999999e-001,-1.3433000000000000e-002,+0.0000000000000000e+000,lc2};
Point(111)={+8.9568999999999999e-001,-1.5013000000000000e-002,+0.0000000000000000e+000,lc2};
Point(112)={+8.8258999999999999e-001,-1.6622999999999999e-002,+0.0000000000000000e+000,lc2};
Point(113)={+8.6907999999999996e-001,-1.8258000000000000e-002,+0.0000000000000000e+000,lc2};
Point(114)={+8.5518000000000005e-001,-1.9914999999999999e-002,+0.0000000000000000e+000,lc2};
Point(115)={+8.4091000000000005e-001,-2.1588000000000000e-002,+0.0000000000000000e+000,lc2};
Point(116)={+8.2628000000000001e-001,-2.3274000000000000e-002,+0.0000000000000000e+000,lc2};
Point(117)={+8.1132000000000004e-001,-2.4969000000000002e-002,+0.0000000000000000e+000,lc2};
Point(118)={+7.9605000000000004e-001,-2.6668000000000001e-002,+0.0000000000000000e+000,lc2};
Point(119)={+7.8047999999999995e-001,-2.8368000000000001e-002,+0.0000000000000000e+000,lc2};
Point(120)={+7.6465000000000005e-001,-3.0062999999999999e-002,+0.0000000000000000e+000,lc2};
Point(121)={+7.4858000000000002e-001,-3.1751000000000001e-002,+0.0000000000000000e+000,lc2};
Point(122)={+7.3228000000000004e-001,-3.3426999999999998e-002,+0.0000000000000000e+000,lc2};
Point(123)={+7.1577999999999997e-001,-3.5085999999999999e-002,+0.0000000000000000e+000,lc2};
Point(124)={+6.9910000000000005e-001,-3.6726000000000002e-002,+0.0000000000000000e+000,lc2};
Point(125)={+6.8227000000000004e-001,-3.8342000000000001e-002,+0.0000000000000000e+000,lc2};
Point(126)={+6.6530000000000000e-001,-3.9928999999999999e-002,+0.0000000000000000e+000,lc2};
Point(127)={+6.4822000000000002e-001,-4.1485000000000001e-002,+0.0000000000000000e+000,lc2};
Point(128)={+6.3105999999999995e-001,-4.3004000000000001e-002,+0.0000000000000000e+000,lc2};
Point(129)={+6.1382000000000003e-001,-4.4482000000000001e-002,+0.0000000000000000e+000,lc2};
Point(130)={+5.9655000000000002e-001,-4.5915999999999998e-002,+0.0000000000000000e+000,lc2};
Point(131)={+5.7925000000000004e-001,-4.7301999999999997e-002,+0.0000000000000000e+000,lc2};
Point(132)={+5.6194999999999995e-001,-4.8634999999999998e-002,+0.0000000000000000e+000,lc2};
Point(133)={+5.4468000000000005e-001,-4.9911999999999998e-002,+0.0000000000000000e+000,lc2};
Point(134)={+5.2744000000000002e-001,-5.1128000000000000e-002,+0.0000000000000000e+000,lc2};
Point(135)={+5.1027000000000000e-001,-5.2280000000000000e-002,+0.0000000000000000e+000,lc2};
Point(136)={+4.9319000000000002e-001,-5.3365000000000003e-002,+0.0000000000000000e+000,lc2};
Point(137)={+4.7621000000000002e-001,-5.4378000000000003e-002,+0.0000000000000000e+000,lc2};
Point(138)={+4.5934999999999998e-001,-5.5316999999999998e-002,+0.0000000000000000e+000,lc2};
Point(139)={+4.4263999999999998e-001,-5.6177999999999999e-002,+0.0000000000000000e+000,lc2};
Point(140)={+4.2609000000000002e-001,-5.6958000000000002e-002,+0.0000000000000000e+000,lc2};
Point(141)={+4.0971999999999997e-001,-5.7653999999999997e-002,+0.0000000000000000e+000,lc2};
Point(142)={+3.9354000000000000e-001,-5.8264000000000003e-002,+0.0000000000000000e+000,lc2};
Point(143)={+3.7758000000000003e-001,-5.8785999999999998e-002,+0.0000000000000000e+000,lc2};
Point(144)={+3.6185000000000000e-001,-5.9218000000000000e-002,+0.0000000000000000e+000,lc2};
Point(145)={+3.4637000000000001e-001,-5.9558000000000000e-002,+0.0000000000000000e+000,lc2};
Point(146)={+3.3115000000000000e-001,-5.9804999999999997e-002,+0.0000000000000000e+000,lc2};
Point(147)={+3.1619999999999998e-001,-5.9957999999999997e-002,+0.0000000000000000e+000,lc2};
Point(148)={+3.0153000000000002e-001,-6.0017000000000001e-002,+0.0000000000000000e+000,lc2};
Point(149)={+2.8717999999999999e-001,-5.9979999999999999e-002,+0.0000000000000000e+000,lc2};
Point(150)={+2.7312999999999998e-001,-5.9850000000000000e-002,+0.0000000000000000e+000,lc2};
Point(151)={+2.5940000000000002e-001,-5.9624999999999997e-002,+0.0000000000000000e+000,lc2};
Point(152)={+2.4601000000000001e-001,-5.9306999999999999e-002,+0.0000000000000000e+000,lc2};
Point(153)={+2.3296000000000000e-001,-5.8896999999999998e-002,+0.0000000000000000e+000,lc2};
Point(154)={+2.2026000000000001e-001,-5.8397999999999999e-002,+0.0000000000000000e+000,lc2};
Point(155)={+2.0791999999999999e-001,-5.7810000000000000e-002,+0.0000000000000000e+000,lc2};
Point(156)={+1.9595000000000001e-001,-5.7135999999999999e-002,+0.0000000000000000e+000,lc2};
Point(157)={+1.8435000000000001e-001,-5.6377999999999998e-002,+0.0000000000000000e+000,lc2};
Point(158)={+1.7313000000000001e-001,-5.5539999999999999e-002,+0.0000000000000000e+000,lc2};
Point(159)={+1.6228999999999999e-001,-5.4625000000000000e-002,+0.0000000000000000e+000,lc2};
Point(160)={+1.5182999999999999e-001,-5.3636000000000003e-002,+0.0000000000000000e+000,lc2};
Point(161)={+1.4176000000000000e-001,-5.2575999999999998e-002,+0.0000000000000000e+000,lc2};
Point(162)={+1.3208000000000000e-001,-5.1450000000000003e-002,+0.0000000000000000e+000,lc2};
Point(163)={+1.2279000000000000e-001,-5.0261000000000000e-002,+0.0000000000000000e+000,lc2};
Point(164)={+1.1389000000000001e-001,-4.9013000000000001e-002,+0.0000000000000000e+000,lc2};
Point(165)={+1.0538000000000000e-001,-4.7711000000000003e-002,+0.0000000000000000e+000,lc2};
Point(166)={+9.7255999999999995e-002,-4.6358000000000003e-002,+0.0000000000000000e+000,lc2};
Point(167)={+8.9519000000000001e-002,-4.4958999999999999e-002,+0.0000000000000000e+000,lc2};
Point(168)={+8.2165000000000002e-002,-4.3519000000000002e-002,+0.0000000000000000e+000,lc2};
Point(169)={+7.5189000000000006e-002,-4.2041000000000002e-002,+0.0000000000000000e+000,lc2};
Point(170)={+6.8587999999999996e-002,-4.0529999999999997e-002,+0.0000000000000000e+000,lc2};
Point(171)={+6.2356000000000002e-002,-3.8989999999999997e-002,+0.0000000000000000e+000,lc2};
Point(172)={+5.6487000000000002e-002,-3.7425000000000000e-002,+0.0000000000000000e+000,lc2};
Point(173)={+5.0973999999999998e-002,-3.5839999999999997e-002,+0.0000000000000000e+000,lc2};
Point(174)={+4.5809999999999997e-002,-3.4237999999999998e-002,+0.0000000000000000e+000,lc2};
Point(175)={+4.0987999999999997e-002,-3.2624000000000000e-002,+0.0000000000000000e+000,lc2};
Point(176)={+3.6498999999999997e-002,-3.1001000000000001e-002,+0.0000000000000000e+000,lc2};
Point(177)={+3.2333000000000001e-002,-2.9373000000000000e-002,+0.0000000000000000e+000,lc2};
Point(178)={+2.8482000000000000e-002,-2.7744000000000001e-002,+0.0000000000000000e+000,lc2};
Point(179)={+2.4936000000000000e-002,-2.6116000000000000e-002,+0.0000000000000000e+000,lc2};
Point(180)={+2.1683000000000001e-002,-2.4494999999999999e-002,+0.0000000000000000e+000,lc2};
Point(181)={+1.8714000000000001e-002,-2.2882000000000000e-002,+0.0000000000000000e+000,lc2};
Point(182)={+1.6017000000000000e-002,-2.1281000000000001e-002,+0.0000000000000000e+000,lc2};
Point(183)={+1.3580000000000000e-002,-1.9694000000000000e-002,+0.0000000000000000e+000,lc2};
Point(184)={+1.1391999999999999e-002,-1.8124999999999999e-002,+0.0000000000000000e+000,lc2};
Point(185)={+9.4412000000000003e-003,-1.6577000000000001e-002,+0.0000000000000000e+000,lc2};
Point(186)={+7.7146000000000003e-003,-1.5051000000000000e-002,+0.0000000000000000e+000,lc2};
Point(187)={+6.1999999999999998e-003,-1.3550000000000000e-002,+0.0000000000000000e+000,lc2};
Point(188)={+4.8849999999999996e-003,-1.2076000000000000e-002,+0.0000000000000000e+000,lc2};
Point(189)={+3.7567999999999998e-003,-1.0632000000000001e-002,+0.0000000000000000e+000,lc2};
Point(190)={+2.8031000000000002e-003,-9.2178999999999994e-003,+0.0000000000000000e+000,lc2};
Point(191)={+2.0111000000000000e-003,-7.8358999999999998e-003,+0.0000000000000000e+000,lc2};
Point(192)={+1.3686000000000000e-003,-6.4863999999999998e-003,+0.0000000000000000e+000,lc2};
Point(193)={+8.6341000000000000e-004,-5.1690000000000000e-003,+0.0000000000000000e+000,lc2};
Point(194)={+4.8376000000000002e-004,-3.8815000000000000e-003,+0.0000000000000000e+000,lc2};
Point(195)={+2.1870000000000000e-004,-2.6178999999999998e-003,+0.0000000000000000e+000,lc2};
Point(196)={+5.8739999999999997e-005,-1.3609000000000000e-003,+0.0000000000000000e+000,lc2};
Point(197)={+0.0000000000000000e+000,+0.0000000000000000e+000,+0.0000000000000000e+000,lc2};
Point(198)={+5.8739999999999997e-005,+1.3609000000000000e-003,+0.0000000000000000e+000,lc2};
Point(199)={+2.1870000000000000e-004,+2.6178999999999998e-003,+0.0000000000000000e+000,lc2};
Point(200)={+4.8376000000000002e-004,+3.8815000000000000e-003,+0.0000000000000000e+000,lc2};
Point(201)={+8.6341000000000000e-004,+5.1690000000000000e-003,+0.0000000000000000e+000,lc2};
Point(202)={+1.3686000000000000e-003,+6.4863999999999998e-003,+0.0000000000000000e+000,lc2};
Point(203)={+2.0111000000000000e-003,+7.8358999999999998e-003,+0.0000000000000000e+000,lc2};
Point(204)={+2.8031000000000002e-003,+9.2178999999999994e-003,+0.0000000000000000e+000,lc2};
Point(205)={+3.7567999999999998e-003,+1.0632000000000001e-002,+0.0000000000000000e+000,lc2};
Point(206)={+4.8849999999999996e-003,+1.2076000000000000e-002,+0.0000000000000000e+000,lc2};
Point(207)={+6.1999999999999998e-003,+1.3550000000000000e-002,+0.0000000000000000e+000,lc2};
Point(208)={+7.7146000000000003e-003,+1.5051000000000000e-002,+0.0000000000000000e+000,lc2};
Point(209)={+9.4412000000000003e-003,+1.6577000000000001e-002,+0.0000000000000000e+000,lc2};
Point(210)={+1.1391999999999999e-002,+1.8124999999999999e-002,+0.0000000000000000e+000,lc2};
Point(211)={+1.3580000000000000e-002,+1.9694000000000000e-002,+0.0000000000000000e+000,lc2};
Point(212)={+1.6017000000000000e-002,+2.1281000000000001e-002,+0.0000000000000000e+000,lc2};
Point(213)={+1.8714000000000001e-002,+2.2882000000000000e-002,+0.0000000000000000e+000,lc2};
Point(214)={+2.1683000000000001e-002,+2.4494999999999999e-002,+0.0000000000000000e+000,lc2};
Point(215)={+2.4936000000000000e-002,+2.6116000000000000e-002,+0.0000000000000000e+000,lc2};
Point(216)={+2.8482000000000000e-002,+2.7744000000000001e-002,+0.0000000000000000e+000,lc2};
Point(217)={+3.2333000000000001e-002,+2.9373000000000000e-002,+0.0000000000000000e+000,lc2};
Point(218)={+3.6498999999999997e-002,+3.1001000000000001e-002,+0.0000000000000000e+000,lc2};
Point(219)={+4.0987999999999997e-002,+3.2624000000000000e-002,+0.0000000000000000e+000,lc2};
Point(220)={+4.5809999999999997e-002,+3.4237999999999998e-002,+0.0000000000000000e+000,lc2};
Point(221)={+5.0973999999999998e-002,+3.5839999999999997e-002,+0.0000000000000000e+000,lc2};
Point(222)={+5.6487000000000002e-002,+3.7425000000000000e-002,+0.0000000000000000e+000,lc2};
Point(223)={+6.2356000000000002e-002,+3.8989999999999997e-002,+0.0000000000000000e+000,lc2};
Point(224)={+6.8587999999999996e-002,+4.0529999999999997e-002,+0.0000000000000000e+000,lc2};
Point(225)={+7.5189000000000006e-002,+4.2041000000000002e-002,+0.0000000000000000e+000,lc2};
Point(226)={+8.2165000000000002e-002,+4.3519000000000002e-002,+0.0000000000000000e+000,lc2};
Point(227)={+8.9519000000000001e-002,+4.4958999999999999e-002,+0.0000000000000000e+000,lc2};
Point(228)={+9.7255999999999995e-002,+4.6358000000000003e-002,+0.0000000000000000e+000,lc2};
Point(229)={+1.0538000000000000e-001,+4.7711000000000003e-002,+0.0000000000000000e+000,lc2};
Point(230)={+1.1389000000000001e-001,+4.9013000000000001e-002,+0.0000000000000000e+000,lc2};
Point(231)={+1.2279000000000000e-001,+5.0261000000000000e-002,+0.0000000000000000e+000,lc2};
Point(232)={+1.3208000000000000e-001,+5.1450000000000003e-002,+0.0000000000000000e+000,lc2};
Point(233)={+1.4176000000000000e-001,+5.2575999999999998e-002,+0.0000000000000000e+000,lc2};
Point(234)={+1.5182999999999999e-001,+5.3636000000000003e-002,+0.0000000000000000e+000,lc2};
Point(235)={+1.6228999999999999e-001,+5.4625000000000000e-002,+0.0000000000000000e+000,lc2};
Point(236)={+1.7313000000000001e-001,+5.5539999999999999e-002,+0.0000000000000000e+000,lc2};
Point(237)={+1.8435000000000001e-001,+5.6377999999999998e-002,+0.0000000000000000e+000,lc2};
Point(238)={+1.9595000000000001e-001,+5.7135999999999999e-002,+0.0000000000000000e+000,lc2};
Point(239)={+2.0791999999999999e-001,+5.7810000000000000e-002,+0.0000000000000000e+000,lc2};
Point(240)={+2.2026000000000001e-001,+5.8397999999999999e-002,+0.0000000000000000e+000,lc2};
Point(241)={+2.3296000000000000e-001,+5.8896999999999998e-002,+0.0000000000000000e+000,lc2};
Point(242)={+2.4601000000000001e-001,+5.9306999999999999e-002,+0.0000000000000000e+000,lc2};
Point(243)={+2.5940000000000002e-001,+5.9624999999999997e-002,+0.0000000000000000e+000,lc2};
Point(244)={+2.7312999999999998e-001,+5.9850000000000000e-002,+0.0000000000000000e+000,lc2};
Point(245)={+2.8717999999999999e-001,+5.9979999999999999e-002,+0.0000000000000000e+000,lc2};
Point(246)={+3.0153000000000002e-001,+6.0017000000000001e-002,+0.0000000000000000e+000,lc2};
Point(247)={+3.1619999999999998e-001,+5.9957999999999997e-002,+0.0000000000000000e+000,lc2};
Point(248)={+3.3115000000000000e-001,+5.9804999999999997e-002,+0.0000000000000000e+000,lc2};
Point(249)={+3.4637000000000001e-001,+5.9558000000000000e-002,+0.0000000000000000e+000,lc2};
Point(250)={+3.6185000000000000e-001,+5.9218000000000000e-002,+0.0000000000000000e+000,lc2};
Point(251)={+3.7758000000000003e-001,+5.8785999999999998e-002,+0.0000000000000000e+000,lc2};
Point(252)={+3.9354000000000000e-001,+5.8264000000000003e-002,+0.0000000000000000e+000,lc2};
Point(253)={+4.0971999999999997e-001,+5.7653999999999997e-002,+0.0000000000000000e+000,lc2};
Point(254)={+4.2609000000000002e-001,+5.6958000000000002e-002,+0.0000000000000000e+000,lc2};
Point(255)={+4.4263999999999998e-001,+5.6177999999999999e-002,+0.0000000000000000e+000,lc2};
Point(256)={+4.5934999999999998e-001,+5.5316999999999998e-002,+0.0000000000000000e+000,lc2};
Point(257)={+4.7621000000000002e-001,+5.4378000000000003e-002,+0.0000000000000000e+000,lc2};
Point(258)={+4.9319000000000002e-001,+5.3365000000000003e-002,+0.0000000000000000e+000,lc2};
Point(259)={+5.1027000000000000e-001,+5.2280000000000000e-002,+0.0000000000000000e+000,lc2};
Point(260)={+5.2744000000000002e-001,+5.1128000000000000e-002,+0.0000000000000000e+000,lc2};
Point(261)={+5.4468000000000005e-001,+4.9911999999999998e-002,+0.0000000000000000e+000,lc2};
Point(262)={+5.6194999999999995e-001,+4.8634999999999998e-002,+0.0000000000000000e+000,lc2};
Point(263)={+5.7925000000000004e-001,+4.7301999999999997e-002,+0.0000000000000000e+000,lc2};
Point(264)={+5.9655000000000002e-001,+4.5915999999999998e-002,+0.0000000000000000e+000,lc2};
Point(265)={+6.1382000000000003e-001,+4.4482000000000001e-002,+0.0000000000000000e+000,lc2};
Point(266)={+6.3105999999999995e-001,+4.3004000000000001e-002,+0.0000000000000000e+000,lc2};
Point(267)={+6.4822000000000002e-001,+4.1485000000000001e-002,+0.0000000000000000e+000,lc2};
Point(268)={+6.6530000000000000e-001,+3.9928999999999999e-002,+0.0000000000000000e+000,lc2};
Point(269)={+6.8227000000000004e-001,+3.8342000000000001e-002,+0.0000000000000000e+000,lc2};
Point(270)={+6.9910000000000005e-001,+3.6726000000000002e-002,+0.0000000000000000e+000,lc2};
Point(271)={+7.1577999999999997e-001,+3.5085999999999999e-002,+0.0000000000000000e+000,lc2};
Point(272)={+7.3228000000000004e-001,+3.3426999999999998e-002,+0.0000000000000000e+000,lc2};
Point(273)={+7.4858000000000002e-001,+3.1751000000000001e-002,+0.0000000000000000e+000,lc2};
Point(274)={+7.6465000000000005e-001,+3.0062999999999999e-002,+0.0000000000000000e+000,lc2};
Point(275)={+7.8047999999999995e-001,+2.8368000000000001e-002,+0.0000000000000000e+000,lc2};
Point(276)={+7.9605000000000004e-001,+2.6668000000000001e-002,+0.0000000000000000e+000,lc2};
Point(277)={+8.1132000000000004e-001,+2.4969000000000002e-002,+0.0000000000000000e+000,lc2};
Point(278)={+8.2628000000000001e-001,+2.3274000000000000e-002,+0.0000000000000000e+000,lc2};
Point(279)={+8.4091000000000005e-001,+2.1588000000000000e-002,+0.0000000000000000e+000,lc2};
Point(280)={+8.5518000000000005e-001,+1.9914999999999999e-002,+0.0000000000000000e+000,lc2};
Point(281)={+8.6907999999999996e-001,+1.8258000000000000e-002,+0.0000000000000000e+000,lc2};
Point(282)={+8.8258999999999999e-001,+1.6622999999999999e-002,+0.0000000000000000e+000,lc2};
Point(283)={+8.9568999999999999e-001,+1.5013000000000000e-002,+0.0000000000000000e+000,lc2};
Point(284)={+9.0834999999999999e-001,+1.3433000000000000e-002,+0.0000000000000000e+000,lc2};
Point(285)={+9.2056000000000004e-001,+1.1886000000000001e-002,+0.0000000000000000e+000,lc2};
Point(286)={+9.3230000000000002e-001,+1.0378999999999999e-002,+0.0000000000000000e+000,lc2};
Point(287)={+9.4355999999999995e-001,+8.9133000000000007e-003,+0.0000000000000000e+000,lc2};
Point(288)={+9.5430999999999999e-001,+7.3866000000000001e-003,+0.0000000000000000e+000,lc2};
Point(289)={+9.6453999999999995e-001,+5.7619000000000004e-003,+0.0000000000000000e+000,lc2};
Point(290)={+9.7423999999999999e-001,+4.2063999999999999e-003,+0.0000000000000000e+000,lc2};
Point(291)={+9.8340000000000005e-001,+2.7246000000000002e-003,+0.0000000000000000e+000,lc2};
Point(292)={+9.9199000000000004e-001,+1.3209999999999999e-003,+0.0000000000000000e+000,lc2};


Circle(1)={2, 1, 3};
Circle(2)={3, 1, 4};
Circle(3)={4, 1, 5};
Circle(4)={5, 1, 2};

Spline(5)={101 : 292, 101};

Curve Loop (1)={1,2,3,4};
Curve Loop (2)={5};


Plane Surface(1)={1, 2};

Rotate{{0,0,1}, {0,1,0}, -aoa*Pi/180}{ Surface{1}; }


Physical Curve("farfield") = {1, 2, 3, 4};
Physical Curve("wall") = {5};
Physical Surface("internal") = {1};

