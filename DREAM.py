# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 15:00:49 2018

@author: Yoshi
"""

import os
import tellurium as te

r=te.loada('''
//Created by libAntimony v2.7.0
model *subchallenge2_missing_connections()

  // Compartments and Species:
  compartment compartment_;
  species $X in compartment_, $__src__ in compartment_, $__waste__ in compartment_;
  species p1 in compartment_, p2 in compartment_, p3 in compartment_, p4 in compartment_;
  species p5 in compartment_, p6 in compartment_, p11 in compartment_, p7 in compartment_;
  species p8 in compartment_, p9 in compartment_, p10 in compartment_;

  // Reactions:
  _J0: $X => p1; (p1_synthesis_rate*(X/r1_Kd)^r1_h)/(1 + (X/r1_Kd)^r1_h);
  _J1: p1 => $__waste__; p1_degradation_rate*p1;
  _J2: $__src__ => p2; (p2_synthesis_rate*(p1/r2_Kd)^r2_h)/(1 + (p1/r2_Kd)^r2_h)/(1 + (p6/r5_Kd)^r5_h);
  _J3: p2 => $__waste__; p2_degradation_rate*p2;
  _J4: $__src__ => p3; (p3_synthesis_rate_1*(p2/r4_Kd)^r4_h)/(1 + (p2/r4_Kd)^r4_h) + (p3_synthesis_rate_2*(p1/r3_Kd)^r3_h)/(1 + (p1/r3_Kd)^r3_h);
  _J5: p3 => $__waste__; p3_degradation_rate*p3;
  _J6: p3 => p4; p4_synthesis_rate_1*p3;
  _J7: $__src__ => p4; p4_synthesis_rate_2/(1 + (p7/r8_Kd)^r8_h);
  _J8: p4 => $__waste__; p4_degradation_rate*p4;
  _J9: $__src__ => p5; r6_basal + (p5_synthesis_rate_1*(p4/r6_Kd)^r6_h)/(1 + (p4/r6_Kd)^r6_h) + (p5_synthesis_rate_3*(p5/r7_Kd)^r7_h)/(1 + (p5/r7_Kd)^r7_h) + p5_synthesis_rate_2;
  _J10: p5 => $__waste__; p5_degradation_rate*p5;
  _J11: $__src__ => p6; p6_synthesis_rate;
  _J12: p6 => $__waste__; p6_degradation_rate*p6;
  _J13: $__src__ => p11; r11_basal + (p11_synthesis_rate_1*(p5/r11_Kd)^r11_h)/(1 + (p5/r11_Kd)^r11_h);
  _J14: p11 => $__waste__; p11_degradation_rate*p11;
  _J15: $__src__ => p7; p7_synthesis_rate*(r11_basal + (p11_synthesis_rate_1*(p5/r11_Kd)^r11_h)/(1 + (p5/r11_Kd)^r11_h));
  _J16: p7 => $__waste__; p7_degradation_rate*p7;
  _J17: $__src__ => p8; p8_synthesis_rate/(1 + (p7/r14_Kd)^r14_h);
  _J18: p8 => $__waste__; p8_degradation_rate*p8;
  _J19: $__src__ => p9; (p9_synthesis_rate*(p8/r15_Kd)^r15_h)/(1 + (p8/r15_Kd)^r15_h);
  _J20: p9 => $__waste__; p9_degradation_rate*p9;
  _J21: $__src__ => p10; p10_synthesis_rate/(1 + (p9/r16_Kd)^r16_h)/(1 + (p7/r13_Kd)^r13_h);
  _J22: p10 => $__waste__; p10_degradation_rate*p10;

  // Species initializations:
  X = 1;
  __src__ = 1;
  __waste__ = 1;
  p1 = 1;
  p2 = 1;
  p3 = 1;
  p4 = 1;
  p5 = 1;
  p6 = 1;
  p11 = 1;
  p7 = 1;
  p8 = 1;
  p9 = 1;
  p10 = 1;

  // Compartment initializations:
  compartment_ = 1;

  // Variable initializations:
  p1_synthesis_rate = 1;
  r1_Kd = 1;
  r1_h = 1;
  p1_degradation_rate = 1;
  p2_synthesis_rate = 1;
  r2_Kd = 1;
  r2_h = 1;
  r5_Kd = 1;
  r5_h = 1;
  p2_degradation_rate = 1;
  p3_synthesis_rate_1 = 1;
  r4_Kd = 1;
  r4_h = 1;
  p3_synthesis_rate_2 = 1;
  r3_Kd = 1;
  r3_h = 1;
  p3_degradation_rate = 1;
  p4_synthesis_rate_1 = 1;
  p4_synthesis_rate_2 = 1;
  r8_Kd = 1;
  r8_h = 1;
  p4_degradation_rate = 1;
  r6_basal = 1;
  p5_synthesis_rate_1 = 1;
  r6_Kd = 1;
  r6_h = 1;
  p5_synthesis_rate_3 = 1;
  r7_Kd = 1;
  r7_h = 1;
  p5_synthesis_rate_2 = 1;
  p5_degradation_rate = 1;
  p6_synthesis_rate = 1;
  p6_degradation_rate = 1;
  r11_basal = 1;
  p11_synthesis_rate_1 = 1;
  r11_Kd = 1;
  r11_h = 1;
  p11_degradation_rate = 1;
  p7_synthesis_rate = 1;
  p7_degradation_rate = 1;
  p8_synthesis_rate = 1;
  r14_Kd = 1;
  r14_h = 1;
  p8_degradation_rate = 1;
  p9_synthesis_rate = 1;
  r15_Kd = 1;
  r15_h = 1;
  p9_degradation_rate = 1;
  p10_synthesis_rate = 1;
  r16_Kd = 1;
  r16_h = 1;
  r13_Kd = 1;
  r13_h = 1;
  p10_degradation_rate = 1;

  // Other declarations:
  const compartment_, p1_synthesis_rate, r1_Kd, r1_h, p1_degradation_rate;
  const p2_synthesis_rate, r2_Kd, r2_h, r5_Kd, r5_h, p2_degradation_rate;
  const p3_synthesis_rate_1, r4_Kd, r4_h, p3_synthesis_rate_2, r3_Kd, r3_h;
  const p3_degradation_rate, p4_synthesis_rate_1, p4_synthesis_rate_2, r8_Kd;
  const r8_h, p4_degradation_rate, r6_basal, p5_synthesis_rate_1, r6_Kd, r6_h;
  const p5_synthesis_rate_3, r7_Kd, r7_h, p5_synthesis_rate_2, p5_degradation_rate;
  const p6_synthesis_rate, p6_degradation_rate, r11_basal, p11_synthesis_rate_1;
  const r11_Kd, r11_h, p11_degradation_rate, p7_synthesis_rate, p7_degradation_rate;
  const p8_synthesis_rate, r14_Kd, r14_h, p8_degradation_rate, p9_synthesis_rate;
  const r15_Kd, r15_h, p9_degradation_rate, p10_synthesis_rate, r16_Kd, r16_h;
  const r13_Kd, r13_h, p10_degradation_rate;
end        
''')
    
a = r.simulate(0,10,100)
r.plot()