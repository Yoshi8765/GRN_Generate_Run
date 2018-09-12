# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:42:03 2018

@author: Kateka
"""

import tellurium as te
import roadrunner

from GetModel import get_model

antimony_str = get_model(5)
r = te.loada(antimony_str)

#s2 = get_model(-1)

#s3 = get_model(5, reg_probs=["hi"])

#s4 = get_model(50, reg_probs=[0.5, 0.5, 0, 0, 0])
#print s4

s5 = get_model(5)

print ("\n\n\ndone!")



