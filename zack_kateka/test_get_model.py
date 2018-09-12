# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:42:03 2018

@author: Kateka
"""

import tellurium as te
import roadrunner

from GetModel import get_model

ant_str = get_model(5)
r = te.loada(ant_str)

