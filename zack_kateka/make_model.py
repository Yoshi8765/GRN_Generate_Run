# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:42:03 2018

@author: Kateka & Zachary McNulty
"""

import tellurium as te
import time
import numpy as np

from GetModel import get_model
from GetModel import convert_to_biotapestry
from Biotapestry import convert_biotapestry_to_antimony




"""
This code will generate an antimony model for an unbroken network, and output this as a .txt file. There are a few settings you can mess with to generate a suitable network, and I describe them below. When you have set them, run this script and a model will be generated in your working directory. It also plots the timecourse data for the model, so you will have a chance to see if you like the behavior it exhibits.

options:

seednum will be the seed for the random processes involved in generating the model. Some models are not as interesting as others, so you may want to play with the seednum until you get something you want.

**You may want to save the seed you use for reproducibility.**

num_genes will be the number of genes in the model

reachability defines the proportion of  genes in the network you want to be connected to the INPUT, directly or indirectly. This is really only important if the INPUT is high relative to the other protein concentrations.

self_feedback_min is the  min number of self feedback loops that the model will have. If you want self feedback, increasing this number will do the trick.

Do NOT change the model name from the default of `pathway`; the other code requires the name to be `pathway`.
"""

seednum = 12345
num_genes = 8
reachability = 0.9
self_feedback_min = 0

ant_str, biotap_str = get_model(num_genes, seed = seednum, model_name ="pathway" , reachability=reachability, self_feedback_min=self_feedback_min, export=True)

r = te.loada(ant_str)
r.simulate(0,200,1000)
r.plot()

print("Done! Model has been generated")
