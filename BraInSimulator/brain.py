import os
import random
import sys
import numpy as np
import scipy as scpy
import scipy.integrate as spin
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

class Neuron:
    # Default constructor
    def __init__(self, tmax, ):

    #self, tmax, I_ampl, unew, vnew, a, b, tau, timescale, coupling, ):



    def defaultfhn(self):
        self.I_ampl = 0.324
        self.tmax = 100
        self.unew = -0.6
        self.vnew = -0.4
        self.a = 0.7
        self.b = 0.8
        self.tau = 12.5

    def singlefhn(self, timescale, coupling, tmax,  ):
        self.
