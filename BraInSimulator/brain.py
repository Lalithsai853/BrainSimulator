import os
import random
import sys
import numpy as np
import scipy as scpy
import scipy.integrate as spin
from matplotlib import pyplot as plt
import collections
from scipy import stats
import random
import matplotlib.patches as mpatches
import pandas
import math
import cmath
import streamlit as st
import networkx as nx
import seaborn as sns

class Brain:
    # Creates the Brain object
    def __init__(self, tmax, a, u0 = None, v0 = None, I_ampl= None, b= None, tau= None, matrixA= None, matrixB= None, timescale= None, coupling= None):
        self.tmax = tmax
        self.a = a
        self.I_ampl = I_ampl
        self.u0 = u0
        self.v0 = v0
        self.b = b
        self.tau = tau
        self.matrixA = matrixA
        self.matrixB = matrixB
        self.timescale = timescale
        self.coupling = coupling
    '''
    parameters:
    inputFile: csv file containing initial connections in form (startregion, endregion)
    start: starting value for each connection
    rng: number of brain regions
    Example:
    matrixA = Brain.csvMatrixA('data.csv'm, 1, 90)
    returns:
    arrA: Matrix with connections
    '''
    @staticmethod
    def csvMatrixA(inputFile, start, rng):
        arrA = np.zeros((rng,rng))

        df = pandas.read_csv(inputFile)
        row_count, column_count = df.shape
        for i in range(1, row_count):
            arrA[df['FromNeuron'][i]][df['ToNeuron'][i]] = start

        return arrA

    '''
    parameters:
    rng: number of brain regions
    Example:
    matrixA = Brain.generatematrixA(90)
    returns:
    arrA: Matrix with random connections
    '''
    @staticmethod
    def generatematrixA(rng):
        arrA = np.zeros((rng,rng))
        for i in range(0, rng):
            num1 = randint(1, rng)
            num2 = randint(1, rng)
            arrA[num1][num2] = 1
            arrA[num2][num1] = 1
        return arrA

    '''
    parameters:
    psi: parameter to control interactions between u and v variables
    Example:
    matrixB = Brain.generatematrixB(psi)
    returns:
    arrB: 2 by 2 matrix to control coupling between u and v variables
    '''
    @staticmethod
    def generatematrixB(psi):
        arrB = np.ones((2,2))
        arrB[0][0] = math.cos(psi)
        arrB[0][1] = math.sin(psi)
        arrB[1][0] = (-1 * math.sin(psi))
        arrB[1][1] = math.cos(psi)
        return arrB
    '''
    parameters:
    u1: Corresponds to u timeseries calculations for every single brain region
    v1: Corresponds to v timeseries calculations for every single brain region
    t: Corresponds to timeseries array
    rng: Number of regions
    Example:
    uarray, varray = new_brain.calculateUVmatrix(uarray, varray, t1, rng)
    returns:
    arrA: U and V arrays containing calculations of timeseries for every single brain region
    '''
    def calculateUVmatrix(self, u1, v1, t, rng):
        uterm = 0
        vterm = 0
        for i in range(0,len(t)-1):
            for n in range(0, rng):
                dudt = u1[n][i] - np.power(u1[n][i], 3)/3 - v1[n][i]
                dvdt = u1[n][i] + self.a
                for j in range(0, rng):
                    uterm = uterm + self.matrixA[n][j] * ((self.matrixB[0][0] * (u1[j][i] - u1[n][i])) + (self.matrixB[0][1] * (v1[j][i] - v1[n][i])))
                    vterm = vterm + self.matrixA[n][j] * ((self.matrixB[1][0] * (u1[j][i] - u1[n][i])) + (self.matrixB[1][1] * (v1[j][i] - v1[n][i])))
                uterm = uterm * self.coupling
                vterm = vterm * self.coupling

                dudt = (dudt + uterm)/(self.timescale)
                dvdt = (dvdt + vterm)

                u1[n][i+1] = u1[n][i] + dudt * (t[i+1] - t[i])
                v1[n][i+1] = v1[n][i] + dvdt * (t[i+1] - t[i])

                uterm = 0
                vterm = 0

        return u1, v1

    '''
    parameters:
    u1: Corresponds to u variable, activator variable
    v1: Corresponds to v variable, inhibitory variable
    t: Corresponds to timeseries array
    Example:
    num: Brain Region to graph (If num = 2, brain region 2's activity will be graphed)
    returns:
    None
    '''
    def graphUVMatrix(self, u1, v1, t, rng, ax, upperBoundary, lowerBoundary, pltNum):
        ax[pltNum].plot(t, u1[rng], 'b.-', t, v1[rng], 'r-')
        ax[pltNum].legend(['Activator Variable', 'Recovery Variable'])
        ax[pltNum].axis([0,self.tmax, lowerBoundary, upperBoundary])
        ax[pltNum].grid(True)

    '''
    parameters:
    m1, m2, m3, m4: Represent a 3x3 matrix with any intiailization that user would like
    Example of m1 array:
    m1 = np.array([[0, 0, 1], [1,0,1], [0,1,1]])
    Example of usage:
    matrixA = Brain.KroneckerMatrix(m1, m2, m3, m4)
    returns:
    81 by 81 kronecker matrix formed from the four 3 by 3 matrices given as input
    '''
    @staticmethod
    def KroneckerMatrix(m1, m2, m3, m4):
        m5 = np.kron(m1, m2)
        m6 = np.kron(m3, m4)
        resmatrix = np.kron(m5, m6)
        return resmatrix

    '''
    parameters:
    rng: number of brain regions
    numConnections: number of initial connections to nearest neighbors
    p: rewiring probability
    Example:
    matrixA = Brain.smallworld(90, 2, 0.2)
    returns:
    rng by rng matrix based on the watts-strogatz algorithm
    '''
    @staticmethod
    def smallworld(rng, numConnections, p):
        graph = nx.watts_strogatz_graph(rng, numConnections, p)
        A1 = nx.to_numpy_array(graph)
        return A1

    '''
    parameters:
    rng: number of brain regions
    lowerLim: lower limit for random generation
    upperLim: upper limit for random generation
    Example:
    matrixA = Brain.randSurrogate(90, 0, 1)
    returns:
    rng by rng matrix with random strength of connections
    '''
    @staticmethod
    def randSurrogate(rng, lowerLim, upperLim):
        resmatrix = np.zeros((rng, rng))
        sum1 = 0
        for i in range(0, rng):
            for j in range(0, rng):
                if (i == j):
                    resmatrix[i][j] = 0
                else:
                    temp = random.uniform(lowerLim, upperLim)
                    resmatrix[i][j] = temp
                    sum1 += temp

        print("AVERAGE WEIGHT")
        print(sum1/(rng * rng))
        return resmatrix

    '''
    parameters:
    rng: number of brain regions
    str: base string to form fractals from
    Example:
    matrixA = Brain.fractalConn(90, "101000101000000000101000101000000000000000000000000000101000101000000000101000101")
    returns:
    rng by rng matrix based on fractal connectivity from base string by shifting base string for "rng" times, where rng is size of the network
    '''
    @staticmethod
    def fractalConn(rng, str):
        resmatrix = np.zeros((rng, rng))
        resarray = np.zeros(rng)
        shifted_list = np.zeros(rng)
        for i in range(0, rng):
            resmatrix[0][i] = int(str[i])
        for i in range(0, rng-1):
            resarray = collections.deque(resmatrix[i])
            resarray.rotate(1)
            shifted_list = list(resarray)
            resmatrix[i+1] = shifted_list
        return resmatrix

    '''
    parameters:
    rng: number of brain regions
    str: base string to form fractals from
    lowerLim: lower limit for random generation
    upperLim: upper limit for random generation
    Example:
    matrixA = Brain.realFractalConn(90, "101000101000000000101000101000000000000000000000000000101000101000000000101000101", 0, 1)
    returns:
    similar matrix that the fractalConn method returns, but with random strength of connections rather than 0s and 1s
    '''
    @staticmethod
    def realFractalConn(rng, str, lowerLim, upperLim):
        resmatrix = np.zeros((rng, rng))
        resarray = np.zeros(rng)
        shifted_list = np.zeros(rng)
        for i in range(0, rng):
            resmatrix[0][i] = int(str[i])
        for i in range(0, rng-1):
            resarray = collections.deque(resmatrix[i])
            resarray.rotate(1)
            shifted_list = list(resarray)
            resmatrix[i+1] = shifted_list

        for i in range(0, rng):
            for j in range(0, rng):
                if (resmatrix[i][j] == 1):
                    # add randomized component
                    resmatrix[i][j] = random.uniform(lowerLim, upperLim)
        return resmatrix

    '''
    parameters:
    rng: number of brain regions
    Example:
    matrixA = Brain.generateStar(90)
    returns:
    rng by rng matrix based on a star formation
    '''
    @staticmethod
    def generateStar(rng):
        graph = nx.star_graph(rng)
        A1 = nx.to_numpy_array(graph)
        return A1

    '''
    parameters:
    rng: number of brain regions
    Example:
    matrixA = Brain.randomTree(90)
    returns:
    rng by rng matrix created based on a uniformly random tree
    '''
    @staticmethod
    def randomTree(rng):
        graph = nx.random_tree(rng)
        A1 = nx.to_numpy_array(graph)
        return A1


    '''
    parameters:
    t: timeseries array
    rng: number of brain regions
    lowerLim: lower limit for random generation
    upperLim: upper limit for random generation
    Example:
    uarray, varray = new_brain.randomUV(t, 90, 0, 1)
    returns:
    u and v matrices with random starting u and v values for every brain region
    '''
    def randomUV(self, t, rng, lowerLim, upperLim):
        u1 = np.ones((rng, len(t)))
        v1 = np.ones((rng, len(t)))
        for i in range(0, rng):
            u1[i][0] = random.uniform(lowerLim, upperLim)
            v1[i][0] = random.uniform(lowerLim, upperLim)
        return u1, v1

    '''
    parameters:
    u: u matrix with calculated u values for all brain regions
    v: u matrix with calculated u values for all brain regions
    t: timeseries array
    rng: number of brain regions
    pltNum: the plot number on which to graph the plot
    ax: axis to graph the kuramoto parameters
    lowerBoundary: lower boundary for graph vertical axis
    upperBoundary: upper boundary for graph vertical axis
    maxt: Maximum value of the x axis for the Kuramoto parameter Graph
    Example:
    new_brain.graphR(uarray, varray, t1, 90, axs, 0.8, 2.56, 0, 1)
    returns:
    None
    '''
    def graphR(self, u, v, t, rng, pltNum, ax, threshold, period, lowerBoundary, upperBoundary, maxt):
        r1 = np.zeros(len(t))
        rterm = 0
        T = period
        angleVec = 0
        newT = 0
        numTerms = 0
        for i in range(0, len(t)):
            for j in range(0, rng):
                angleVec = ((2 * math.pi * (math.atan(v[j][i]/u[j][i])))/T)
                newT = (np.e ** (1j * angleVec))
                rterm += newT
            r1[i] = (np.abs(rterm))/rng
            if (r1[i] >= threshold):
                numTerms = numTerms + 1
            rterm = 0

        print("Percent over threshold: ", ((numTerms/len(t)) * 100))
        ax[pltNum].plot(t, r1, 'b.-')
        ax[pltNum].axis([0, maxt, lowerBoundary, upperBoundary])
        ax[pltNum].grid(True)

    '''
    parameters:
    numplots: number of graphs on display
    figSize: size of the figure
    title: title for the figure
    Example:
    fig, axs = new_brain.setUpStreamLit(2, 15, "Network and Graph")
    returns:
    None
    '''
    def setUpStreamLit(self, numplots, figSize, title):
        fig, axs = plt.subplots(numplots, figsize = (figSize,figSize))
        fig.suptitle(title)
        return fig, axs

    '''
    parameters:
    adjMatrix: Adjacency matrix representation of the network
    width: width of the boxes for the heatmap
    axs: The axis to graph the heatmap upon
    colorMap: Color scheme for heat map
    pltNum: Plot Number to graph upon
    Example:
    new_brain.createHeatMap(matrixA, 0.5, axs, "YlGnBu", 1)
    returns:
    None
    '''
    def createHeatMap(self, adjMatrix, width, axs, colorMap, pltNum):
        sns.heatmap(adjMatrix, linewidth = width, ax = axs[pltNum], cmap = colorMap)

def main():
    # Various filenames change the structure of the adjacency matrix
    # data.csv - Chain of brain regions connected together
    # data2.csv - One brain region connected to all other brain regions
    # NOTE" FOR TESTING: Please plot either the heatmap or the r value, unfortunately streamlit is having
    # trou"ble displaying 3 plots at once.
    fileName = 'data.csv'
    matrixB = Brain.generatematrixB((math.pi/2)- 1)
    # Placeholder values if you would like to convert the multiFHN brain region model to a FHN system with a single brain region
    I_ampl1 = 0
    b1 = 0
    tau1 = 0
    u01 = 0
    v01 = 0
    # Sliders that change the values of various parameters
    a = st.slider("A value", 0.0, 1.0, 0.5, 0.01)
    timescale = st.slider("Timescale Value", 0.0, 1.0, 0.05, 0.01)
    coupling = st.slider("Coupling", 0.000, 10.000, 0.500, 0.001)
    tmax = st.slider("Maximum Time", 5, 1000, 10, 5)
    steps = st.slider("Number of steps", 100, 20000, 400, 20)
    num = st.slider("Brain Region Number", 0, 89, 0, 1)
    p = st.slider("Rewiring probability", 0.0000, 1.0000, 0.5000, 0.0001)
    rng = 90
    t1 = np.linspace(0, tmax, steps)
    # Setting up adjacency matrix
    #matrixA = Brain.realFractalConn(rng, "101000101000000000101000101000000000000000000000000000101000101000000000101000101", 0.0001, 0.001)
    matrixA = Brain.smallworld(rng, 2, p)
    # Creates the Brain object
    new_brain = Brain(tmax, a, u01, v01, I_ampl1, b1, tau1, matrixA, matrixB, timescale, coupling)

    #Gets random starting values for u and v
    uarray, varray = new_brain.randomUV(t1, rng, 0, 1)

    # Calculates u and v values for entire timeseries
    uarray, varray = new_brain.calculateUVmatrix(uarray, varray, t1, rng)

    # Sets up streamlit deck
    fig, axs = new_brain.setUpStreamLit(2, 15, "Network and Graph")

    # Graphs both variables
    #new_brain.graphUVMatrix(uarray, varray, t1, 2, axs, -3, 3, 0)

    # Graphs global kuramoto parameter variables
    new_brain.graphR(uarray, varray, t1, rng, 0, axs, 0.8, 2.70, 0, 1, tmax)

    # Creates heatmap, you can comment line 299 and uncomment
    # line 302 to display heatmap (don't run both at once due to streamlit rendering issues)
    new_brain.createHeatMap(matrixA, 0.5, axs, "YlGnBu", 1)

    st.pyplot(fig)

if __name__ == "__main__":
    main()
