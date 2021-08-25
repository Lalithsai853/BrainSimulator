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
import streamlit as st
import networkx as nx
import seaborn as sns
from kuramoto import Kuramoto
class Neuron:
    # Creates the neuron object
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

    # parameters: csv file and starting value for all connections
    # returns: adj matrix
    @staticmethod
    def csvmatrixA(inputFile, start):
        arrA = np.zeros((90,90))

        df = pandas.read_csv(inputFile)
        row_count, column_count = df.shape
        for i in range(1, row_count):
            arrA[df['FromNeuron'][i]][df['ToNeuron'][i]] = 1
        for i in range(0,89):
            for j in range(0,89):
                if arrA[i][j] != 1:
                    arrA[i][j] = start
        return arrA

    @staticmethod
    # parameters: number of neurons
    # returns: randomly connected adj matrix
    def generatematrixA(n):
        arrA = np.zeros((90,90))
        for i in range(0, n):
            num1 = randint(1, 90)
            num2 = randint(1, 90)
            arrA[num1][num2] = 1
            arrA[num2][num1] = 1
        return arrA


    @staticmethod
    # parameters: Psi
    # returns: B array
    def generatematrixB(psi):
        arrB = np.ones((2,2))
        arrB[0][0] = math.cos(psi)
        arrB[0][1] = math.sin(psi)
        arrB[1][0] = (-1 * math.sin(psi))
        arrB[1][1] = math.cos(psi)
        return arrB


    # parameters: empty u and v arrays, time array, and range
    # range chooses up to which neuron to go to, (Ex: for kronecker matrix, you
    # should only go to 81 while for other adj. matrices, 90 should be the range)
    # returns: U and v arrays with calculated values
    def calculateUVmatrix(self, u1, v1, t, range1):
        uterm = 0
        vterm = 0
        for i in range(0,len(t)-1):
            for n in range(0, range1):
                dudt = u1[n][i] - np.power(u1[n][i], 3)/3 - v1[n][i]
                dvdt = u1[n][i] + self.a
                #print("U and V derivatives after power")
                #print(str(dudt) + "  " + str(dvdt))
                #assert(dudt <= 50 and dudt >= -50)
                for j in range(0, range1):
                    uterm = uterm + self.matrixA[n][j] * ((self.matrixB[0][0] * (u1[j][i] - u1[n][i])) + (self.matrixB[0][1] * (v1[j][i] - v1[n][i])))
                    vterm = vterm + self.matrixA[n][j] * ((self.matrixB[1][0] * (u1[j][i] - u1[n][i])) + (self.matrixB[1][1] * (v1[j][i] - v1[n][i])))
                    #print("A values and u/B values")
                    #print(str(self.matrixA[n][j]) + "  " + str(((self.matrixB[0][0] * (u1[j][i] - u1[n][i])) + (self.matrixB[0][1] * (v1[j][i] - v1[n][i])))))

                uterm = uterm * self.coupling
                vterm = vterm * self.coupling

                #print("U and V Terms after Coupling")
                #print(str(uterm) + "  " + str(vterm))

                dudt = (dudt + uterm)/(self.timescale)
                dvdt = (dvdt + vterm)

                #print("Derivatives after Timescale and adding u and v terms")
                #print(str(dudt) + "  " + str(dvdt))

                u1[n][i+1] = u1[n][i] + dudt * (t[i+1] - t[i])
                v1[n][i+1] = v1[n][i] + dvdt * (t[i+1] - t[i])

                #print("U and V values")
                #print(str(u1[n][i+1]) + "  " + str(v1[n][i+1]))
                uterm = 0
                vterm = 0

        return u1, v1


    # parameters: U array, v array, time array, neuron number, and ax, which is
    # the object used to graph the subplot
    # returns: streamlit plot (upper and lower boundaries can be changed)
    def graphUVmatrix(self, u1, v1, t, n, ax):
        ax[0].plot(t, u1[n], 'b.-', t, v1[n], 'r-')
        upperBoundary = 3
        lowerBoundary = -3
        ax[0].legend(['Activator Variable', 'Recovery Variable'])
        ax[0].axis([0,self.tmax, lowerBoundary, upperBoundary])
        ax[0].grid(True)


    # parameters: a, b, c (integers) depending on the values of a, b, and c, an overflow error may occur
    # The structure can be changed by changing the orderings of a, b, and c below
    # returns: a 81x81 matrix
    @staticmethod
    def KroneckerMatrix(a, b, c):
        matrix1 = np.array([[b, c, b], [c,b,a], [b,a,b]])
        matrix2 = np.array([[a, a, c], [b,b,c], [b,c,a]])
        matrix3 = np.array([[b, a, b], [c,c,a], [b,b,b]])
        matrix4 = np.array([[c, b, b], [a,b,b], [b,a,b]])

        matrix5 = np.kron(matrix1, matrix2)
        matrix6 = np.kron(matrix3, matrix4)

        resmatrix = np.kron(matrix5, matrix6)
        return resmatrix


    # Don't have empirical data currently
    def empiricalNetwork():
        pass

    # Generates completely random matrix with random weights
    @staticmethod
    def randSurrogate():
        resmatrix = np.zeros((90, 90))
        lst = [0.00001, 0.0001, 0.0001, 0.01, 0.01, 1]
        for i in range(0, 90):
            for j in range(0, 90):
                if (i == j):
                    resmatrix[i][j] = 0
                else:
                    temp = random.randrange(5)
                    temp = lst[temp]
                    resmatrix[i][j] = temp
        return resmatrix

    # Creates fractal connectivity topology with weights only 1s and 0s
    @staticmethod
    def fractalConn(a):
        b = str(a)
        resmatrix = np.zeros((81, 81))
        resarray = np.zeros(81)
        shifted_list = np.zeros(81)
        b = "101000101000000000101000101000000000000000000000000000101000101000000000101000101"
        for i in range(0, 81):
            resmatrix[0][i] = int(b[i])
        for i in range(0, 80):
            resarray = collections.deque(resmatrix[i])
            resarray.rotate(1)
            shifted_list = list(resarray)
            resmatrix[i+1] = shifted_list
        return resmatrix

    # Creates fractal connectivity topology with real weights
    @staticmethod
    def realfractalConn():
        lst = [0.00001, 0.0001, 0.0001, 0.01, 0.01, 1]
        b = str(a)
        resmatrix = np.zeros((81, 81))
        resarray = np.zeros(81)
        shifted_list = np.zeros(81)
        b = "101000101000000000101000101000000000000000000000000000101000101000000000101000101"
        for i in range(0, 81):
            resmatrix[0][i] = int(b[i])
        for i in range(0, 80):
            resarray = collections.deque(A[i])
            resarray.rotate(1)
            shifted_list = list(resarray)
            resmatrix[i+1] = shifted_list

        for i in range(0, 81):
            for j in range(0, 81):
                if (resmatrix[i][j] == 1):
                    temp = random.randrange(5)
                    temp = lst[temp]
                    resmatrix[i][j] = temp
                # add randomized component
        return resmatrix

    # Uses Watts-Strogatz Algorithm
    @staticmethod
    def smallworld(p):
        graph = nx.watts_strogatz_graph(90, 3, p)
        A1 = nx.to_numpy_array(graph)
        return A1

    # Generates star topology
    @staticmethod
    def generateStar():
        graph = nx.star_graph(90)
        A1 = nx.to_numpy_array(graph)
        return A1

    # Generates random tree topology
    @staticmethod
    def randomTree():
        graph = nx.random_tree(90)
        A1 = nx.to_numpy_array(graph)
        return A1


    # parameters: time array
    # returns: u and v arrays with random starting values
    def randomUV(self, t):
        u1 = np.ones((90, len(t)))
        v1 = np.ones((90, len(t)))
        for i in range(0, 89):
            if (i < 20):
                u1[i][0] = random.uniform(-1 , 0)
                v1[i][0] = random.uniform(-1 , 0)
            elif (i >= 20 and i < 40):
                u1[i][0] = random.uniform(-1 , 0)
                v1[i][0] = random.uniform(-1 , 0)
            elif (i >= 40 and i < 60):
                u1[i][0] = random.uniform(-1 , 0)
                v1[i][0] = random.uniform(-1 , 0)
            elif (i >= 60 and i <= 89):
                u1[i][0] = random.uniform(-1 , 0)
                v1[i][0] = random.uniform(-1 , 0)
        return u1, v1

    # parameters: U array, with values for u
    # returns: correlation matrix
    def corrMatrix(self, u, axs):
        corr_matrix = np.corrcoef(u)
        print(corr_matrix)

    # parameters: u array
    # returns: correlation parameter
    def pCorrelation(self, u, n1, n2):
        r, p = stats.pearsonr(u[n1], u[n2])
        print("Correlation Coefficient", r)

    # Calculates global kuramoto paramater R given u and v arrays
    def graphR(self, maxt, u, v, t, ax):
        r1 = np.zeros(len(t))
        rterm = 0
        T = 2.56

        for i in range(0, len(t)):
            for j in range(0, 90):
                rterm += (((math.e**((2 * math.pi)/T * (((math.atan(v[j][i]/u[j][i])))*1j).real))))
                #print(((2 * math.pi)/T * (((math.atan(v[j][i]/u[j][i])))*1j)))
                #print("Divider")
                #print((math.e**((2 * math.pi)/T * (((math.atan(v[j][i]/u[j][i])))*1j))).real)
            print(rterm)
            rterm = rterm/90
            r1[i] = rterm
            rterm = 0

        ax[1].plot(t, r1, 'b.-')
        upperBoundary = 1
        lowerBoundary = 0
        ax[1].axis([0,maxt, lowerBoundary, upperBoundary])
        ax[1].grid(True)




def main():
    # Various filenames change the structure of the adjacency matrix
    # data.csv - Chain
    # data2.csv - One neuron connected to all other neurons
    # KroneckerMatrix and Jellyfish Methods also change the structure of the adjacency matrix
    # NOTE FOR TESTING: Please plot either the heatmap or the r value, unfortunately streamlit is having
    # trouble displaying 3 plots at once.
    fileName = 'data.csv'
    matrixB = np.ones((2,2))
    matrixB = Neuron.generatematrixB((math.pi/2)- 1) #Psi Paramater according to Chaos Paper

    #These values are for a single neuron system (usually not used, should keep these values as default)
    u01 = 0
    v01 = 0
    I_ampl1= None
    b1 = None
    tau1 = None

    # Sliders that change the values of various parameters according to the Chaos Paper
    a = st.slider("A value", 0.0, 1.0, 0.5, 0.01)
    timescale = st.slider("Timescale Value", 0.0, 1.0, 0.05, 0.01)
    coupling = st.slider("Coupling", 0.0, 1.0, 0.5, 0.01)
    tmax = st.slider("Maximum Time", 5, 1000, 10, 5)
    steps = st.slider("Number of steps", 100, 20000, 400, 20)
    num = st.slider("Neuron Number", 0, 89, 0, 1)
    p = st.slider("Rewiring probability", 0.0000, 1.0000, 0.5000, 0.0001)
    t1 = np.linspace(0, tmax, steps)

    matrixA = Neuron.smallworld(p)

    # Creates the neuron object
    new_neuron = Neuron(tmax, a, u01, v01, I_ampl1, b1, tau1, matrixA, matrixB, timescale, coupling)

    #Gets random values for u and v
    uarray, varray = new_neuron.randomUV(t1)

    # Calculates u and v values for entire timeseries
    uarray, varray = new_neuron.calculateUVmatrix(uarray, varray, t1, 90)

    fig, axs = plt.subplots(2, figsize = (15,15))
    fig.suptitle("Network and Graph")

    # Graphs both variables
    new_neuron.graphUVmatrix(uarray, varray, t1, num, axs)

    new_neuron.graphR(tmax, uarray, varray, t1, axs)
    # Creates heatmap, you can uncomment the line 333 and uncomment
    # line 336 to display heatmap (don't run both at once due to streamlit rendering issues)
    #sns.heatmap(matrixA, linewidth = 0.5, ax = axs[1])

    #new_neuron.pCorrelation(uarray, 32, 54)
    #new_neuron.corrMatrix(uarray, axs)

    # Draws the networks, not very useful
    #G = nx.from_numpy_matrix(matrixA)
    #nx.draw(G, with_labels = True, font_weight = 'bold')

    st.pyplot(fig)


if __name__ == "__main__":
    main()
