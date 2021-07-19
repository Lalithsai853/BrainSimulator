import os
import random
import sys
import numpy as np
import scipy as scpy
import scipy.integrate as spin
from matplotlib import pyplot as plt
import random
import matplotlib.patches as mpatches
import pandas
import math
import streamlit as st
import networkx as nx
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
        for i in range(0,len(t)-2):
            for n in range(0, range1):
                dudt = u1[n][i] - np.power(u1[n][i], 3)/3 - v1[n][i]
                dvdt = u1[n][i] + self.a

                for j in range(0, range1):
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


    @staticmethod
    # parameters: starting value
    # returns: a 90x90 matrix that is shaped like a Jellyfish
    def JellyfishMatrix(startingval):
        resmatrix = np.zeros((90, 90))
        for i in range(0, 89):
            for j in range(0, 89):
                resmatrix[i][j] = startingval
        for i in range(35, 54):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        resmatrix[35][54] = 1
        resmatrix[54][35] = 1

        for i in range(20, 34):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(55, 69):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        resmatrix[20][69] = 1
        resmatrix[69][20] = 1
        resmatrix[34][55] = 1
        resmatrix[55][34] = 1
        for i in range(0, 4):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(5, 9):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(10, 14):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(14, 19):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(70, 74):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(75, 79):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(80, 84):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        for i in range(85, 89):
            resmatrix[i][i+1] = 1
            resmatrix[i+1][i] = 1
        resmatrix[0][23] = 1
        resmatrix[23][0] = 1
        resmatrix[5][26] = 1
        resmatrix[26][5] = 1
        resmatrix[10][29] = 1
        resmatrix[29][10] = 1
        resmatrix[14][32] = 1
        resmatrix[32][14] = 1
        resmatrix[70][58] = 1
        resmatrix[58][70] = 1
        resmatrix[75][61] = 1
        resmatrix[61][75] = 1
        resmatrix[80][64] = 1
        resmatrix[64][80] = 1
        resmatrix[85][67] = 1
        resmatrix[67][85] = 1
        return resmatrix

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

def main():
    # Various filenames change the structure of the adjacency matrix
    # data.csv - Chain
    # data2.csv - One neuron connected to all other neurons
    # KroneckerMatrix and Jellyfish Methods also change the structure of the adjacency matrix
    fileName = 'data.csv'
    matrixA = Neuron.csvmatrixA(fileName, 0.001)
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
    tmax = st.slider("Maximum Time", 5, 100, 10, 5)
    steps = st.slider("Number of steps", 100, 1000, 400, 20)
    num = st.slider("Neuron Number", 0, 89, 0, 1)
    t1 = np.linspace(0, tmax, steps)

    # Creates the neuron object
    new_neuron = Neuron(tmax, a, u01, v01, I_ampl1, b1, tau1, matrixA, matrixB, timescale, coupling)

    #Gets random values for u and v
    uarray, varray = new_neuron.randomUV(t1)

    # Calculates u and v values for entire timeseries
    uarray, varray= new_neuron.calculateUVmatrix(uarray, varray, t1, 89)

    fig, axs = plt.subplots(2)
    fig.suptitle("Network and Graph")

    # Graphs both variables
    new_neuron.graphUVmatrix(uarray, varray, t1, num, axs)

    # Draws the network
    G = nx.from_numpy_matrix(matrixA)
    nx.draw(G, with_labels = True, font_weight = 'bold')
    st.pyplot(fig)


if __name__ == "__main__":
    main()
