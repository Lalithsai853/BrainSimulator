# BrainSimulator
This is a python script and a few other files that can simulate neuron activity based on various adjacency matrices to see if we can simulate seizures. There are still a few features that need to be updated, such as corrletation parameters to check for synchronization among neurons. If ytou have any questions, you can open up an issue and I'll make sure to look at it.

# External Packages
- Matplotlib - For graphing the neuron activity
- Streamlit - Interface to view and change parameters
- Numpy - Numerical calculation
- Scipy - Needed as a prerequisite for other packages
- Pandas - Needed to read files for adjacency matrix

# Methods
- csvMatrix(inputFile, start) - Creates an adjacency matrix based on a CSV file (Ex: data.csv)
- generatematrixA(n) - Generates a random adjacency Matrix
- generatematrixB(psi) - Generates the matrix B
- calculateUVmatrix(self, u1, v1, t, range1) - Completes the numerical calculations for the u and v arrays for each neuron
- graphUVmatrix(self, u1, v1, t, n, ax) - Graphs the u and v values for a given neuron
- KroneckerMatrix(a, b, c) - Creates the Kronecker Adjacency Matrix
- JellyfishMatrix(startingval) - Creates an adjacency matrix based on the structure of a jellyfish
- randomUV(self, t) - Sets random starting values for u and v for all neurons

# Future Developments
- Include parameters to measure synchronization easily (detection for seizures)
- Include a method to visualize the adjacency matrix (spyplot)
- Test further network topologies to see if we can better results that is similar to real data
- Use the simulator to attempt to reverse engineer and be able to predict seizures
- Increase the display to make it easier to analyze graphs and networks

