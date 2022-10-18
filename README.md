# BrainSimulator
This is a python script and a few other files that can simulate brain activity based on various brain network topologies to observe the impact of network structure on seizure activity. Several features have been updated which can be used to check for synchronization among neurons. If you have any questions, you can open up an issue and I'll make sure to look at it.

# External Packages
- Matplotlib - For graphing the neuron activity
- Streamlit - Interface to view and change parameters
- Numpy - Numerical calculation
- Scipy - Needed as a prerequisite for other packages
- Pandas - Needed to read files for adjacency matrix
- NetworkX - Used to visualize and generate some of the networks

# Methods (Detailed documentation can be found in the comments in the source file)
- csvMatrixA(inputFile, start) - Creates an adjacency matrix based on a CSV file (Ex: data.csv)
- generatematrixA(n) - Generates an adjacency Matrix with randomly chosen connections
- generatematrixB(psi) - Generates the matrix B that is used for network 
- calculateUVmatrix(self, u1, v1, t, range1) - Completes the numerical calculations for the u and v arrays for each brain region
- graphUVmatrix(self, u1, v1, t, n, ax) - Graphs the u and v values for a given brain region
- KroneckerMatrix(m1, m2, m3, m4): - Creates an 81 by 81 Kronecker Adjacency Matrix based on four 3 by 3 matrices
- smallworld(rng, numConnections, p): - 
- randSurrogate(rng, lowerLim, upperLim): - Creates an adjacency matrix with random connection strengths
- fractalConn(rng, str): - Uses a base string and then shifts this right "rng" times to form a fractal network
- realFractalConn(rng, str, lowerLim, upperLim): -Uses a base string and then shifts this right "rng" times to form a fractal network, but the strength of the connection is randomly assigned
- generateStar(rng): - Generates a rng by rng star network
- randomTree(rng): - Generates a rng by rng uniform tree network
- randomUV(self, t, rng, lowerLim, upperLim): - Sets random starting values for u and v for all neurons
- graphR(self, u, v, t, rng, pltNum, ax, threshold, period, lowerBoundary, upperBoundary, maxt): - Graphs the Kuramoto parameter to track synchronization across all brain regions
- setUpStreamLit(self, numplots, figSize, title): - Sets up the streamlit deck to show sliders + visualizations
- createHeatMap(self, adjMatrix, width, axs, colorMap, pltNum): - Creates the heatmap of the brain network showing strength of connections and general structure
# Files Included 
- data.csv - Contains the simple chain adjacency matrix
- data2.csv - Contains the adjacency matrix in which one neuron is connected to all others
- neuron.py - Contains all methods for graphing, visualization, and data analysis
- regular.txt - Real-life data of a rat without seizures
- seizure.txt - Real-life data of a rat with seizures
- samplegraph.py - Graphs the real-life data to find similarities

# Guide to Run
- Start by downloading all of the files above
- Navigate to the correct directory through the terminal
- To run the file, please use the following command: streamlit run neuron.py in the terminal
- Once the interface pops up, the sliders can be used to change the various parameters

# Future Developments
- Create a downloadable package so it is easier for future researchers to use the software
- Expand the documentation even further to prevent issues for future developers
- Use the simulator to attempt to reverse engineer and be able to predict seizures
