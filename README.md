# BrainSimulator
This is a python script and a few other files that can simulate brain activity based on various brain network topologies to observe the impact of network structure on seizure activity. The base of the simulator is built based upon the Fitz-Hugh Nagumo oscillator, which was chosen for its ability to model excitement and inhibition in a brain network. Several features have been updated which can be used to check for synchronization among neurons, and this synchronization is used to tell whether or not a seizure occurred. If you have any questions, you can open up an issue and I'll make sure to look at it.

# External Packages
- Matplotlib - For graphing the neuron activity
- Streamlit - Interface to view and change parameters
- Numpy - Numerical calculation
- Scipy - Needed as a prerequisite for other packages
- Pandas - Needed to read files for adjacency matrix
- NetworkX - Used to visualize and generate some of the networks

# Features
- Graphs brain activity based on the Fitz-Hugh Nagumo system and connections between brain regions
- Tracks synchronization of brain regions through the Kuramoto Parameter to observe whether or not a seizure occurred
- Plots the heatmap of the network for visualization
- Ready-to-go networks such as the Watts-Strogatz, Kronecker Adjacency Matrices/Network Structures, and more
- Gives the user the ability to add custom networks through csv files
- Clean and simple UI with sliders through streamlit for user experience
- Can save the FHN, Kuramoto Parameter, and heatmap plots automatically given the location to store the file

# Methods (Detailed documentation about usage, parameters, and return objects can be found in the comments in the source file)
- csvMatrixA(inputFile, start) - Creates an adjacency matrix based on a CSV file (Ex: data.csv)
- generatematrixA(n) - Generates an adjacency Matrix with randomly chosen connections
- generatematrixB(psi) - Generates the matrix B that is used for network 
- calculateUVmatrix(self, u1, v1, t, range1) - Completes the numerical calculations for the u and v arrays for each brain region
- KroneckerMatrix(m1, m2, m3, m4): - Creates an 81 by 81 Kronecker Adjacency Matrix based on four 3 by 3 matrices
- smallworld(rng, numConnections, p): - Creates a small world network based on the Wattz-Strogatz algorithm 
- randSurrogate(rng, lowerLim, upperLim): - Creates an adjacency matrix with random connection strengths
- fractalConn(rng, str): - Uses a base string and then shifts this right "rng" times to form a fractal network
- realFractalConn(rng, str, lowerLim, upperLim): - Uses a base string and then shifts this right "rng" times to form a fractal network, but the strength of the connection is randomly assigned
- generateStar(rng): - Generates a rng by rng star network
- randomTree(rng): - Generates a rng by rng uniform tree network
- randomUV(self, t, rng, lowerLim, upperLim): - Sets random starting values for u and v for all neurons
- graphUVMatrix(self, u1, v1, t, rng, ax, upperBoundary, lowerBoundary, pltNum, fig=None, save_path=None): - Graphs the u and v values for a given brain region
- graphR(self, u, v, t, rng, pltNum, ax, threshold, period, lowerBoundary, upperBoundary, maxt, fig=None, save_path=None): - Graphs the Kuramoto parameter to track synchronization across all brain regions
- createHeatMap(self, adjMatrix, width, ax, colorMap, pltNum, fig=None,save_path=None): - Creates the heatmap of the brain network showing strength of connections and general structure
- setUpStreamLit(self, numplots, figSize, title): - Sets up the streamlit deck to show sliders + visualizations

# Networks included in source file
- Kronecker matrices
- Random connection network
- Random strength of connection network
- Watts-Strogatz Network
- Uniform Tree Network
- Star Network
- Fractal Connectivity Network
- Fractal Connectivity Network with random strength of connections

# Files Included 
- data.csv - Contains the simple chain adjacency matrix
- data2.csv - Contains the adjacency matrix in which one neuron is connected to all others
- neuron.py - Contains all methods for graphing, visualization, and data analysis
- regular.txt - Real-life data of a rat without seizures
- seizure.txt - Real-life data of a rat with seizures
- samplegraph.py - Graphs the real-life data to find similarities

# Slider Information
The following is information about each of the sliders
- A value: represents the strength of excitement of the variables
- Timescale Value: represents the scale between the sudden excitement of the activator variable and inhibition of the inhibitor variable 
- Coupling: represents the strength of interaction between the activator and inhibitor variables
- Maximum Time: The maximum amount of time units for the simulation
- Number of Steps: The number of steps within the maximum amount of time units 
- Brain Region Number: The brain region to graph (For the graph of the Kuramoto parameter, it is the same regardless of brain region selection
- Rewiring Probability: The probability to rewire the connections to a random node in the Watts-Strogatz algorithm

# Guide to Run
- Start by downloading all of the files above
- Navigate to the correct directory through the terminal
- To run the file, please use the following command: streamlit run neuron.py in the terminal
- Once the interface pops up, the sliders can be used to change the various parameters

# Future Developments
- Create a downloadable package so it is easier for future researchers to use the software
- Expand the documentation even further to prevent issues for future developers
- Use the simulator to attempt to reverse engineer and be able to predict seizures

