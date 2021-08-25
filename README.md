# BrainSimulator
This is a python script and a few other files that can simulate neuron activity based on various adjacency matrices to see if we can simulate seizures. There are still a few features that need to be updated, such as correlation parameters to check for synchronization among neurons. If you have any questions, you can open up an issue and I'll make sure to look at it.

# External Packages
- Matplotlib - For graphing the neuron activity
- Streamlit - Interface to view and change parameters
- Numpy - Numerical calculation
- Scipy - Needed as a prerequisite for other packages
- Pandas - Needed to read files for adjacency matrix
- NetworkX - Used to visualize the network

# Methods
- csvMatrix(inputFile, start) - Creates an adjacency matrix based on a CSV file (Ex: data.csv)
- generatematrixA(n) - Generates a random adjacency Matrix
- generatematrixB(psi) - Generates the matrix B
- calculateUVmatrix(self, u1, v1, t, range1) - Completes the numerical calculations for the u and v arrays for each neuron
- graphUVmatrix(self, u1, v1, t, n, ax) - Graphs the u and v values for a given neuron
- KroneckerMatrix(a, b, c) - Creates the Kronecker Adjacency Matrix
- JellyfishMatrix(startingval) - Creates an adjacency matrix based on the structure of a jellyfish
- randomUV(self, t) - Sets random starting values for u and v for all neurons

# Files Included 
- data.csv - Contains the simple chain adjacency matrix
- data2.csv - Contains the adjacency matrix in which one neuron is connected to all others
- neuron.py - Contains all methods for graphing, visualization, and data analysis
- regular.txt - Real-life data of a rat without seizures
- seizure.txt - Real-life data of a rat with seizures
- samplegraph.py - Graphs the real-life data to find similarities

# Guide to Run
- Start by downloading all of the files above
- To run the file, please use the following command: streamlit run neuron.py in the terminal
- Once the interface pops up, the sliders can be used to change the various parameters

# Future Developments
- Include parameters to measure synchronization easily (detection for seizures) (DONE)
- Include a method to visualize the adjacency matrix (spyplot) (DONE)
- Test further network topologies to see if we can better results that are similar to real data (DONE)
- Use the simulator to attempt to reverse engineer and be able to predict seizures
- Increase the display to make it easier to analyze graphs and networks (DONE)
- Create a downloadable package so it is easier for future researchers to use the software

