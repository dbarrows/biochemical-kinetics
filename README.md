Modeling Arrays of Reactions Software
=====================================

Ryerson University, MTH 40 A/B Project - Dexter Barrows


MARS will take a Systems Biology Markup Language (SBML) file and model individual species vs time relationships on a single graph. The underlying algorithm used is the Stochastic Simulation Algorithm (SSA), also known as the Gillespie algorithm, popularized by Dan Gillespie in 1977. It is also capable of implementing the stochastic modeling methods: Tau-leaping adaptation of SSA, and the Chemical Lengevin Equation (CLE). Additionally, it can implement the  Reaction Rate Equation (RRE) method using differential equations. 

MARS also has the ability to generate a histogram of the final amounts of each species in the system based on about 10,000 simulated trajectories. This capability is enabled with the use of a function argument detailed below, and is compatible with the SSA, Tau-leaping, and CLE methods. MARS can also use a CUDA-supported GPU to rapidly simulate multiple trajectories. Any Nvidia GPU with Compute Capability 1.3 or higher will work. If you are running MARS with a GPU on a Windows platform, it will be subject to the time limitation mentioned [here](https://developer.nvidia.com/cuda-faq) under "Q: What Is The Maximum Kernel Execution Time?". If you do need to run simulations longer than 5 seconds, and are running Windows 7, you can disable or extend the watchdog timer (also called the Timeout Display Recovery or TDR) by creating/editing the registy key using the instructions found [here](http://msdn.microsoft.com/en-us/library/windows/hardware/ff569918.aspx).

MARS is written in MATLAB code. This version includes libraries to run on OS X or Windows based platforms with Intel x86 processors.


Behaviour
---------

To obtain plotting data for a SBML file named 'filename' using the default values:

* Time (to run the simulation for) = 1000
* Values recorded every 20 steps

Simply enter the following:

        >> [Y, t] = MARS( 'filename' )

Where _Y_ is the name of the matrix to contain the species amount data, and _time_ is the vector containing the time-step data, And MARS will do the rest. The matrix _Y_ will be formatted such that each row contains the species amount data for a particular species after each step (so the first row of Y might look like [245 244 243 243 242 ...], where each entry shows the number of species #1 present after step 1, then step 2, etc.). The vector _time_ will be the system's time elapsed up to the end of that step (so it might look like [0.256 0.5822 0.9567 ...]).

To run 10,000 trajectories and obtain the resulting data (perhaps to plot using a histogram):

        >> [Y] = MARS( 'filename', 'Hist')
        
Where _Y_ is the name of the matrix to contain the species amount data. It will be formatted such that each column contains the species amount data for a particular species after each step (so the first column of Y might look like [245 244 243 243 242 ...]', where each entry shows the number of species #1 present after trajectory #1 runs it course, then after trajectory #2, etc.). 

 
Optional single arguments:

*   To run the program in 'verbose' mode, which will enable diagnostic outputs while    setting up and running the simulation:

        'Verbose'

    For example, to run a simulation on a system in the SBML file 'example.xml':

        >> MARS( 'example.xml' , 'Verbose' )

*   To split the results graph into individual graphs for each species present:

        >> MARS( 'example.xml' , 'Split' )

*   To use the RRE method:

        >> MARS( 'example.xml' , 'RRE' )

*   To use the system's CUDA-enabled GPU (to generate or plot histogram data):

        >> MARS( 'example.xml' , 'GPU' )
        
    Note that the 'GPU' argument automatically implies the 'Hist' argument, so you need not also include it. You can certainly have it there too, but it's not necessary.


MARS implements the name-value pair paradigm found in many MATLAB programs for several argument types:

*   To run the simulation for a specific period of time:

        'Time', time_final

    Where 'time_final' is an integer representing the time you want the the simulation to run for. For example, for the SBML file 'example.xml' and a running time of 10^6:

        >> MARS( 'example.xml' , 'Time' , 10^6 )

*   To record the system state (which will also affect how fine-grained the graph will appear) every certain number of steps:

        'Record', record_step_size

    Where 'record_step_size' is an integer representing how many steps you want to allow the algorithm to take before next recording the state of the system. For example, for the SBML file 'example.xml' and recording results every 20 steps:

        >> MARS( 'example.xml' , 'Record' , 20 )

    Note that this argument is incompatible with the 'Hist' argument, so 'Record' will be ignored if 'Hist' is also used.

*   To use the Tau-leaping method, with an automatic selection for the value of tau (note that this estimate will be extremely crude, use with caution):

        >> MARS( 'example.xml' , 'Tau')

*   To use the Tau-leaping method, with a specified value for tau of 1.3*10^-3:

        >> MARS( 'example.xml' , 'Tau', 1.3*10^-3 )

*   To use the CLE method, with an automatic selection for the value of tau (note again that this estimate will be extremely crude, use with caution):

        >> MARS( 'example.xml' , 'CLE')

*   To use the CLE method, with a specified value for tau of 2.7*10^-3:

        >> MARS( 'example.xml' , 'CLE', 2.7*10^-3 )

*   To produce a graph of all species data in the system, either with single-trajectory generation or histogram generation:

        >> MARS( 'example.xml' , 'Graph')

*   To graph the data only for specific species in the system, follow the 'Graph' argument with a vector containing the index numbers of the species you want to graph using the order of species from the SBML file you provided. For example to graph data for species #1, #4, and #5 in a system with 5 or more species present:

        >> MARS( 'example.xml' , 'Graph', [1 4 5])


Any valid combinations of arguments are also allowed. Order is not important except for the name-value pair sets of arguments.

For example, to run a simulation with the default final time and record step size but with verbose output and a histogram for the system in 'example.xml', enter:

        >> MARS( 'example.xml' , 'Verbose' , 'Hist' )

Or, to run a simulation with a running time of 30,000 and default record step size and with a histogram for the system in 'example.xml', enter:

        >> MARS( 'example.xml' , 'Hist' , 'Time' , 3*10^4 )

Or, to be really fancy to run a simulation with a running time of 30,000, record step size of 30, with verbose output, using CLE with a tau value of 4.567*10^-8 and with a histogram for the system in 'example.xml', then graph the system only for species #4 and #7, do all the trajectories on the system's GPU, and suppress the output of the returned data, enter:

        >> MARS( 'example.xml' , 'Time' , 3*10^4 , 'CLE', 4.567*10^-8, 'Hist' , 'Verbose' , 'Record' , 30 , 'Graph', [4 7], 'GPU');

...or some other ordered combination of the given arguments.


Citations / Acknowledgements
----------------------------

This software makes use of (and I therefore make no attempt to take credit for) the excellent software listed here:

* libSBML http://sbml.org/Software/libSBML
* SBMLToolbox http://sbml.org/Software/SBMLToolbox

As per the the request of the authors listed, these pieces of software will be formally cited in the paper to be authored by myself when this project manifests itself as a thesis paper.
