# polyDMD
Run a discontinuous potential Molecular Dynamics Simulation of any kind of system!

EXECUTION:

make distclean: Removes earlier made executables and logfiles

make: Creates new executables

./MD -h : Shows the help menu using command line parser (readinput.cpp)

./MD -<program options> : Executes the program for your specified variables

FILES:

*.H : Header files contain the class inclusions and the global variable definitions

SYSTEM.CPP : Initialises the system variables like length of box, number of particles etc, their coordinates, velocities and energies. Also contains the overlap functions called later in the code

TIMECALC.CPP : Formulates the timelist and decides which particle pair collide first

COLLISION.CPP : Proceeds with the collision, moves the particles by the collision time, then redistributes the velocities. Also contains the thermostat function

RDF.CPP : Formulates the radial distribution function of the system

MAIN.CPP : Executes all the files and the code main file

RDF_FIN.CPP : Function to calculate the radial distribution function of the coordinates in case the system crashes (NOT READY YET)


For the parameter file, the format should be very specific. Follows NAMD procedure. 

For multiple nonbonded potentials for the same molecule pair, inner length of next well should be the same as the outer length of the previous well (this is there in the parameter file ".par" extension)

Now we are implementing the false positioning method as well. So the coordinates stored at any time are not the actual coordinates, rather the coordinates at the last update.

The input files needed are PDB (Protein Data Bank), PSF (Protein Structure Format) and PAR (Parameter File). Use ./MD -h for input instructions.
