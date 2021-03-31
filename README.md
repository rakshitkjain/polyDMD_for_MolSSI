# polyDMD
Simulation system are polydisperse square well chains in a rectangular box

EXECUTION:

make distclean: Removes earlier made executables and logfiles

make: Creates new executables

./MD -h : Shows the help menu using command line parser (readinput.cpp)

./MD -program options : Executes the program for your specified variables

FILES:

Header files contain the class inclusions and the global variable definitions

System.cpp : Initialises the system variables like length of box, number of particles etc, their coordinates, velocities and energies. Also contains the overlap functions called later in the code

Timecalc.cpp: Formulates the timelist and decides which particle pair collide first

Collision.cpp: Proceeds with the collision, moves the particles by the collision time, then redistributes the velocities. Also contains the thermostat function

Rdf.cpp: Formulates the radial distribution function of the system

Main.cpp : Executes all the files and the code main file

Rdf_fin.cpp: Function to calculate the radial distribution function of the coordinates in case the system crashes (NOT READY YET)


For the parameter file, the format should be very specific. Follows a bit of NAMD procedure. 

For multiple nonbonded potentials for the same molecule pair, inner length of next well should be the same as the outer length of the previous well (this is there in the parameter file ".par" extension)

Now we are implementing the false positioning method as well. So the coordinates stored at any time are not the actual coordinates, rather the coordinates at the last update
