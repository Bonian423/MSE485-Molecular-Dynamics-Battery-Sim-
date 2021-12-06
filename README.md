# MSE485-Molecular-Dynamics-Battery-Sim-

MSE 485 Final Project: simulating Lithium based battery 
Purpose: simulate charging process using monte carlo method


Simulation structure:
    TWO BOXES:  one for anode and one for cathode
    Universal electric field and drives lithium ions in one direction.
    LiCoO2 lattice:
        fixed CO and O atoms
        free Li ions
        reintroduction of Li ions during the charging process
    LiC lattice:
        Fixed C atoms
        free Li ions
        Gradual introduction of Li atoms over time
    No implementation of separator.


Simulation calculation construction:
    atomic interactions: 
        LJ interaction:
            Traditional LJ force calculation for non-charged particles
        Electrostatic interaction:
            Implementation through Ewald summation
            kspace calculation for Ewal sum
            coloumb self interaction


Observables:  
    Diffusion coefficient
    Kinetic and Potential Energy
    Eward Energy
    Charging temperature
    Structure factor 
    Velocity autocorrelation


Python code files documentation:
    Lattice_CoO2.py:
        python class, construction of lattice of LiCoO2 
    Lattice_Graphite.py:
        python class, construction of lattice of LiC 
    PLOTTING_DATA.py:
        visualization of simulation observables in python
    Trail 1~18.spydata:
        simulation observables data storage
    total_MD_code_V2:
        final version of calculation fucntions of forces and observables
    total_MD_implement_CoO2_V2:
        final version of simulation executing code for LiCoO2 lattice
    total_MD_implement_Graphite_V2:
        final version of simulation executing code for LiC lattice
    
