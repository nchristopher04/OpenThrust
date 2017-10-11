OpenThrust
========

OpenThrust is a software created by members of the Waterloo Rocketry design team at the University of Waterloo. The aim of the project is to provide an accurate thrust curve and apogee prediction for hybrid engine rockets that use NOS as their oxidizer.

Report: https://www.overleaf.com/read/kfprjrjbkmbh

Features
--------

- Accurate thrust curve prediction
- Apogee prediction
- More to be added...

Planned Future Additions
------------------------

- Easier to change rocket parameters (through cfg file maybe)
- Non constant OF ratio
- Documentation of current code
- Code cleanup, making everything more readable to users (maybe following a style guide if its not a huge hassle
- Bugfixes
- Optimization

Using the Program
-----------------

Right now to use the program, you should first generate an RPA table using RPA software found here: http://www.propulsion-analysis.com/index.htm. For our purposes, the free version worked fine.
The table should have increments of 10 Psi for the pressure and 0.1 for the OF ratio.
Then you should go through the code files and switch the values to correspond to your rocket. After that build and run the program.

Contribute
----------

- Issue Tracker: github.com/nchristopher04/OpenThrust/issues
- Source Code: github.com/nchristopher04/OpenThrust

License
-------
