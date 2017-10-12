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
- Adding a plotter to program
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

Disclaimer
----------
This program merely models the thrust and apogee of a rocket using research and data. It can not guarantee the safety or functionality of any rocket being built. We do not take any responsibility for any mistakes, accidents, disasters, etc. that can occur while building or testing your rocket. This program is meant as an aide and should not be used as the sole source of data. Double check all your calculations and always follow proper safety procedures. 

License
-------
OpenThrust is licensed under GPL 3. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
