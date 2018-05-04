# Meniscus

This despository was created to share a Matlab program which I initially created in 2014-2015 
and continued to improve in 2017/2018 to construct equilibrium shapes of axisymmetric free surface problems 
using continuation method.
(pending drops, attached bubbles, sessile drops, etc, all geometries generically called menisci) 

The whole program is contained in "meniscus.m" which is written as a "class" in matlab-object style.

The current distribution includes five examples of scripts which demonstrate the software :

Example_PendingDrop.m
Example_PendingDrop_ContactAngle.m  (working but not 100% satisfactory, convergence issues, to be improved).
Example_SessileDrops_Pined.m
Example_SessileDrops_Angle.m
Example_SessileDrops_Water_PMMS.m


Main functionalities are :

1/ 
m = meniscus('type',Npoint,[parameters])  
    -> Constuctor for an initial profile.
The types currently operational are :
    'flat' - > initial flat surface ; here [parameters] = [Radius] 
    'flatinv' - > initial flat surface with inverted orientation (see NB below); here parameter = [Radius] 
    'cap' - > initial flat surface ; here [parameters] = [Radius,thetas] ; thetas in DEGREES
    'capinv' - > initial flat surface with inverted orientation (see NB below); here [parameters] = [Radius,thetas] 
    'cylinder' -> cylindrical surface (use this for liquid bridge). Here [parameters] = [Radius,Length].

NB according to a common convention, the contour encircles the drop/bubble/inclusion in the trigonometric sense.
Accordingly :  
    -   A meniscus in "direct" definition starts from the pined point and ends on the axis 
        (use this if the fluid is "below" the meniscus contour ; for sessile drops or attached bubbles)
    -   A a meniscus in "inverted" definition starts at the axis and ends at the pinned point 
        (use this if the fluid is "above" the meniscus contour ; for hanging drops and bubbles at ceiling)

2/ 
m=m.step('type',value)
    -> Calculation of a new equilibrium shape, starting from a previous shape, and specifying a variation of one parameter.
    'type' can be (currently) :
        'P' -> specify a new value of P
        'V' -> specify a new value of V
        'dP' -> specify an INCREMENT in P (new value is P+dP)
        'dV' -> specify an INCREMENT in V (new value is V+dV)
        'dS' -> specify an increment in ARCLENGTH in the P-V plane.

3/ 
m = m.loop('type',dX,N)
    -> Loop over "m.step" to construct a family of menisci shapes by continuation.
    'type' can be (currently) :
           'dP' (continuation by varying the pressure)
           'dV' (continuation by varying the volume)
           'dS' (arglength-continuation in the P-V plane).


NOTE on graphical output.
the software includes an number of possibilities for graphical output, and is designed to be highly customizable. 
the "basic" outputs are :
   - figure 10 : representation of the menisci shapes in "R-Z" plane 
        (by default the software will plot each meniscus calculated in "step" mode, and only the cases corresponding to 
            "start/end/stability change" in "loop" mode) 
   - figure 20 : representation of the menisci characteristic in a P-V diagram 
                    (the program will plots points in "step" mode and curves in "loop" mode). 
   A number of other possibilities are hidden in the program... for instance figure 121 will plot thetas as function of V.
   To add a new plot simply add the corresponding number to m.whichfigures 
   (for instance    m.whichfigures = [m.whichfigures 121]  to add figure 121 as described).

HISTORY :

current version is "meniscus.m" in the root repository.

Two previous versions of the software are provided here
- Version 1.4 in "old-style" matlab, including several scripts
- Verson 2.5 in object-oriented matlab, with only two short example scripts 

An explanation of the methods under the form of an unfinished paper "JFM_Menisci" is provided in PDFs directory. 

All this is kindly shared without any guarranty ; anyone who wants to help me continue this project is warmly thanked !

D. Fabre, october 4, 2017.

How to install and use this software ?

    If you just want to install the current stable version, simply type the following command in terminal (after making sure the git command is available on your system)

git clone https://github.com/erbafdavid/Meniscus

    If you want to participare to the project you should create a git account (...)
