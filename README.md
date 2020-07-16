# Ellipsometry
Code for data anlyses on data obtained through an ellipsometry measurement

The  program  is  build  in  python  version  3.7,  the  required  packages  are:  mat-plotlib.pyplot, Numpy, math, cmath, Pandas and lmfit.
To run, two seperatate files need to be present:
A file named ‘materials.txt’ that contains information on the materials that make up the sample. It needs to conatian the name, index of refraction and thickness of all layers making up the sample.
Because there is not necessarily one unique solution that produces a good fit, an initial educated guess needs to be defined for any unknown variable. This will also make sure the algorithm knows where to start with the fit. For the thickness of the substrate and ambient layers an arbitrary number can be added, since this value is not used in any of the calculations.  (It is assumedto be infinite)
The other file needed is the measured data.  See "example_data.phi" for an example on the format the data needs to be in. 
To  run  the  program,  fist  the  path  to  the  ‘materials.txt’  file  and  the  path  tot he measured data needs to be defined.  this can be done by defining a variable ‘locdata’ and ‘locmaterials’ in the code.  (There is a space for this in the firstfew lines of code)
Now the program can be run.  In the command line the user will be prompted for some information:
Fist the program needs to know the number of layers that make up the sample. This is counting the substrate layer but not counting the ambient air.  (So for asilicon waver the number of layers is two:  The silicon substrate and the silicon-oxide on top.)  So far the program supports a system of 1,2 or 3 layers.
Next it will prompt you for the names of the layers and whether to fit any of the variables.  The names of the layers are those inputted in the ‘materials.txt’ 
file and are case-sensitive.  The layers are counted from the substrate up. To fit for a parameter type ‘y’ when prompted.
Now the program will ask fora range over the initial guess where it needs to look for the best value.  Make sure this range stays withing reasonable values (usually a few nanometer) but a big range is no problem.  Try to only fit for one variable at a time.  Asking the program to fit the thickness of the ambient or substrate layer will not have any effect, since these are treated as infinite in the model.
The output will be the values of all parameters as a print statement, and some graphs to show the fit.  The program only fits to the real part of ρ, the program can in principle be adapted to make it fit to one of the other plotted graphs, but this requires some coding.
The outputted data in the command line begins with a table of all variables, then the name of the used model and some statistics about the fit. The last part  is  the  final  value  of  all  variables,  counted  from  the  top  down  (1  for  the ambient, the highest number for the substrate).  n is the real part of the index of refraction, k the imaginary part, d is the thickness of the layer
