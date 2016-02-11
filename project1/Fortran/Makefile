# ***************************************
F90     = gfortran
FFLAGS  = -llapack -lblas
LIBRARY = -L/usr/local/lib
INCLUDE = -I/usr/local/include
DFLAGS  =

# ***************************************
OBJECTS = math.o \
          main.o

MODULES = math.mod

.PHONY: clean

project1 : $(MODULES) $(OBJECTS)
	          $(F90) $(FFLAGS) $(LIBRARY) $(INCLUDE) $(OBJECTS) -o project1.exe

%.o   : %.f90
	       $(F90) -cpp $(DFLAGS) -c $<

%.mod : %.f90
	       $(F90) -cpp $(DFLAGS) -c $<

clean :
	       rm -rf $(OBJECTS) $(MODULES) project1.exe
