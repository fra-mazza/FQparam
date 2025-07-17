# Variabili
FC = gfortran
FCFLAGS = -O3 -Wall -llapack -lblas # -fopenmp -funroll-loops  -fbackslash 
# LDFLAGS = -llapack

# File sorgente e oggetti
MODS = misc.o  fq_util.o  # parsing.o
OBJS = $(MODS) main.o
EXEC = optimizer 

.PHONY : clean


optimizer :  $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXEC) $^ $(LDFLAGS)

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean :
	-rm *.o *.mod $(EXEC)


