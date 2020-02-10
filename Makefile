
FC = gfortran 
#FC = ifort
FLAGS = -O3
#FC = ifort 
#FLAGS = -g -CB -check all 
 
TARGET1 = md 

SOURCES1 = constants.f90 \
	   parameters.f90 \
	   list.f90 \
	   boxes.f90 \
	   input.f90 \
	   forces.f90 \
	   simulations.f90 \
	   dynamics.f90 \
	   clock.f90 \
	   md.f90 
OBJS1 = $(SOURCES1:.f90=.o)

%.o: %.f90 
	$(FC) $(FLAGS) -c  $<

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(FC) -o $(TARGET1) $(OBJS1) $(LIBS)


test: test.o
	$(FC) $(FLAGS) -c test.f90
	$(FC) -o test  test.o list.o clock.o



clean:
	rm *.o *.mod #$(TARGET1) 

parameters.o : constants.o
boxes.o : constants.o
simulations.o : constants.o parameters.o forces.o dynamics.o clock.o
forces.o : constants.o list.o boxes.o 	
dynamics.o : constants.o parameters.o	
md.o : constants.o parameters.o list.o boxes.o input.o simulations.o
