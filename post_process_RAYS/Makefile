# Find all source files, create a list of corresponding object files
F90_SRCS=$(wildcard *.f90)
F90_OBJS=$(patsubst %.f90,%.o,$(F90_SRCS))

F77_SRCS=$(wildcard *.f)
F77_OBJS=$(patsubst %.f,%.o,$(F77_SRCS))


# Ditto for mods (They will be in both lists)
MODS=$(wildcard *_m.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

# Compiler/Linker settings
FC = gfortran
FLFLAGS = -g
FCFLAGS = -g -c -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
PROGRAM = post_process_RAYS
PRG_OBJ = $(PROGRAM).o
		
# make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(F90_OBJS) : %.o : %.f90
	$(FC) $(FCFLAGS) -o $@ $<

$(F77_OBJS) : %.o : %.f
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(F90_OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

depend .depend:
	makedepf90  -o $(PROGRAM) *.f90 *.f > .depend
	
-include .depend

clean:
	rm $(F90_OBJS)
	rm $(F77_OBJS)
	rm $(MOD_OBJS)
	rm *.mod

distclean: clean
	rm .depend
	
# If something doesn't work right, have a 'make debug' to 
# show what each variable contains.
debug:
	@echo "F90_SRCS = $(F90_SRCS)"
	@echo "F90_OBJS = $(F90_OBJS)"
	@echo "F77_SRCS = $(F77_SRCS)"
	@echo "F77_OBJS = $(F77_OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"