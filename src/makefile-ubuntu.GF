#-----------------------------------------------------------------
#
# tmseGR
#
#-----------------------------------------------------------------
#
# $Header: /home/cvs/phsht/GrapheneTMM/src/Restart/ALL/makefileGR,v 1.1 2011/07/22 17:49:19 ccspam Exp $
#
# $Log: makefileGR,v $
# Revision 1.1  2011/07/22 17:49:19  ccspam
# Programs for ZZ and AC with the restart routine for francesca
#
# Revision 1.1  2011/05/31 13:55:07  ccspam
# *** empty log message ***
#
# Revision 1.1  2011/05/06 08:13:09  phsht
# 1st installement
#
#-----------------------------------------------------------------#

F90 = gfortran -std=legacy
F77 = gfortran -ffixed-form -std=legacy

F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer # debug
#F90FLAGS = -pg -fbounds-check -fno-range-check -fimplicit-none # profiling
F77FLAGS = -fbounds-check -fno-range-check # debug

LIBDIR = $(HOME)/f77/lib
LIBFLAGS  =  -lm #-llapack -lblas

OBJFILES = tm1dNNN_modules.o random.o main.o util.o inout.o 

tm1dNNN.GF: $(OBJFILES)
	$(F90) $(F90FLAGS) -o $@ \
		$(OBJFILES) $(LIBFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -DF90 -c $< 

#tm1dNNN_modules.o:     tm1dNNN_modules.f90
#	$(F90) $(F90FLAGS) -c tm1dNNN_modules.f90

#random.o:     random.f90
#	$(F90) $(F90FLAGS) -c random.f90

#main.o:     main.f90
#	$(F90) $(F90FLAGS) -c main.f90

#util.o:     util.f90
#	$(F90) $(F90FLAGS) -c util.f90

#inout.o:     inout.f90
#	$(F90) $(F90FLAGS) -c inout.f90

clean:
	rm -f core *.mod *.o




