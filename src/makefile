#-----------------------------------------------------------------
#
# TMSE3DAO - Transfer matrix method for the Anderson
# model in three dimensions
#
#-----------------------------------------------------------------
 
#-----------------------------------------------------------------
#
# $Header: /home/kerekou1/rar/f77/tmODse2d/RCS/makefile.LX,v 1.3 1997/05/05 07:15:37 rar Exp $
#
# $Log: makefile.LX-PGI,v $
# Revision 1.3  1997/05/05  07:15:37  rar
# added *.for and *.l for m clean
#
# Revision 1.2  1996/12/10  11:28:09  rar
# ftn->for
#
# Revision 1.1  1996/11/14  10:40:57  rar
# Initial revision
#
#
#-----------------------------------------------------------------

F77 = pgf77

#F77FLAGS = -C -g -Mdlines -Mextend		# For debugging
F77FLAGS = -O2 -tp p6 -Munroll -Mnoframe -Mdclchk -Mextend	# For optimizing
#F77FLAGS = -xpg -C -Mprof		# For profiling

F77FLAGS2 = -O2 -tp p6 -Munroll -Mnoframe	# For the F77 lib routines

CC = gcc
CCFLAGS = -O2

LIBDIR = $(HOME)/f77/lib

CONVERT   = $(HOME)/bin/D2CDCtmp

PRINT1 = a2ps
PRINT2 = mpage -2 -A | fleet -odouble

dfcc3d:		main.o util.o inout.o random.o f77fflush.o makefile
	$(F77) $(F77FLAGS) -o dfcc3d \
		main.o util.o inout.o random.o f77fflush.o -llapack -lblas

f77fflush.o:	$(LIBDIR)/f77fflush.c
	$(CC) $(CCFLAGS) -c $(LIBDIR)/f77fflush.c

main.o:		main.f common.h
	$(F77) $(F77FLAGS) -c main.f

util.o:		util.f  common.h
	$(F77) $(F77FLAGS) -c util.f

inout.o:	inout.f common.h
	$(F77) $(F77FLAGS) -c inout.f

random.o:	random.f
	$(F77) $(F77FLAGS2) -c random.f

tar:		dfcc3d.tar
dfcc3d.tar:	makefile.* main.f inout.f util.f random.f common.h \
		dfcc3d.inp
	tar -cvf dfcc3d.tar makefile.* main.f inout.f util.f \
	random.f common.h \
	dfcc3d.inp

rcs:
	if !(-r main.f) co main.f

print:
	$(PRINT1) common.h main.f util.f inout.f random.f \
		dfcc3d.inp makefile |$(PRINT2)


clean:
	rm -f core \
		*.o *.for *.l


