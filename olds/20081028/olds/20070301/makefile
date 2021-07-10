#Commands: scp *.f90 lammpi@143.107.225.110:/home/lammpi/version36
#          ssh lammpi@143.107.225.110  
#          ./resetdir; ./resetfile; make shear; mpirun C prg;
#          scp -r lammpi@143.107.225.110:/home/lammpi/version36/visual /mnt/dados/usp/temp/version36/
#          rm -rf dat2plt *.plt ; ifort dat2plt.f -o dat2plt /usr/local/tecplot8/lib/tecio.a /usr/local/tecplot8/lib/ctype.o ; dat2plt ;

#DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON:conform_check=ON:varargs_prototypes=ON

#------------------------------------------------------------------------------------------------------------------------------
FFLAG = -r8 -O3 -fixed -132
#------------------------------------------------------------------------------------------------------------------------------
#FCOMPL = mpif90
FCOMPL = ifort
#------------------------------------------------------------------------------------------------------------------------------
#LIBS  = /usr/local/tecplot8/lib/tecio.a /usr/local/tecplot8/lib/ctype.o -lios3
#------------------------------------------------------------------------------------------------------------------------------
obj   = prg
#------------------------------------------------------------------------------------------------------------------------------
# List of other directories for source files
.PREFIXES: .

.SUFFIXES:
.SUFFIXES: .f90 .o 

.f90.o:
	${FCOMPL} ${FFLAG} -c $<
#------------------------------------------------------------------------------------------------------------------------------
MODULES = snmpi.o snmethod.o sninitial.o snbound.o snfilter.o snconf.o 

FILES   = snttcl.o sngeneral.o          

OBJECTS = ${MODULES} ${FILES}
#------------------------------------------------------------------------------------------------------------------------------
shear:	${OBJECTS}
	${FCOMPL} ${FFLAG} -o ${obj} ${OBJECTS} ${LIBS}
	chmod 710 ${obj}
#------------------------------------------------------------------------------------------------------------------------------
clean:
	@rm -f ${OBJECTS} *.mod
	@rm -f ${OBJECTS} *.o
	@rm -f ${OBJECTS} *~
#------------------------------------------------------------------------------------------------------------------------------
