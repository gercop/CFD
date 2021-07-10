#Commands: scp *.f90 lammpi@143.107.225.110:/home/lammpi/version36
#          ssh lammpi@143.107.225.110  
#          ./resetdir; ./resetfile; make shear; mpirun C prg;
#          scp *.dat gercop@143.107.236.123:/mnt/dados/usp/temp/version36/visual
#          scp -r lammpi@143.107.225.110:/home/lammpi/version36/visual /mnt/dados/usp/temp/version36/
#          rm -rf dat2plt *.plt ; ifort dat2plt.f -o dat2plt /usr/local/tecplot8/lib/tecio.a /usr/local/tecplot8/lib/ctype.o ; dat2plt ;
#          ls visual ; resetdir ; make dat2plt ; prg ;

#DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON:conform_check=ON:varargs_prototypes=ON

#------------------------------------------------------------------------------------------------------------------------------
FFLAG = -O3 -fixed -132
#------------------------------------------------------------------------------------------------------------------------------
FCOMPL = mpif90 -q
#------------------------------------------------------------------------------------------------------------------------------
LIBS01  = 
LIBS10  = 
LIBS11  = 
LIBS50  = /usr/local/tecplot10/lib/tecio.a /usr/local/tecplot10/lib/ctype.o
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
MODULES01 = nsmpi.o nsmethod.o
MODULES10 = nsmpi.o nsmethod.o nsexactly.o nsinitial.o nsfilter.o nsconf.o
MODULES11 = nsmpi.o nsmethod.o nsinitial.o nsfilter.o nsconf.o
MODULES50 = nsmethod.o nsdat2asc.o nsdat2plt.o nsdat2mat.o nsdat2fft.o

FILES01   = nsprecess.o          
FILES10   = nsmain.o nsgeneral.o          
FILES11   = nsmain.o nsgeneral.o                    
FILES50   = nsgeneral.o nsprocess.o 

OBJECTS01 = ${MODULES01} ${FILES01}
OBJECTS10 = ${MODULES10} ${FILES10}
OBJECTS11 = ${MODULES11} ${FILES11}
OBJECTS50 = ${MODULES50} ${FILES50}
#------------------------------------------------------------------------------------------------------------------------------
grid:	${OBJECTS01}
	${FCOMPL} ${FFLAG} -o ${obj} ${OBJECTS01} ${LIBS01}
	chmod 710 ${obj}
#------------------------------------------------------------------------------------------------------------------------------
acoustic:${OBJECTS10}
	${FCOMPL} ${FFLAG} -o ${obj} ${OBJECTS10} ${LIBS10}
	chmod 710 ${obj}
#------------------------------------------------------------------------------------------------------------------------------
shear:	${OBJECTS11}
	${FCOMPL} ${FFLAG} -o ${obj} ${OBJECTS11} ${LIBS11}
	chmod 710 ${obj}
#------------------------------------------------------------------------------------------------------------------------------
afterpro:${OBJECTS50}
	 ${FCOMPL} ${FFLAG} -o ${obj} ${OBJECTS50} ${LIBS50}
	 chmod 710 ${obj}
#------------------------------------------------------------------------------------------------------------------------------
clean:
	@rm -f ${OBJECTS} *.mod
	@rm -f ${OBJECTS} *.o
	@rm -f ${OBJECTS} *~
#------------------------------------------------------------------------------------------------------------------------------
