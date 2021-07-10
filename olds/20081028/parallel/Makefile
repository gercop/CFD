#-----------------------------------------------------------------------
# MAKEFILE for ONSS
#-----------------------------------------------------------------------
# FORTRAN FLAGS for ifort:
#FFLAG  = -O3 -i8 -r8 -fpp -tpp7
# FORTRAN FLAGS for AURA:
#FFLAG  = -w -fpp -O3 -i8 -r8 -lmpi -lfmpi
# FORTRAN FLAGS for SOL:
#FFLAG = -w -O3 -64 -r12000 -mips4 -mp -ftpp -r8 -DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON:conform_check=ON:varargs_prototypes=ON
#FFLAG =   -64 -ftpp -r8 -DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON:conform_check=ON:varargs_prototypes=ON
#FFLAG = -w -O3 -64 -r12000 -mips4 -mp -ftpp -r8
#FFLAG = -w -O3 -64 -r12000 -mips4 -ftpp -r8
FFLAG = -r8 -O3
#FFLAG = -r8 -ftpp -DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON:conform_check=ON:varargs_prototypes=ON
#-O3 -tpp7
#-----------------------------------------------------------------------
# LIBRARIES:
#LIBS	= -I ${EAS3IOMOD_DIR} -L ${HOME}/LIB -leas3iomod_ia32r8i8 \
#        -I ${BASELIB_DIR} -L ${HOME}/LIB -lbaselib_ia32r8i8
#LIBS   = -L/home1/u21/alaible/lib/  -lios
#LIBS   = -lios -lmpi
LIBS   = -lios2 
#-----------------------------------------------------------------------
obj     = test_onss.x

# FORTRAN compiler for echigo 
#FCOMPL  = ifort
# FORTRAN compiler for aura and SOL:
#FCOMPL = f90
FCOMPL = mpif90
#-----------------------------------------------------------------------
# list of other directories for source files
.PREFIXES: .

.SUFFIXES:
.SUFFIXES: .F90 .o 

.F90.o:
	${FCOMPL} ${FFLAG} -c ${LIBS} $<
#-----------------------------------------------------------------------
MODULES = params.o output_ios.o parallel_mpi.o grid.o sutherland.o \
          interfaces.o space_integration.o \
          stresses.o make_rhs.o time_integration.o whole_integration.o 

FILES = main.o bc_free.o bc_inflow.o \
        bc_wall.o bc_outflow.o tools.o

OBJECTS = ${MODULES} ${FILES}
#-----------------------------------------------------------------------
APPLIC:	${OBJECTS}
	${FCOMPL} ${FFLAG} -o ${obj} ${OBJECTS} ${LIBS}
	chmod 710 ${obj}
#-----------------------------------------------------------------------
clean:
	@rm -f ${OBJECTS} *.mod
#--END---END---END---END---END---END---END---END---END---END---END------
