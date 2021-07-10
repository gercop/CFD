#-----------------------------------------------------------------------
# MAKEFILE for ONSS
#-----------------------------------------------------------------------
# FORTRAN FLAGS:
FFLAG  = -O3 -i8 -r8 -fpp -tpp7
#-O3 -tpp7
#-----------------------------------------------------------------------
obj     = test_onss.x

FCOMPL  = ifort
#-----------------------------------------------------------------------
# list of other directories for source files
.PREFIXES: .

.SUFFIXES:
.SUFFIXES: .f90 .o 

.f90.o:
	${FCOMPL} ${FFLAG} -c ${LIBS} $<
#-----------------------------------------------------------------------
MODULES = params.o grid.o space_integration.o \
        make_rhs.o time_integration.o whole_integration.o

FILES = main.o

OBJECTS = ${MODULES} ${FILES}
#-----------------------------------------------------------------------
APPLIC:	${OBJECTS}
	${FCOMPL} ${FFLAG} -o ${obj} ${OBJECTS}
	chmod 710 ${obj}
#-----------------------------------------------------------------------
clean:
	@rm -f ${OBJECTS} *.mod
#--END---END---END---END---END---END---END---END---END---END---END------
