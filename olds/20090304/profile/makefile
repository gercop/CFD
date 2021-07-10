#-----------------------------------------------------------------------
# MAKEFILE for PC
#
# last updated: 02/01/05 Dominic von Terzi
#
#-----------------------------------------------------------------------
# FORTRAN flags:
FFLAGC = -i8 -r8
#-----------------------------------------------------------------------
# Libraries:
#LIBDIR  = ${HOME}/LIB
LIBDIR  = .
LIBS    = -L${LIBDIR} 
#-----------------------------------------------------------------------
opt     = 
profile = 
obj     = profcom
ildir   =  
input   = 

FCOMPL  = ifort 
#-----------------------------------------------------------------------
# list of other directories for source files
.PREFIXES: .

.SUFFIXES:
.SUFFIXES: .f .c .o .fil

.f.o:
	$(FCOMPL) -c $(ildir) $(FFLAGC) $(opt) $(profile) $<

.f.fil:
	$(FCOMPL) -il $(FFLAGC) $<
#-----------------------------------------------------------------------
FOBJECTS = profile.o
#-----------------------------------------------------------------------
OBJECTS = $(FOBJECTS)
INLINE  = 
#-----------------------------------------------------------------------
APPLIC:	$(INLINE) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) $(profile) -o $(obj) $(OBJECTS) $(LIBS)
	chmod 710 $(obj)
#-----------------------------------------------------------------------
clean:
	@rm -f $(INLINE) $(OBJECTS)

fclean:
	@rm -f $(INLINE) $(FOBJECTS)
#--END---END---END---END---END---END---END---END---END---END---END------
