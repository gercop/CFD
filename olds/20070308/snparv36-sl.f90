cccccc *******************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition
cccccc                             for turbulence in a compressible 
cccccc                             boundary layer on a plane plate.
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos
cccccc Programmer\Guiding: Marcello A. Faraco de Medeiros
cccccc Date      :
cccccc Version	 : 3.6
cccccc File      : snparv36.f90
cccccc *******************************************************************
      character*7 dirgraf

      integer
     &  qtimes,imax_phy,imax,jmax,kmax,tmin,tmax

      real*8
     &  Pi,c,A,alpha,kapa,beta,ffreq,omega,
     &  x_0,y_0,z_0,Lx,Ly,Lz,dx,dy,dz,dcx,dcy,dcz,dt,
     &  CFL,D,sigma,gamma,R,cv,cp,Ma,Re,Pr,M1,M2

      !MS$DEFINE tp_s    = 05            !Initial condition type:
                                         !  1 = Progressive Wave;
                                         !  2 = Standing Wave;
                                         !  3 = Shear Layer - Temporal development;
                                         !  4 = Shear Layer - Spatial development - Buffer Zone;
                                         !  5 = Shear Layer - Spatial development - Stretching;
                                         !  6 = Boundary Layer - Spatial development;
                                 
      !MS$DEFINE tp_proc = 01            !Process type:
                                         !  1 = Full energy equation flow; 
                                         !  2 = Isentropic flow (Adiabatic);

      !MS$DEFINE cancel_alongvis = 0     !Cancel alongment viscous in y-direction

      !MS$DEFINE filter_in_x     = 1     !Compact filter in x-direction
      !MS$DEFINE filter_in_y     = 1     !Compact filter in y-direction
      !MS$DEFINE filter_in_z     = 0     !Compact filter in z-direction

      !MS$DEFINE stretching_in_x = 0	 !Stretching in x-direction
      !MS$DEFINE stretching_in_y = 1     !Stretching in y-direction
      !MS$DEFINE stretching_in_z = 0     !Stretching in z-direction

      parameter (dirgraf = 'visual/')

      parameter (
     &  qtimes        = 500,
     &  imax          = 000000100+1, !Maximum size in x-direction 
     &  jmax          = 000000100+1, !Maximum size in y-direction 
     &  kmax          = 000000000+1, !Maximum size in z-direction 
     &  tmin          = 00000013000, !Maximum size in time
     &  tmax          = 01000000000, !Maximum size in time
     &  imax_phy      = imax-0100) !Physical domain - Points from 0..imax_phy

      parameter (
     &  Pi       = 3.141592653589793d0,       !Pi
     &  Ma       = 0.5d0,                     !Mach Number 
     &  Re       = 1.d+02,                    !Reynolds number    $$  
      \left [
        \begin{array}{c}
          \frac{1}{12dx}(f_{M-1} + 28 f_{M} - 28 f_{1} - f_{2})       \\
          \frac{1}{12dx}(f_{M} + 28 f_{0} - 28 f_{2} - f_{3})         \\
          \vdots                                                      \\
          \frac{1}{12dx}(f_{i+2} + 28 f_{i+1} - 28 f_{i-1} - f_{i-2}) \\
          \vdots                                                      \\
          \frac{1}{12dx}(f_{M-3} + 28 f_{M-2} - 28 f_{M} - f_{0})     \\
          \frac{1}{12dx}(f_{M-2} + 28 f_{M-1} - 28 f_{0} - f_{1})     \\
        \end{array}
    \right],
    \label{esquema_6ordem_comp_2derivada}
  $$

     &  Pr       = 0.71d0,                    !Prandth number
     &  c        = 340.210683840878d0,        !Speed of sound - m/s
     &  gamma    = 1.4d0,                     !Ratio of specific heats
     &  R        = 287.d0,                    !Gas constant - m^2/s^2
     &  cv       = R/(gamma-1.d0),            !Spec. heat at const. volume - J/(Kg K)
     &  cp       = gamma*R/(gamma-1.d0),      !Spec. heat at const. pressure - J/(Kg K)

ccccc ***********************************
ccccc Parameters for shear layer problem
ccccc ***********************************
     &  M1       = Ma,                        !Mach number for y>0
     &  M2       = Ma-0.1d0,                  !Mach number for y<0 

     &  A        = 0.10000000d0,              !Wave amplitude (Fluctuation) - m
     &  x_0      = 0.0d0,                     !Initial domain in x direction - m
     &  y_0      = 0.0d0,                     !Initial domain in y direction - m
     &  z_0      = 0.0d0,                     !Initial domain in z direction - m
     &  Lx       = 100.000000d0-x_0,          !Comp. domain in x-direction - m
     &  Ly       = 050.000000d0-y_0,          !Comp. domain in y-direction - m
     &  Lz       = 016.000000d0-z_0,          !Comp. domain in z-direction - m
     &  alpha    = 4.d0*Pi/Lx,                !Wave number in x-direction - 1/m
     &  kapa     = 2.d0*Pi/Ly,                !Wave number in y-direction - 1/m
     &  beta     = 2.d0*Pi/Lz,                !Wave number in z-direction - 1/m
     &  ffreq    = 0.132d0*(M1+M2)/(2.d0*Ma), !Fundamental Frequency - 1/s
     &  omega    = 2.d0*Pi*ffreq,             !Frequency - 1/s
     &  sigma    = 2.d0,                      !Gaussian curvature

     &  CFL      = 1.0d0,                     !CFL condition - convective term
     &  D        = 1.0d0,                     !CFL condition - viscous term
     &  dx       = Lx/dble(imax-1),           !Step size in the x-direction
     &  dy       = Ly/dble(jmax-1),           !Step size in the y-direction
     &  dz       = 0.d0, !Lz/dble(kmax-1),           !Step size in the z-direction
     &  dcx      = 1.d0/dble(imax-1),         !Step size in the x-direction
     &  dcy      = 1.d0/dble(jmax-1),         !Step size in the y-direction
     &  dcz      = 0.d0, !1.d0/dble(kmax-1),         !Step size in the z-direction
     &  dt       = 1.0d-02)                   !Step size in the time     

cccc *********************************************************************
cccc *********************************************************************
cccc Rules: 
cccc   Definition of Variables
cccc     |-> se      : Exact solution
cccc     |-> u       : Velocity component in x-direction 
cccc     |-> v       : Velocity component in y-direction
cccc     |-> w       : Velocity component in z-direction  
cccc     |-> p       : Pressure 
cccc     |-> rho     : Density 
cccc     |-> Et      : Total energy
cccc     |-> e       : Internal energy
cccc     |-> Ma      : Mach number
cccc     |-> Re      : Reynolds number
cccc     |-> x       : Reference point in x-direction
cccc     |-> y       : Reference point in y-direction
cccc     |-> z       : Reference point in z-direction
cccc     |-> t       : Reference point in time
cccc     |-> dx      : Step size in x-direction
cccc     |-> dy      : Step size in y-direction
cccc     |-> dz      : Step size in z-direction
cccc     |-> dt      : Step size in time
cccc     |-> imax    : Maximum position of comput. domain in x-direction
cccc     |-> jmax    : Maximum position of comput. domain in y-direction
cccc     |-> kmax    : Maximum position of comput. domain in z-direction
cccc     |-> tmax    : Maximum position of the mesh in time
cccc ---------------------------------------------------------------------
cccc Literal operators: 
cccc     |-> // catenation
cccc Arithmetic operators:
cccc     |-> +  Sum
cccc     |-> -  Subtraction
cccc     |-> *  Multiplication
cccc     |-> /  Division
cccc     |-> ** Potentiation
cccc Relate operators:
cccc     |-> .LT.  Less than
cccc     |-> .LE.  Less and equal than
cccc     |-> .EQ.  Equal
cccc     |-> .NE.  Not Equal 
cccc     |-> .GT.  Great than
cccc     |-> .GE.  Great and equal than
cccc Intrisic function: 
cccc     |-> Int(x)     : Convert the variable x for a integer type
cccc     |-> Real(x)    : Convert the variable x for a real type
cccc     |-> Dble(x)    : Convert the variable x for a double precision type.
cccc     |-> char(x)    : Convert the variable x integer for charater
cccc     |-> IFix(x)    : Convert the variable x real for integer
cccc Convertions - General:
cccc     |-> Char       : Create character array (string).
cccc     |-> Double     : Convert string to numeric character codes.
cccc     |-> Cellstr    : Create cell array of strings from character array.
cccc     |-> Blanks     : String of blanks.
cccc     |-> Deblank    : Remove trailing blanks.
cccc     |-> Eval       : Execute string with MATLAB expression.
cccc Convertions - String tests:
cccc     |-> Ischar     : True for character array (string).
cccc     |-> Iscellstr  : True for cell array of strings.
cccc     |-> Isletter   : True for letters of the alphabet.
cccc     |-> Isspace    : True for white space characters. 
cccc Convertions - String operations:
cccc     |-> Strcat     : Concatenate strings.
cccc     |-> Strvcat    : Vertically concatenate strings.
cccc     |-> Strcmp     : Compare strings.
cccc     |-> Strncmp    : Compare first N characters of strings.
cccc     |-> Strcmpi    : Compare strings ignoring case.
cccc     |-> Strncmpi   : Compare first N characters of strings 
cccc                      ignoring case.
cccc     |-> Findstr    : Find one string within another.
cccc     |-> Strjust    : Justify character array.
cccc     |-> Strmatch   : Find possible matches for string.
cccc     |-> Strrep     : Replace string with another.
cccc     |-> Strtok     : Find token in string.
cccc     |-> Upper      : Convert string to uppercase.
cccc     |-> Lower      : Convert string to lowercase.
cccc Convertions - String to number conversion:
cccc     |-> Num2str    : Convert number to string.
cccc     |-> Int2str    : Convert integer to string.
cccc     |-> Mat2str    : Convert matrix to eval'able string.
cccc     |-> Str2double : Convert string to double precision value.
cccc     |-> Str2num    : Convert string matrix to numeric array.
cccc     |-> Sprintf    : Write formatted data to string.
cccc     |-> Sscanf     : Read string under format control.
cccc Convertions - Base number conversion:
cccc     |-> Hex2num    : Convert IEEE hexadecimal to double precision number.
cccc     |-> Hex2dec    : Convert hexadecimal string to decimal integer.
cccc     |-> Dec2hex    : Convert decimal integer to hexadecimal string.
cccc     |-> Bin2dec    : Convert binary string to decimal integer.
cccc     |-> Dec2bin    : Convert decimal integer to binary string.
cccc     |-> Base2dec   : Convert base B string to decimal integer.
cccc     |-> Dec2base   : Convert decimal integer to base B string.
cccc Blocks of commands:
cccc     |-> select case (variable)
cccc     |->   case (condition 1)
cccc     |->   case (condition 2)
cccc     |-> end select
cccc ---------------------------------------------------------------------
cccc Rule of implementation: 
cccc     |-> Evitar comandos "IF" dentro de blocos de comandos "DO"
cccc *********************************************************************
cccc *********************************************************************
