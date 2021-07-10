cccccc ********************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition *
cccccc                             for turbulence in a compressible       *
cccccc                             boundary layer on a plane plate.       *
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos.              *            
cccccc Programmer\Guiding: Marcello Augusto Faraco de Medeiros.           *
cccccc Date      :                                                        *
cccccc Update    :                                                        *
cccccc Version	 : 3.7                                                    *
cccccc File      : snparv.f90                                             *
cccccc ********************************************************************
      character*7 dirgraf

      integer
     &  qtimes,qtimes_x,qtimes_y,qtimes_z,imax,imax_phy,jmax,kmax,tmin,tmax

      real*8
     &  Pi,Ma,Re,Pr,M1,M2,c,A_1,A_2,A_3,A_4,alpha,kapa,beta,ffreq,omega,
     &  x_0,y_0,z_0,Lx,Ly,Lz,Lx_spo,dx,dy,dz,dcx,dcy,dcz,dt,CFL,D,sigma,gamma

      !MS$DEFINE parallel        = 01

      !MS$DEFINE tp_dim          = 03

      !MS$DEFINE tp_s            = 05    !Simulation type:
                                         !  1 = Model Equation;
                                         !  2 = Progressive Wave;
                                         !  3 = Standing Wave;
                                         !  4 = Shear Layer    - Temporal Development;
                                         !  5 = Shear Layer    - Spatial Development;
                                         !  6 = Boundary Layer - Spatial Development;
                                 
      !MS$DEFINE tp_proc         = 01    !Process type:
                                         !  1 = Full energy equation flow; 
                                         !  2 = Isentropic flow (Adiabatic);

      !MS$DEFINE tp_viscous      = 01    !Viscous type:
                                         !  1 = Suntherland law; 
                                         !  2 = Dynamic Viscosity constant;

      !MS$DEFINE cancel_alongvis = 00    !Cancel alongment viscous in the y-direction

      !MS$DEFINE filter_in_x     = 01    !Compact filter in the x-direction
      !MS$DEFINE filter_in_y     = 01    !Compact filter in the y-direction
      !MS$DEFINE filter_in_z     = 01    !Compact filter in the z-direction

      !MS$DEFINE stretching_in_x = 00	 !Stretching in the x-direction
      !MS$DEFINE stretching_in_y = 01    !Stretching in the y-direction
      !MS$DEFINE stretching_in_z = 00    !Stretching in the z-direction

      !MS$DEFINE show_msg        = 00

      parameter (dirgraf = 'visual/')

      parameter (
     &  qtimes        = 500,
     &  qtimes_x      = 2,
     &  qtimes_y      = 1,
     &  qtimes_z      = 1,

     &  imax          = 00000000081, !Maximum size in the x-direction 
     &  jmax          = 00000000081, !Maximum size in the y-direction 
     &  kmax          = 00000000081, !Maximum size in the z-direction 
     &  tmin          = 00000000000, !Maximum size in the time
     &  tmax          = 00100000000, !Maximum size in the time

     &  imax_phy      = imax-000080) !Physical domain - Points from 0..imax_phy

      parameter (
     &  Pi       = 3.141592653589793d0,       !Pi

     &  Ma       = 0.8d0,                     !Mach Number, dimensionless 
     &  Re       = 5.d+02,                    !Reynolds number, dimensionless
     &  Pr       = 0.71d0,                    !Prandth number, dimensionless
     &  c        = 340.210683840878d0,        !Speed of sound - m/s
     &  gamma    = 1.4d0,                     !Ratio of specific heats, dimensionless

ccccc ************************************
ccccc Parameters for mixing layer problem
ccccc ************************************
     &  M1       = Ma,                        !Mach number for y>0, dimensionless
     &  M2       = 0.001d0,                   !Mach number for y<0, dimensionless 

     &  A_1      = 0.10000000d0,              !Wave amplitude (Fluctuation) - m
     &  A_2      = 0.01000000d0,              !Wave amplitude (Fluctuation) - m
     &  A_3      = 0.00000000d0,              !Wave amplitude (Fluctuation) - m
     &  A_4      = 0.00000000d0,              !Wave amplitude (Fluctuation) - m
     &  x_0      = 0.0d0,                     !Initial domain in the x-direction - m
     &  y_0      = 0.0d0,                     !Initial domain in the y-direction - m
     &  z_0      = 0.0d0,                     !Initial domain in the z-direction - m
     &  Lx       = 050.000000d0-x_0,          !Comp. domain in the x-direction - m
     &  Ly       = 060.000000d0-y_0,          !Comp. domain in the y-direction - m
     &  Lz       = 025.000000d0-z_0,          !Comp. domain in the z-direction - m
     &  Lx_spo   = 400.d0,                    !Sponge domain extend from Lx to Lx_spo - m
     &  alpha    = 4.d0*Pi/Lx,                !Wavenumber in the x-direction - 1/m
     &  kapa     = 2.d0*Pi/Ly,                !Wavenumber in the y-direction - 1/m
     &  beta     = 4.d0*Pi/Lz,                !Wavenumber in the z-direction - 1/m
     &  ffreq    = 0.186d0*(M1+M2)/(2.d0*Ma), !Fundamental Frequency - 1/s
     &  omega    = 2.d0*Pi*ffreq,             !Frequency - 1/s
     &  sigma    = 2.d0,                      !Gaussian curvature
     &  CFL      = 1.0d0,                     !CFL condition - convective term
     &  D        = 1.0d0,                     !CFL condition - viscous term
     &  dt       = 1.0d-02)                   !Step size in the time     

      !MS$IF (tp_dim.EQ.1)
        parameter (
     &    dx       = Lx/dble(imax-1),         !Step size in the x-direction
     &    dy       = 0.d0,                    !Step size in the y-direction
     &    dz       = 0.d0,                    !Step size in the z-direction
     &    dcx      = 1.d0/dble(imax-1),       !Step size in the x-direction
     &    dcy      = 0.d0,                    !Step size in the y-direction
     &    dcz      = 0.d0)                    !Step size in the z-direction
      !MS$ELSEIF (tp_dim.EQ.2)
        parameter (
     &    dx       = Lx/dble(imax-1),         
     &    dy       = Ly/dble(jmax-1),         
     &    dz       = 0.d0,                    
     &    dcx      = 1.d0/dble(imax-1),       
     &    dcy      = 1.d0/dble(jmax-1),       
     &    dcz      = 0.d0)                    
      !MS$ELSEIF (tp_dim.EQ.3)
        parameter (
     &    dx       = Lx/dble(imax-1),           
     &    dy       = Ly/dble(jmax-1),           
     &    dz       = Lz/dble(kmax-1),           
     &    dcx      = 1.d0/dble(imax-1),         
     &    dcy      = 1.d0/dble(jmax-1),         
     &    dcz      = 1.d0/dble(kmax-1))         
      !MS$ENDIF

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
