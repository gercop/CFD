cccccc ********************************************************************
cccccc Doctorate Project USP/EESC: Numerical simulation of the transition *
cccccc                             for turbulence in a compressible       *
cccccc                             boundary layer on a plane plate.       *
cccccc Programmer\Guide  : Ricardo Alberto Coppola Germanos.              *            
cccccc Programmer\Guiding: Marcello Augusto Faraco de Medeiros.           *
cccccc Date      : 04/12/2007                                             *
cccccc Update    :                                                        *
cccccc Version	 : 3.7                                                    *
cccccc File      : nspar.f90                                              *
cccccc ********************************************************************
      character*07 dirgraf
      character*05 fnamegraf

      integer
     &  imax,imax_phy,jmax,kmax,tmin,tmax,qtimes,
     &  qtimes_x,qtimes_y,qtimes_z,extra_x,extra_y,extra_z

      real*8
     &  Pi,Ma,Re,Pr,M1,M2,A_1,A_2,A_3,A_4,alpha,kapa,beta,ffreq,omega,phi,
     &  x_0,y_0,z_0,Lx,Ly,Lz,Lx_spo,dx,dy,dz,dcx,dcy,dcz,dt,
     &  CFL,D,sigma,gamma,sp_x,sp_y

      !MS$DEFINE parallel        = 00

      !MS$DEFINE tp_dim          = 03

      !MS$DEFINE tp_s            = 05    !Type of simulation:
                                         !  1 = Model Equation;
                                         !  2 = Progressive Wave;
                                         !  3 = Standing Wave;
                                         !  4 = Shear Layer    - Temporal Development;
                                         !  5 = Shear Layer    - Spatial Development;
                                         !  6 = Boundary Layer - Spatial Development;
                                 
      !MS$DEFINE tp_proc         = 01    !Type of process:
                                         !  1 = Full energy equation flow; 
                                         !  2 = Isentropic flow (Adiabatic);

      !MS$DEFINE tp_viscous      = 01    !Law for viscosity:
                                         !  1 = Suntherland law; 
                                         !  2 = Dynamic Viscosity constant;

      !MS$DEFINE cancel_alongvis = 00    !Cancel alongment viscous in the y-direction

      !MS$DEFINE filter_in_x     = 01    !Compact filter in the x-direction
      !MS$DEFINE filter_in_y     = 01    !Compact filter in the y-direction
      !MS$DEFINE filter_in_z     = 01    !Compact filter in the z-direction

      !MS$DEFINE stretching_in_x = 01	 !Stretching in the x-direction
      !MS$DEFINE stretching_in_y = 01    !Stretching in the y-direction
      !MS$DEFINE stretching_in_z = 00    !Stretching in the z-direction

      !MS$DEFINE precision_plot  = 01    !Precision of data for graphics:
                                         !  1 - Single Precision;
                                         !  2 - Double Precision;

      !MS$DEFINE show_msg        = 00

      !MS$DEFINE tp_dat2plt      = 01    !Convert binary files to tecplot pattern;
      !MS$DEFINE tp_dat2mat      = 00    !Convert binary files to matlab pattern;
      !MS$DEFINE tp_dat2fft      = 00    !Carry out Fourier Spectral Analysis;

      parameter (dirgraf = 'visual/')
      !MS$IF (tp_s.EQ.2).OR.(tp_s.EQ.3).OR.(tp_s.EQ.4)
        parameter (fnamegraf = 'SN__T')
      !MS$ELSEIF (tp_s.EQ.5).OR.(tp_s.EQ.6)
        parameter (fnamegraf = 'SN__S')
      !MS$ENDIF

      parameter (
     &  qtimes        = 030,                  !Number of steps stored in time
     &  qtimes_x      = 001,                  !Number of points stored in the x coordinate
     &  qtimes_y      = 001,                  !Number of points stored in the y coordinate
     &  qtimes_z      = 001,                  !Number of points stored in the z coordinate

     &  extra_x       = 001,                  !Put together different simulations in the x coordinate
     &  extra_y       = 001,                  !Put together different simulations in the y coordinate
     &  extra_z       = 001,                  !Put together different simulations in the z coordinate

     &  imax          = 000000800+1,          !Maximum size in the x-direction
     &  jmax          = 000000200+1,          !Maximum size in the y-direction 
     &  kmax          = 000000030+1,          !Maximum size in the z-direction 
     &  tmin          = 00000000000,          !Maximum size in the time
     &  tmax          = 00000010000,          !Maximum size in the time

     &  imax_phy      = imax-000050)          !Physical domain - Points from 0..imax_phy

      parameter (
     &  Pi       = 3.141592653589793d0,       !Pi

     &  Ma       = 0.8d0,                     !Mach Number, dimensionless 
     &  Re       = 5.d+02,                    !Reynolds number, dimensionless
     &  Pr       = 1.d0,                      !Prandth number, dimensionless
     &  gamma    = 1.4d0,                     !Ratio of specific heats, dimensionless

ccccc ************************************
ccccc Parameters for mixing layer problem
ccccc ************************************
     &  M1       = Ma,                        !Mach number for y>0, dimensionless
     &  M2       = 0.0001d0,                  !Mach number for y<0, dimensionless 

     &  A_1      = 0.00100000d0,              !Wave amplitude (Fluctuation) - m
     &  A_2      = 0.00010000d0,              !Wave amplitude (Fluctuation) - m
     &  A_3      = 0.00000000d0,              !Wave amplitude (Fluctuation) - m
     &  A_4      = 0.00000000d0,              !Wave amplitude (Fluctuation) - m
     &  phi      = 0.0d0,                     !Phase difference - m
     &  x_0      = 0.0d0,                     !Initial domain in the x-direction - m
     &  y_0      = 0.0d0,                     !Initial domain in the y-direction - m
     &  z_0      = 0.0d0,                     !Initial domain in the z-direction - m
     &  Lx       = 250.000000d0-x_0,          !Comp. domain in the x-direction - m
     &  Ly       = 060.000000d0-y_0,          !Comp. domain in the y-direction - m
     &  Lz       = 015.000000d0-z_0,          !Comp. domain in the z-direction - m
     &  Lx_spo   = 600.d0,                    !Sponge domain extend from Lx to Lx_spo - m
     &  alpha    = 2.d0*Pi/Lx,                !Wavenumber in the x-direction - 1/m
     &  kapa     = 2.d0*Pi/Ly,                !Wavenumber in the y-direction - 1/m
     &  beta     = 2.d0*Pi/Lz,                !Wavenumber in the z-direction - 1/m
     &  ffreq    = 0.264d0*(M1+M2)/(2.d0*Ma), !Fundamental Frequency - 1/s
     &  omega    = 2.d0*Pi*ffreq,             !Frequency - 1/s
     &  sigma    = 2.d0,                      !Gaussian curvature

     &  sp_x     = 03.d0,                     !Stretching parameter in the x-coordinate
                                              !  For slayer - sp_x = 3.00
                                              !  For blayer - sp_x = 3.00
     &  sp_y     = 06.d0,                     !Stretching parameter in the y-coordinate
                                              !  For slayer - sp_x = 6.00
                                              !  For blayer - sp_x = 1.08
     &  CFL      = 1.0d0,                     !CFL condition - convective term
     &  D        = 1.0d0,                     !CFL condition - viscous term
     $  dt       = 1.0d-02)                   !Step size in the time     

      !MS$IF (tp_dim.EQ.1)
        parameter (
     &    dx       = Lx/dble(imax-1),         !Computational grid size in the x-direction
     &    dy       = 0.d0,                    !Computational grid size in the y-direction
     &    dz       = 0.d0,                    !Computational grid size in the z-direction
     &    dcx      = 1.d0/dble(imax-1),       !Physical grid size in the x-direction
     &    dcy      = 0.d0,                    !Physical grid size in the y-direction
     &    dcz      = 0.d0)                    !Physical grid size in the z-direction
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
