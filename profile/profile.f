c     03/27/90
c     Diese Version des PROFKOM Programmes ist im maerz 1990 per
c     Direktuebertragung in Tucson eingetroffen.
c     So sollte diese Version endlich die Faehigkeiten der alten
c     Programme vereinigen. D.h. sowohl die korrekte v-Verteilung
c     als auch die richtigen Grenzschichtparameter ausgeben.
c
c     04/12/90
c     Installation auf SUN mit anschluss an das IOS ein- ausgabesystem.
c     Variable HFELD entfernt
c
c     01/16/90
c     upgraded to real*8
c
      PROGRAM PROFN
C
C     GRUNDSTROEMUNGSPROFILE FUER DIE EBENE KOMPRESSIBLE
C     PLATTENGRENZSCHICHT (MACH=0-10)
C
C     ******************************************************************
C     *                                                                *
C     *  PARAMETERLISTE:                                               *
C     *                                                                *
C     *  ICPPRV=0    CP KONSTANT,PRANDTLZAHL KONSTANT                  *
C     *  ICPPRV=1    CP KONSTANT,PRANDTLZAHL VARIABEL                  *
C     *  ICPPRV=2    CP VARIABEL,PRANDTLZAHL VARIABEL                  *
C     *                                                                *
C     *  LSTR=0      SCHRITTWEITE KONSTANT                             *
C     *  LSTR=1      STRECKUNG ASINH                                   *
C     *  LSTR=2      STUTZSTELLEN NACH CHEBYSHEV                       *
C     *                                                                *
C     *  ICOM=0      ADIABATE WAND                                     *
C     *  ICOM=1      KOSTANTE WANDTEMPERATUR                           *
C     *                                                                *
C     *  RM1         MACHZAHL DER AUSSENSTROEMUNG                      *
C     *  T1          TEMPERATUR DER AUSSENSTROEMUNG                    *
C     *  CP          SPEZIFISCHE WAERMEKAPAZITAET BEI KONSTANTEM DRUCK *
C     *  CV          SPEZIFISCHE WAERMEKAPAZITAET BEI KONSTANTEM       *
C     *              VOLUMEN                                           *
C     *  DN          SCHRITTWEITE BEI LSTR=0                           *
C     *  ETA         MAXIMALER WERT DER AEHNLICHKEITSKOORDINATE ETA    *
C     *              (BEGRENZT DAS INTEGRATIONSGEBIET)                 *
C     *  ETAS        STELLE, AB DER SCHRITTWEITE KONSTANT              *
C     *  ASTR        STRECKFAKTOR FUER STRECKUNG ASINH                 *
C     *  YSTR        STRECKFAKTOR CHEBYSCHEFF                          *
C     *  NSTR        STRECKFAKTOR CHEBYSCHEFF                          *
C     *  DUDHS       SCHAETZWERT FUER 1.ABLEITUNG VON U NACH ETA       *
C     *  THETA       SCHAETZWERT FUER THETA (-> MACK)                  *
C     *              (NUR FUER RECHNUNG MIT ADIABATER WAND; KANN       *
C     *              NAHEZU UNVERAENDERT GELASSEN WERDEN)              *
C     *  GSCH        SCHAETZWERT FUER 1.ABLEITUNG VON THETA NACH ETA   *
C     *              (NUR FUER RECHNUNG MIT KONSTANTER WANDTEMPERATUR; *
C     *              KANN NAHEZU UNVERAENDERT GELASSEN WERDEN)         *
C     *  TW          WANDTEMPERATUR BEI RECHNUNG MIT KONSTANTER        *
C     *              WANDTEMPERATUR                                    *
C     *  PR0         PRANDTLZAHL FUER RECHNUNG BEI ICPPRV=0            *
C     *  IAUT=1      PROGRAMM WAEHLT SELBSTSTAENDIG SCHAETZWERTE       *
C     *  IAUT=0      SCHAETZWERTE MUESSEN IM INPUT ANGEGEBEN WERDEN    *
C     *              (AB MACH 6)                                       *
C     *                                                                *
C     *  GENAUIGKEIT EPS MUSS IN DER PARAMETERANWEISUNG EINGEG. WERDEN.*
C     *  SCHRITTWEITE DN WIRD BEI RECHNUNG MIT KONSTANTER SHRITTWEITE  *
C     *  FUER DN GROESSER 0.1 VOM PROGRAMM VERKLEINERT,ANSCHLIESSEND   *
C     *  ABER NUR DIE IM INPUT BESTIMMTEN STELLEN AUSGEDRUCKT.         *
C     *  DUDHS UND THETAS SIND SCHAETZWERTE MIT DENEN DIE ITERATION    *
C     *  BEGINNT. BEI HOHER MACHZAHL (AB MACH=2) SIND FUER DIE BE-     *
C     *  RECHNUNG GENAUERE SCHAETZWERTE ERFORDERLICH. DIESE            *
C     *  SCHAETZWERTE SIND BIS MACH 6.0 IN EINEM DATA-BLOCK IN         *
C     *  ZEILE     DES HAUPTPROGRAMMS MIT DEN ZUGEHOERIGEN MACH-       *
C     *  ZAHLEN VORGEGEBEN. WIRD DER INPUT-PARAMETER IAUT=1 GEWAEHLT,  *
C     *  (SINNVOLL BIS CA. MACH 6), SO WERDEN DIESE SCHAETZWERTE       *
C     *  AUTOMATISCH VOM PROGRAMM GEWAEHLT.                            *
C     *                                                                *
C     *  LSTR DIENT ZUR FESTLEGUNG OB MIT KONSTANTER SCHRITTWEITE GE-  *
C     *  ARBEITET WERDEN SOLL ODER NICHT.                              *
C     *  NSTR LEGT DIE ANZAHL DER PUNKTE FUER DIE CHEBYSHEV'SCHE       *
C     *  KOLLOKATIONSVERTEILUNG FEST.                                  *
C     *  YSTR IST DER STRECKFAKTOR FUER RECHNUNG MIT CHEBYSHEV POLYN.  *
C     *  SCHRITTWEITE DER PROFILE MUSS MIT NACHFOLGENDEM PROGRAMM      *
C     *  UEBEREINSTIMMEN, EBENSO YSTR,NSTR,T1(ANSTROEMTEMP.)  UND      *
C     *  RM1 (MACHZAHL)                                                *
C     *  ICOM LEGT FEST OB MIT ADIABATER WAND ODER MIT KONSTANTER      *
C     *  TEMPERATUR AN DER WAND (TW,ICOM=1) GERECHNET WERDEN SOLL      *
C     *  DIE BERECHNUNG FUER MACHZAHL 0 (INKOMPRESSIBLER FALL) MUSS    *
C     *  FUER ADIABATE WAND DURCHGEFUEHRT WERDEN                       *
C     *  BEI RECHNUNG MIT ADIABATER WAND SIND SCHAETZWERTE FUER        *
C     *  THETA(1)(THETAS) UND DUDH(1)(DUDHS) VORZUGEBEN,BEI RECHNUNG   *
C     *  MIT KONSTANTER WANDTEMPERATUR MUESSEN TW UND SCHAETZWERTE     *
C     *  FUER DUDH(1) UND G(1)(ABLEITUNG VON THETA) VORGEGEBEN WERDEN  *
C     *  UNBEKANNTE SCHAETZWERTE FUER HOEHERE MACHZAHLEN KOENNEN DURCH *
C     *  VORAUSGEHENDE RECHNUNGEN MIT NIEDEREN MACHZAHLEN ERMITTELT    *
C     *  WERDEN. GRUNDSAETZLICH SOLLTE DER SCHAETZWERT  FUER DUDH      *
C     *  IMMER ETWAS KLEINER ALS DER TATSAECHLICHE WERT GEWAEHLT       *
C     *  WERDEN.                                                       *
C     *  DHMAX LEGT MAXIMALE SCHITTWEITE FUER RECHNUNG MIT NICHTAEQU.  *
C     *  GITTER FEST. FALLS KEINE KONVERGENZ ERZIELT WIRD IST DIESE    *
C     *  U.U. ZU VARIIEREN.                                            *
C     *  BEI RECHNUNG MIT NICHAEQUIDISTANTEM GITTER UND ICPPRV=2       *
C     *  KANN EVENTUELL KEINE KONVERGENZ ERZIELT WERDEN ODER DAS       *
C     *  PROGRAMM STUERZT AB, DA WERTE DES MITTLEREN CP DES VORIGEN    *
C     *  ITERATIONSSCHRITTS MITBENUTZT WERDEN. IN DIESEM FALL MUSS     *
C     *  LSTR=0 GEWAEHLT WERDEN.                                       *
C     *                                                                *
C     *  AUSGABEPARAMETER :                                            *
C     *                                                                *
C     *       U = DIMENSIONSLOSE LAENGSGESCHWINDIGKEIT                 *
C     *       DUDH = 1. ABLEITUNG VON U NACH ETA                       *
C     *       D2UDH = 2.ABLEITUNG VON U NACH ETA                       *
C     *       TDL = DIMENSIONSLOSE GESCHWINDIGKEIT (Tempetrure ???)    *
C     *       DTDH = 1.ABLEITUNG VON TDL NACH ETA                      *
C     *       D2TDH = 2.ABLEITUNG VON TDL NACH ETA                     *
C     *       FV = GROESSE ZUR BERECHNUNG VON V                        *
C     *       RHODL = DIMENSIONSLOSE DICHTE                            *
C     *       DRHODH = 1.ABLEITUNG VON RHODL NACH ETA                  *
C     *       D2RDH = 2.ABLEITUNG VON RHODL NACH ETA                   *
C     *       DRMUE = DIMENSIONSLOSE VISKOSITAET                       *
C     *       DMUEDT = 1.ABLEITUNG VON DRMUE NACH TDL                  *
C     *       D2MDT = 2.ABLEITUNG VON DRMUE NACH TDL                   *
C     *       DRK = DIMENSIONSLOSE WAERMELEITFAEHIGKEIT                *
C     *       DRKDT = 1.ABLEITUNG VON DRK NACH TDL                     *
C     *       D2RKDT = 2.ABLEITUNG VON DRK NACH TDL                    *
C     *       PR = PRANDTLZAHL                                         *
C     *       DRD2U = VERALLGEMEINERTER WENDEPUNKT                     *
C     *       DETA                                                     *
C     *       D2ETA                                                    *
C     *                                                                *
C     *  DIMENSIONSLOSE GROESSEN WURDEN MIT DEN ZUGEHOERIGEN WERTEN    *
C     *  DER AUSSENSTROEMUNG NORMIERT.                                 *
C     *                                                                *
C     *  DAS WAR'S , VIEL ERFOLG !                                     *
C     *                                                                *
C     ******************************************************************
C
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      implicit real*8 (a-h,o-z)
C
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION DELTA1(MYG),DELTA2(MYG),DELTA3(MYG),DELTAU(MYG),
     &          DELTAH(MYG)
      DIMENSION DGDH(MYG),DFDH(MYG)
      DIMENSION PRES1(MYG),PRES2(MYG),PREGS(MYG)
      DIMENSION DRMS1(MYG),DRMS2(MYG),DRMGS(MYG)
      DIMENSION RKY1(MYG),RKY2(MYG),RKY3(MYG),RKY4(MYG),RKY5(MYG),
     &          RKY6(MYG),RKY7(MYG),RKY8(MYG)
      DIMENSION GS(MYG),G(MYG),THETA(MYG),F(MYG)
      DIMENSION RDLS1(MYG),RDLS2(MYG),RDLGS(MYG)
      DIMENSION HS1(MYG),HS2(MYG),HGS(MYG)
      DIMENSION TDLS1(MYG),TDLS2(MYG),TDLGS(MYG)
      DIMENSION US1(MYG),US2(MYG),UGS(MYG)
      DIMENSION T(MYG),RK(MYG),RMUE(MYG),H(MYG)
      DIMENSION DH(MYG),DH2(MYG),DH4(MYG),DH6(MYG),SDH(MYG)
      DIMENSION CPMS1(MYG),CPMS2(MYG),CPMGS(MYG)
      DIMENSION RVS1(MYG),RVS2(MYG),RVGS(MYG),RV(MYG),ETAV(MYG),
     &          ETAV2(MYG)
      DIMENSION DTDES1(MYG),DTDES2(MYG),DTDEGS(MYG)
      DIMENSION GES1(MYG),GES2(MYG),GEGS(MYG)
      DIMENSION DUDES1(MYG),DUDES2(MYG),DUDEGS(MYG)
      DIMENSION TEMP(IMAX),CPDR(IMAX),C(IMAX-1,3),CPST(IMAX),
     &          CPST2(MYG),CPM(MYG)
      DIMENSION VIKOEF(8),RSCH(IRSCH),RMSCH(IRSCH),WENDEP(MYG)
C
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
      COMMON / CDELT / DELTA1,DELTA2,DELTA3,DELTAU,DELTAH
      COMMON / CPR / PRES1,PRES2,PREGS
      COMMON / CVI / DRMS1,DRMS2,DRMGS
      COMMON / CURK / US1,US2,UGS
      COMMON / CTDLRK / TDLS1,TDLS2,TDLGS
      COMMON / CRHORK / RDLS1,RDLS2,RDLGS
      COMMON / CHRK / HS1,HS2,HGS
      COMMON / CDUDH / DUDES1,DUDES2,DUDEGS
      COMMON / CDTDH / DTDES1,DTDES2,DTDEGS
      COMMON / CRV / RVS1,RVS2,RVGS,ETAV,ETAV2,RV
      COMMON / CG / GES1,GES2,GEGS
      COMMON / CABL / DGDH,DFDH
      COMMON / CY / RKY1,RKY2,RKY3,RKY4,RKY5,RKY6,RKY7,RKY8
      COMMON / CINT / GS,G,THETA,F
      COMMON / CDIM / T,RK,RMUE,H
      COMMON / CCP / CPDR,CPST,TEMP,C,CPST2,CPM
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      COMMON / CST / DH,DH2,DH4,DH6,SDH
      character*72 inf(100)
      dimension ltimes(100)
      character*124 filen
C
      DATA RSCH  / 0.12,0.1,0.09,0.08,0.07,0.06,0.05/
      DATA RMSCH / 0.,1.01,2.1,3.01,4.01,5.01,6.01/
C
C *** LESEN DES INPUTS ***
C
      CALL LES(RMSCH,RSCH)
C
C *** KONSTANTE BERECHNEN ***
C
      R=CP-CV
      PI=4.*ATAN(1.)
      GRENZ=0.999
C
C *** ANZAHL DER HERAUSGESCHRIEBENEN FELDER ***
C
      IWINF=19
C
C *** BEI RECHNUNG MIT STRECKUNG WERDEN 2 ZUS. FELDER HERAUSGESCHR. ***
C
      IF(LSTR.EQ.1)  IWINF=21
C
C *** FESTLEGUNG EINER STANDARDTEMPERATUR ***
C
      IF (ICOM.EQ.1) THEN
         T0 = MIN(T1,TW) - 10.
      ELSE IF (ICOM.EQ.0) THEN
         T0 = T1 - 10.
      END IF
C
C *** BERECHNUNG DER SPLINE-KOEFFIZIENTEN ***
C
      IF(ICPPRV.EQ.2) THEN
        CALL BECPST (CPDR,R,IMAX,CPST)
        CALL SPKOEF (TEMP,CPST,IMAX,C)
      END IF
C
C *** FESTLEGUNG DER HERAUSGESCHRIEBENEN FELDER ***********************
C
      LTIMES(1)=1
      INF(1)=' U '
      INF(2)=' T '
      INF(3)=' RHO '
      INF(4)=' FV '
      INF(5)=' DUDH '
      INF(6)=' D2UDH '
      INF(7)=' DTDH '
      INF(8)=' D2TDH '
      INF(9)=' DRHODH '
      INF(10)='D2RHODH '
      INF(11)=' DRMUE '
      INF(12)='DMUEDT '
      INF(13)='D2MUEDT'
      INF(14)='DRD2U '
      INF(15)='PR '
      INF(16)='DRK '
      INF(17)='DRKDT '
      INF(18)='D2RKDT'
      INF(19)='PARAM'
      INF(20)='DETA '
      INF(21)='D2ETA'
C
C *** KENNSATZ SCHREIBEN **********************************************
C
c      CALL WRITEKS (12,1,1,MYG,1,IWINF,LTIMES,INF,100)
c     aufruf IOS :
c          1  -  1.ios-kanal
c          12 -  betriebssystem i/o-unit
c          'profout' - filename (ohne endung)
c          1         - ausgabeart unformattiert,wahlfrei,limittierte blo
c          1,1,1,myg - felddimension
c          iwinf     - anzahl der ausgegebenen parameter
c          ltimes    - zeitfeld (hier nur 1. element belegt)
c          inf       - informationsfeld
c          0         - anzahl zusaetzlicher informationszeilen
      filen = 'profout'
      imach = 2
      !!!!!!call writecd (1,12,filen,imach,1,1,1,myg,iwinf,ltimes,inf,0)
C
      ITEIL=1
C
      IF (LSTR.EQ.0) THEN
C
C *** KONSTANTE SCHRITTWEITE ***
C
        CALL LSTR0(ETA,DN,ITST,ITEIL,ITST1,DN1)
C
      ELSE IF(LSTR.EQ.1) THEN
C
C *** STRECKUNG MIT ASINH ***
C
        IF(ETAS.GT.ETA)                STOP 'ETAS'
        CALL LSTR1(ITST,ITSTC,DETA,D2ETA)
C
      ELSE IF(LSTR.EQ.2) THEN
C
C *** STRECKUNG NACH CHEBYSCHEFF ***
C
        CALL LSTR2(PI,ITST,NSC,ITSTC,NTSTR,NCHEB)
C
      END IF
C
C *** BERECHNUNG KONSTANTER GROESSEN ***
C
      CALL KONST(C5)
C
C *** VORGABE VON ANFANGSBEDINGUNGEN ***
C
      CALL ANKE(CPMS1,CPMS2,CPMGS,ITST)
C
C
C *** BEGINN DER ITERATIONSSCHLEIFE ***
C
      DO 10 K=1,ITMAX
C
C  ** SCHLEIFE UEBER STUETZSTELLEN **
C
        DO 20 I=2,ITST
          M=I-1
C
C *** BERECHNUNG VON ENTHALPIE,TEMP. UND VISKOSITAET AN DER STELLE M
C
          IF(ICPPRV.EQ.0) THEN
            CALL HTM0(H(M),T(M),TDL(M),RMUE(M),DRMUE(M),F(M),DUDH(M),
     &                THETA(M),PR(M))
          ELSE IF(ICPPRV.EQ.1) THEN
            CALL HTM1(H(M),T(M),TDL(M),RMUE(M),DRMUE(M),F(M),DUDH(M),
     &                THETA(M),PR(M))
          ELSE IF(ICPPRV.EQ.2) THEN
            CALL HTM2(H(M),T(M),TDL(M),RMUE(M),DRMUE(M),F(M),DUDH(M),
     &                THETA(M),PR(M),CPM(M),M)
          END IF
          CALL FUGS(DGSDHM,U(M),TDL(M))
          CALL FUF(DFDH(M),GS(M),DRMUE(M),F(M))
          CALL FUTET(DTEDHM,PR(M),DRMUE(M),G(M))
          CALL FUG(DGDH(M),PR(M),DRMUE(M),GS(M),G(M),DUDH(M))
          CALL FUY1(DY1DHM,DRMUE(M),RKY5(M))
          CALL FUY2(DY2DHM,DRMUE(M),RKY6(M))
          CALL FUY3(DY3DHM,PR(M),DRMUE(M),RKY7(M))
          CALL FUY4(DY4DHM,PR(M),DRMUE(M),RKY8(M))
          CALL FUY5(DY5DHM,GS(M),DRMUE(M),RKY5(M))
          CALL FUY6(DY6DHM,GS(M),DRMUE(M),RKY6(M))
          CALL FUY7(DY7DHM,DRMUE(M),F(M),PR(M),GS(M),RKY5(M),
     &              RKY7(M))
          CALL FUY8(DY8DHM,DRMUE(M),F(M),PR(M),GS(M),RKY6(M),
     &              RKY8(M))
C
C *** RUNGE KUTTA VERFAHREN 1.HALBSCHRITT *****************************
C
          CALL RK12G(GSS1,GS(M),DH2(M),DGSDHM)
          CALL RK12G(US1(I),U(M),DH2(M),DUDH(M))
          CALL RK12G(FS1,F(M),DH2(M),DFDH(M))
          CALL RK12G(TETS1,THETA(M),DH2(M),DTEDHM)
          CALL RK12G(GS1,G(M),DH2(M),DGDH(M))
          GES1(I)=GS1
          CALL RK12G(Y1S1,RKY1(M),DH2(M),DY1DHM)
          CALL RK12G(Y2S1,RKY2(M),DH2(M),DY2DHM)
          CALL RK12G(Y3S1,RKY3(M),DH2(M),DY3DHM)
          CALL RK12G(Y4S1,RKY4(M),DH2(M),DY4DHM)
          CALL RK12G(Y5S1,RKY5(M),DH2(M),DY5DHM)
          CALL RK12G(Y6S1,RKY6(M),DH2(M),DY6DHM)
          CALL RK12G(Y7S1,RKY7(M),DH2(M),DY7DHM)
          CALL RK12G(Y8S1,RKY8(M),DH2(M),DY8DHM)
C
C *** BERECHNUNG VON ENTHALPIE,VISKOSITAET UND TEMPERATUR FUER
C     2.HALBSCHRITT **************************************************
C
          IF(ICPPRV.EQ.0) THEN
            CALL HTM0(HS1(I),TS1,TDLS1(I),RMUES1,DRMS1,FS1,DUDHS1,
     &                TETS1,PRS1)
          ELSE IF(ICPPRV.EQ.1) THEN
            CALL HTM1(HS1(I),TS1,TDLS1(I),RMUES1,DRMS1,FS1,DUDHS1,
     &                TETS1,PRS1)
          ELSE IF(ICPPRV.EQ.2) THEN
            CALL HTM2(HS1(I),TS1,TDLS1(I),RMUES1,DRMS1,FS1,DUDHS1,
     &                TETS1,PRS1,CPMS1(I),M)
          END IF
          DUDES1(I)=DUDHS1
          PRES1(I)=PRS1
          CALL FUGS(DGSDH1,US1(I),TDLS1(I))
          CALL FUF(DFDH1,GSS1,DRMS1,FS1)
          CALL FUTET(DTEDH1,PRS1,DRMS1,GS1)
          CALL FUG(DGDH1,PRS1,DRMS1,GSS1,GS1,DUDHS1)
          CALL FUY1(DY1DH1,DRMS1,Y5S1)
          CALL FUY2(DY2DH1,DRMS1,Y6S1)
          CALL FUY3(DY3DH1,PRS1,DRMS1,Y7S1)
          CALL FUY4(DY4DH1,PRS1,DRMS1,Y8S1)
          CALL FUY5(DY5DH1,GSS1,DRMS1,Y5S1)
          CALL FUY6(DY6DH1,GSS1,DRMS1,Y6S1)
          CALL FUY7(DY7DH1,DRMS1,FS1,PRS1,GSS1,Y5S1,Y7S1)
          CALL FUY8(DY8DH1,DRMS1,FS1,PRS1,GSS1,Y6S1,Y8S1)
C
C *** RUNGE KUTTA VERFAHREN 2.HALBSCHRITT *****************************
C
          CALL RK12G(GSS2,GS(M),DH2(M),DGSDH1)
          CALL RK12G(US2(I),U(M),DH2(M),DUDHS1)
          CALL RK12G(FS2,F(M),DH2(M),DFDH1)
          CALL RK12G(TETS2,THETA(M),DH2(M),DTEDH1)
          CALL RK12G(GS2,G(M),DH2(M),DGDH1)
          GES2(I)=GS2
          CALL RK12G(Y1S2,RKY1(M),DH2(M),DY1DH1)
          CALL RK12G(Y2S2,RKY2(M),DH2(M),DY2DH1)
          CALL RK12G(Y3S2,RKY3(M),DH2(M),DY3DH1)
          CALL RK12G(Y4S2,RKY4(M),DH2(M),DY4DH1)
          CALL RK12G(Y5S2,RKY5(M),DH2(M),DY5DH1)
          CALL RK12G(Y6S2,RKY6(M),DH2(M),DY6DH1)
          CALL RK12G(Y7S2,RKY7(M),DH2(M),DY7DH1)
          CALL RK12G(Y8S2,RKY8(M),DH2(M),DY8DH1)
C
C *** BERECHNUNG VON ENTHALPIE,VISKOSITAET UND TEMPERATUR FUER
C     1 VOLLSCHRITT **************************************************
C
          IF(ICPPRV.EQ.0) THEN
            CALL HTM0(HS2(I),TS2,TDLS2(I),RMUES2,DRMS2,FS2,DUDHS2,
     &                TETS2,PRS2)
          ELSE IF(ICPPRV.EQ.1) THEN
            CALL HTM1(HS2(I),TS2,TDLS2(I),RMUES2,DRMS2,FS2,DUDHS2,
     &                TETS2,PRS2)
          ELSE IF(ICPPRV.EQ.2) THEN
            CALL HTM2(HS2(I),TS2,TDLS2(I),RMUES2,DRMS2,FS2,DUDHS2,
     &                TETS2,PRS2,CPMS2(I),M)
          END IF
          DUDES2(I)=DUDHS2
          PRES2(I)=PRS2
          CALL FUGS(DGSDH2,US2(I),TDLS2(I))
          CALL FUF(DFDH2,GSS2,DRMS2,FS2)
          CALL FUTET(DTEDH2,PRS2,DRMS2,GS2)
          CALL FUG(DGDH2,PRS2,DRMS2,GSS2,GS2,DUDHS2)
          CALL FUY1(DY1DH2,DRMS2,Y5S2)
          CALL FUY2(DY2DH2,DRMS2,Y6S2)
          CALL FUY3(DY3DH2,PRS2,DRMS2,Y7S2)
          CALL FUY4(DY4DH2,PRS2,DRMS2,Y8S2)
          CALL FUY5(DY5DH2,GSS2,DRMS2,Y5S2)
          CALL FUY6(DY6DH2,GSS2,DRMS2,Y6S2)
          CALL FUY7(DY7DH2,DRMS2,FS2,PRS2,GSS2,Y5S2,Y7S2)
          CALL FUY8(DY8DH2,DRMS2,FS2,PRS2,GSS2,Y6S2,Y8S2)
C
C *** RUNGE KUTTA VERFAHREN 1.VOLLSCHRITT *****************************
C
          CALL RK12G(GSGS,GS(M),DH(M),DGSDH2)
          CALL RK12G(UGS(I),U(M),DH(M),DUDHS2)
          CALL RK12G(FGS,F(M),DH(M),DFDH2)
          CALL RK12G(TETGS,THETA(M),DH(M),DTEDH2)
          CALL RK12G(GGS,G(M),DH(M),DGDH2)
          GEGS(I)=GGS
          CALL RK12G(Y1GS,RKY1(M),DH(M),DY1DH2)
          CALL RK12G(Y2GS,RKY2(M),DH(M),DY2DH2)
          CALL RK12G(Y3GS,RKY3(M),DH(M),DY3DH2)
          CALL RK12G(Y4GS,RKY4(M),DH(M),DY4DH2)
          CALL RK12G(Y5GS,RKY5(M),DH(M),DY5DH2)
          CALL RK12G(Y6GS,RKY6(M),DH(M),DY6DH2)
          CALL RK12G(Y7GS,RKY7(M),DH(M),DY7DH2)
          CALL RK12G(Y8GS,RKY8(M),DH(M),DY8DH2)
C
C *** BERECHNUNG VON ENTHALPIE,VISKOSITAET UND TEMPERATUR FUER
C     1.HALBSCHRITT **************************************************
C
          IF(ICPPRV.EQ.0) THEN
            CALL HTM0(HGS(I),TGS,TDLGS(I),RMUEGS,DRMGS,FGS,DUDHGS,
     &                TETGS,PRGS)
          ELSE IF(ICPPRV.EQ.1) THEN
            CALL HTM1(HGS(I),TGS,TDLGS(I),RMUEGS,DRMGS,FGS,DUDHGS,
     &                TETGS,PRGS)
          ELSE IF(ICPPRV.EQ.2) THEN
            CALL HTM2(HGS(I),TGS,TDLGS(I),RMUEGS,DRMGS,FGS,DUDHGS,
     &                TETGS,PRGS,CPMGS(I),M)
          END IF
          DUDEGS(I)=DUDHGS
          PREGS(I)=PRGS
          CALL FUGS(DGSDHG,UGS(I),TDLGS(I))
          CALL FUF(DFDHG,GSGS,DRMGS,FGS)
          CALL FUTET(DTEDHG,PRGS,DRMGS,GGS)
          CALL FUG(DGDHG,PRGS,DRMGS,GSGS,GGS,DUDHGS)
          CALL FUY1(DY1DHG,DRMGS,Y5GS)
          CALL FUY2(DY2DHG,DRMGS,Y6GS)
          CALL FUY3(DY3DHG,PRGS,DRMGS,Y7GS)
          CALL FUY4(DY4DHG,PRGS,DRMGS,Y8GS)
          CALL FUY5(DY5DHG,GSGS,DRMGS,Y5GS)
          CALL FUY6(DY6DHG,GSGS,DRMGS,Y6GS)
          CALL FUY7(DY7DHG,DRMGS,FGS,PRGS,GSGS,Y5GS,Y7GS)
          CALL FUY8(DY8DHG,DRMGS,FGS,PRGS,GSGS,Y6GS,Y8GS)
C
C *** RUNGE-KUTTA SCHRITT Y(N+1) ***
C
          CALL RKYNP1(GS(I),GS(M),DH6(M),DGSDHM,DGSDH1,DGSDH2,DGSDHG)
          CALL RKYNP1(U(I),U(M),DH6(M),DUDH(M),DUDHS1,DUDHS2,DUDHGS)
          CALL RKYNP1(F(I),F(M),DH6(M),DFDH(M),DFDH1,DFDH2,DFDHG)
          CALL RKYNP1(THETA(I),THETA(M),DH6(M),DTEDHM,DTEDH1,
     &                DTEDH2,DTEDHG)
          CALL RKYNP1(G(I),G(M),DH6(M),DGDH(M),DGDH1,DGDH2,DGDHG)
          CALL RKYNP1(RKY1(I),RKY1(M),DH6(M),DY1DHM,DY1DH1,DY1DH2,
     &                DY1DHG)
          CALL RKYNP1(RKY2(I),RKY2(M),DH6(M),DY2DHM,DY2DH1,DY2DH2,
     &                DY2DHG)
          CALL RKYNP1(RKY3(I),RKY3(M),DH6(M),DY3DHM,DY3DH1,DY3DH2,
     &                DY3DHG)
          CALL RKYNP1(RKY4(I),RKY4(M),DH6(M),DY4DHM,DY4DH1,DY4DH2,
     &                DY4DHG)
          CALL RKYNP1(RKY5(I),RKY5(M),DH6(M),DY5DHM,DY5DH1,DY5DH2,
     &                DY5DHG)
          CALL RKYNP1(RKY6(I),RKY6(M),DH6(M),DY6DHM,DY6DH1,DY6DH2,
     &                DY6DHG)
          CALL RKYNP1(RKY7(I),RKY7(M),DH6(M),DY7DHM,DY7DH1,DY7DH2,
     &                DY7DHG)
          CALL RKYNP1(RKY8(I),RKY8(M),DH6(M),DY8DHM,DY8DH1,DY8DH2,
     &                DY8DHG)
C
C *** BERECHNUNG VON ENTHALPIE,VISKOSITAET UND TEMPERATUR FUER
C     SCHRITT Y(N+1) *******************************************
C
          IF(I.EQ.ITST) THEN
            IF(ICPPRV.EQ.0) THEN
              CALL HTM0(H(I),T(I),TDL(I),RMUE(I),DRMUE(I),F(I),DUDH(I),
     &                  THETA(I),PR(I))
            ELSE IF(ICPPRV.EQ.1) THEN
              CALL HTM1(H(I),T(I),TDL(I),RMUE(I),DRMUE(I),F(I),DUDH(I),
     &                THETA(I),PR(I))
            ELSE IF(ICPPRV.EQ.2) THEN
              CALL HTM2(H(I),T(I),TDL(I),RMUE(I),DRMUE(I),F(I),DUDH(I),
     &                  THETA(I),PR(I),CPM(I),M)
            END IF
            CALL FUF(DFDH(I),GS(I),DRMUE(I),F(I))
            CALL FUG(DGDH(I),PR(I),DRMUE(I),GS(I),G(I),DUDH(I))
          END IF
C
   20   CONTINUE
C
C *** SCHLEIFE FUER STUETZSTELLEN BEENDET *****************************
C
C *** NEUE ANFANGSBEDINGUNGEN BERECHNEN ***********************
C
        CALL NEUAB(CPMS1,CPMS2,CPMGS,ITEND,ITST,DA,DB,TDLS1,TDLS2,
     &             TDLGS)
C
        IF(ITEND.EQ.1) GOTO 100
   10 CONTINUE
C
C *** ITERATIONSSCHLEIFE BEENDET ***************************************
C
      PRINT *
      PRINT *,'     KEINE KONVERGENZ ERZIELT !'
      STOP 'KONV'
  100 CONTINUE
C
      KIT=K-1
C
C *** LETZTMALIGE INTEGRATION ZUR BERECHNUNG VON FV UND DICHTE *********
C
      CALL FVDI(DH6,ITST)
C
C *** BERECHNUNG DER GRENZSCHICHTDICKE DELTA RUNGE KUTTA O(4) **********
C
C *** BERECHNUNG DER GRENZSCHICHTKONSTANTEN ****************************
C
      CALL DELTI(DELTA,DELT1G,DELT2G,DELT3G,DELTUG,DELTHG,ITST,GRENZ)
C
C *** BERECHNUNG DER ABLEITUNGEN VON U,T UND RHO **********************
C
      CALL ABLEIT(ITST)
C
C *** BERECHNUNG DES VERALLGEMEINERTEN WENDEPUNKTS ********************
C
      CALL WENDE(ITST,ICOM,WENDEP,JW)
C
C *** RUECKTRANSFORMATION DER AUSZUDRUCKENDEN FELDER BEI SCHRITTWEITE
C     GROESSER 0.1 ****************************************************
C
      IF((LSTR.EQ.0).AND.(ITEIL.GT.1)) THEN
        CALL FEAEND(ITST1,ITST,ITEIL,DN,DN1)
      END IF
C
C *** AUSDRUCKEN DER ERRECHNETEN WERTE ********************************
C
      CALL DRUCK(DELTA,DELT1G,DELT2G,DELT3G,DELTUG,DELTHG,WENDEP,ITST,
     &           KIT,DA,DB,JW)
C
      !!!!!!CALL writed (1,U)
      !!!!!!CALL writed (1,TDL)
      !!!!!!CALL writed (1,RHODL)
      !!!!!!CALL writed (1,RV)
      !!!!!!CALL writed (1,DUDH)
      !!!!!!CALL writed (1,D2UDH)
      !!!!!!CALL writed (1,DTDH)
      !!!!!!CALL writed (1,D2TDH)
      !!!!!!CALL writed (1,DRHODH)
      !!!!!!CALL writed (1,D2RDH)
      !!!!!!CALL writed (1,DRMUE)
      !!!!!!CALL writed (1,DMUEDT)
      !!!!!!CALL writed (1,D2MDT)
      !!!!!!CALL writed (1,DRD2U)
      !!!!!!CALL writed (1,PR)
      !!!!!!CALL writed (1,DRK)
      !!!!!!CALL writed (1,DRKDT)
      !!!!!!CALL writed (1,D2RKDT)
      !!!!!!CALL writed (1,pr)

      IF (LSTR.EQ.1) THEN
        !!!!!!CALL writed (1,DETA)
        !!!!!!CALL writed (1,D2ETA)
      END IF

      open(100,file='profout.dat',status='unknown')
      do i=1,1000
        write(100,'(f30.16,f30.16,f30.16,f30.16)') 
     &    U(i),TDL(i),RHODL(i),RV(i)
      end do
      close(100)
 
      STOP 'NORMAL'
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG DER ABLEITUNGEN VON U,T,RHO UND MUE
C *********************************************************************
C
      SUBROUTINE ABLEIT(ITST)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION DGDH(MYG),DFDH(MYG)
      DIMENSION RVS1(MYG),RVS2(MYG),RVGS(MYG),RV(MYG),ETAV(MYG),
     &          ETAV2(MYG)
      DIMENSION DTDHS1(MYG),DTDHS2(MYG),DTDHGS(MYG)
      DIMENSION DH(MYG),DH2(MYG),DH4(MYG),DH6(MYG),SDH(MYG)
      DIMENSION TDLS1(MYG),TDLS2(MYG),TDLGS(MYG)
      DIMENSION DUDHS1(MYG),DUDHS2(MYG),DUDHGS(MYG)
      DIMENSION TEMP(IMAX),CPDR(IMAX),C(IMAX-1,3),CPST(IMAX),
     &          CPST2(MYG),CPM(MYG)
      DIMENSION T(MYG),RK(MYG),RMUE(MYG),H(MYG)
      DIMENSION GES1(MYG),GES2(MYG),GEGS(MYG)
      DIMENSION GS(MYG),G(MYG),THETA(MYG),F(MYG)
      DIMENSION DRMS1(MYG),DRMS2(MYG),DRMGS(MYG)
      DIMENSION UHS1(MYG),UHS2(MYG),UGS(MYG)
      DIMENSION PRS1(MYG),PRS2(MYG),PRGS(MYG)
      DIMENSION VIKOEF(8)
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
      COMMON / CABL / DGDH,DFDH
      COMMON / CDUDH / DUDHS1,DUDHS2,DUDHGS
      COMMON / CDTDH / DTDHS1,DTDHS2,DTDHGS
      COMMON / CRV / RVS1,RVS2,RVGS,ETAV,ETAV2,RV
      COMMON / CURK / UHS1,UHS2,UGS
      COMMON / CST / DH,DH2,DH4,DH6,SDH
      COMMON / CG / GES1,GES2,GEGS
      COMMON / CPR / PRS1,PRS2,PRGS
      COMMON / CINT / GS,G,THETA,F
      COMMON / CTDLRK / TDLS1,TDLS2,TDLGS
      COMMON / CDIM / T,RK,RMUE,H
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CCP / CPDR,CPST,TEMP,C,CPST2,CPM
      COMMON / CVI/ DRMS1,DRMS2,DRMGS
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
C
C *** KONSTANTEN ******************************************************
C
      C6=0.75*2.454E2
      C7=T1*C1/2./RMUE1
      C8=T1*C7
      C9=DLOG(1.d1)
      C10=T1/RK1*C2
      C11=C10*T1
C
      TN1=1.104D2
      TN2=110.4*TN1
      TN3=110.4*TN2
      TN4=110.4*TN3
      TN5=110.4*TN4
      TN6=110.4*TN5
      TN7=110.4*TN6
C
      C25=T1/RMUE1*EMUE
      C26=T1*C25
C
C *** BERECHNUNG DER ABLEITUNGEN VON MUE UND KAPPA NACH DER TEMPERATUR
C
      DO 51 I=1,ITST
        IF(ICPPRV.GE.1) THEN
          CALL RKFT(T(I),C2,C4,RK(I))
        ELSE IF(ICPPRV.EQ.0) THEN
          CALL RKFMUE(RMUE(I),RK(I))
        END IF
        DRK(I)=RK(I)/RK1
        IF((T(I).GT.80.).AND.(ICPPRV.GE.1)) THEN
          C12=10.**(-12./T(I))
          C13=245.4*C12
          C14=T(I)+C13
          C15=12.*C9*C13
          DRKDTN=1.5*T(I)*SQRT(T(I))+1.5*SQRT(T(I))*C13-T(I)*
     &           SQRT(T(I))-C15/SQRT(T(I))
          DRKDT(I)=C10*DRKDTN/C14/C14
          D2RKDT(I)=C11/C14/C14/C14*((1.+C15/T(I)/T(I))*DRKDTN-C14*
     &              (2.25*SQRT(T(I))+1.5*245.4*(0.5/SQRT(T(I))*C12+
     &              SQRT(T(I))*C12*C9*12./T(I)/T(I))-1.5*SQRT(T(I))-
     &              245.4*C9*12.*(C12*C9*12./SQRT(T(I))/T(I)/T(I)-0.5*
     &              C12/SQRT(T(I))/T(I))))
        ELSE IF((T(I).LE.80.).AND.(ICPPRV.GE.1)) THEN
          DRKDT(I)=C4*T1/RK1
          D2RKDT(I)=0.
        END IF
        IF(T(I).GT.TV2) THEN
          C19=T(I)+110.4
          C20=T(I)+331.2
          DMUEDT(I)=C7*SQRT(T(I))*C20/C19/C19
          D2MDT(I)=C8*((1.5*SQRT(T(I))+165.6/SQRT(T(I)))*C19-2.*
     &             SQRT(T(I))*C20)/C19/C19/C19
        ELSE IF(T(I).GT.TV1) THEN
          DMUEDT(I)=C25*(7.*VIKOEF(1)*T(I)**6/TN7+6.*VIKOEF(2)*T(I)**5/
     &              TN6+5.*VIKOEF(3)*T(I)**4/TN5+4.*VIKOEF(4)*T(I)**3/
     &              TN4+3.*VIKOEF(5)*T(I)**2/TN3+2.*VIKOEF(6)*T(I)/TN2+
     &              VIKOEF(7)/TN1)
          D2MDT(I)=C26*(42.*VIKOEF(1)*T(I)**5/TN7+30.*VIKOEF(2)*
     &             T(I)**4/TN6+20.*VIKOEF(3)*T(I)**3/TN5+12.*VIKOEF(4)
     &             *T(I)**2/TN4+6.*VIKOEF(5)*T(I)/TN3+2.*VIKOEF(6)/
     &             TN2)
        ELSE
          DMUEDT(I)=C3*T1/RMUE1
          D2MDT(I)=0.
        END IF
        IF(ICPPRV.EQ.0) THEN
          DRKDT(I)=DMUEDT(I)
          D2RKDT(I)=D2MDT(I)
        END IF
   51 CONTINUE
      C21=DH01H1/T1
      IF(ICPPRV.EQ.0) THEN
        DO 53 I=1,ITST
          DTDH(I)=C21*PR(I)*G(I)/CP/DRMUE(I)
          IF(I.GE.2) THEN
            CT=TDLS1(I)*T1
            CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,CT,CM)
            DRMS1(I)=CM/RMUE1
            CT=TDLS2(I)*T1
            CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,CT,CM)
            DRMS2(I)=CM/RMUE1
            CT=TDLGS(I)*T1
            CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,CT,CM)
            DRMGS(I)=CM/RMUE1
            DTDHS1(I)=C21*PRS1(I)*GES1(I)/CP/DRMS1(I)
            DTDHS2(I)=C21*PRS2(I)*GES2(I)/CP/DRMS2(I)
            DTDHGS(I)=C21*PRGS(I)*GEGS(I)/CP/DRMGS(I)
          END IF
          D2TDH(I)=C21*(DGDH(I)*DRMUE(I)-G(I)*DMUEDT(I)*DTDH(I))/
     &             CP/DRMUE(I)/DRMUE(I)*PR(I)
          D2UDH(I)=(DFDH(I)*DRMUE(I)-F(I)*DMUEDT(I)*DTDH(I))/DRMUE(I)/
     &             DRMUE(I)
   53   CONTINUE
      ELSE
        DO 55 I=1,ITST
          IF(ICPPRV.EQ.2) THEN
            CALL CPINT(T(I),C,IMAX,TEMP,CPST,CPST2(I))
            CT=TDLS1(I)*T1
            CALL CPINT(CT,C,IMAX,TEMP,CPST,CPS1)
            CT=TDLS2(I)*T1
            CALL CPINT(CT,C,IMAX,TEMP,CPST,CPS2)
            CT=TDLGS(I)*T1
            CALL CPINT(CT,C,IMAX,TEMP,CPST,CPGS)
          ELSE
            CPST2(I)=CP
            CPS1=CP
            CPS2=CP
            CPGS=CP
          END IF
          DTDH(I)=C21*PR(I)*G(I)/DRMUE(I)/CPST2(I)
          IF(I.GE.2) THEN
            CT=TDLS1(I)*T1
            CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,CT,CM)
            DRMS1(I)=CM/RMUE1
            CT=TDLS2(I)*T1
            CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,CT,CM)
            DRMS2(I)=CM/RMUE1
            CT=TDLGS(I)*T1
            CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,CT,CM)
            DRMGS(I)=CM/RMUE1
            DTDHS1(I)=C21*PRS1(I)*GES1(I)/CPS1/DRMS1(I)
            DTDHS2(I)=C21*PRS2(I)*GES2(I)/CPS2/DRMS2(I)
            DTDHGS(I)=C21*PRGS(I)*GEGS(I)/CPGS/DRMGS(I)
          END IF
          D2TDH(I)=RMUE1*C21/RK(I)/RK(I)*(DGDH(I)*RK(I)-G(I)*DRKDT(I)
     &             *DTDH(I)*RK1)
          D2UDH(I)=(DFDH(I)*DRMUE(I)-F(I)*DMUEDT(I)*DTDH(I))/DRMUE(I)/
     &             DRMUE(I)
   55   CONTINUE
      END IF
C
C *** BERECHNUNG VON RV ***
C
      ETAV(1)=0.
      DO 54 I=2,ITST
        M=I-1
        ETAV(I)=DH(M)+ETAV(M)
        ETAV2(M)=ETAV(M)+DH2(M)
   54 CONTINUE
      RV(1)=0.
      DO 50 I=2,ITST
        M=I-1
        RVS1(I)=RV(M)+DH2(M)*( DUDH  (M)/TDL  (M)-    U(M) *
     &                         DTDH  (M)/TDL  (M)/TDL  (M))*ETAV (M)
        RVS2(I)=RV(M)+DH2(M)*( DUDHS1(I)/TDLS1(I)- UHS1(I) *
     &                         DTDHS1(I)/TDLS1(I)/TDLS1(I))*ETAV2(M)
        RVGS(I)=RV(M)+DH (M)*( DUDHS2(I)/TDLS2(I)- UHS2(I) *
     &                         DTDHS2(I)/TDLS2(I)/TDLS2(I))*ETAV2(M)
        RV(I)=RV(M)+  DH6(M)*(
     &                       ( DUDH  (M)/TDL  (M)-    U(M) *
     &                         DTDH  (M)/TDL  (M)/TDL  (M))*ETAV (M) +
     &                    2.*( DUDHS1(I)/TDLS1(I)- UHS1(I) *
     &                         DTDHS1(I)/TDLS1(I)/TDLS1(I))*ETAV2(M) +
     &                    2.*( DUDHS2(I)/TDLS2(I)- UHS2(I) *
     &                         DTDHS2(I)/TDLS2(I)/TDLS2(I))*ETAV2(M) +
     &                       ( DUDHGS(I)/TDLGS(I)-  UGS(I) *
     &                         DTDHGS(I)/TDLGS(I)/TDLGS(I))*ETAV (I)  )
   50 CONTINUE
      DO 57 I=1,ITST
        RV(I)=RV(I)/RHODL(I)/2.
   57 CONTINUE
C
C *** ABLEITUNGEN VON RHO NACH ETA ***
C
      DO 58 I=1,ITST
        DRHODH(I)=-DTDH(I)/TDL(I)/TDL(I)
        D2RDH(I)=(2.*DTDH(I)*DTDH(I)-TDL(I)*D2TDH(I))/TDL(I)/TDL(I)/
     &             TDL(I)
   58 CONTINUE
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR VORGABE VON ANFANGS- UND RANDBEDINGUNGEN
C *********************************************************************
C
      SUBROUTINE ANKE(CPMS1,CPMS2,CPMGS,ITST)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION TEMP(IMAX),CPDR(IMAX),C(IMAX-1,3),CPST(IMAX),
     &          CPST2(MYG)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION RKY1(MYG),RKY2(MYG),RKY3(MYG),RKY4(MYG),RKY5(MYG),
     &          RKY6(MYG),RKY7(MYG),RKY8(MYG)
      DIMENSION GS(MYG),G(MYG),THETA(MYG),F(MYG)
      DIMENSION CPMS1(MYG),CPMS2(MYG),CPMGS(MYG),CPM(MYG)
      DIMENSION T(MYG),RK(MYG),RMUE(MYG),H(MYG)
      DIMENSION VIKOEF(8)
      COMMON / CDIM / T,RK,RMUE,H
      COMMON / CY / RKY1,RKY2,RKY3,RKY4,RKY5,RKY6,RKY7,RKY8
      COMMON / CINT / GS,G,THETA,F
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      COMMON / CCP / CPDR,CPST,TEMP,C,CPST2,CPM
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
C
      U(1)=0.
      GS(1)=0.
      RKY1(1)=0.
      RKY2(1)=0.
      RKY3(1)=0.
      RKY5(1)=1.
      RKY6(1)=0.
      RKY7(1)=0.
      DUDH(1)=DUDHS
C
      IF (ICOM.EQ.1) THEN
        G(1)=GSCH
        RKY4(1)=0.
        RKY8(1)=1.
C
        IF((ICPPRV.EQ.1).OR.(ICPPRV.EQ.0)) THEN
          THETA(1)=CP*(TW-T1)/DH01H1
          H(1)=H1+DH01H1*THETA(1)
          T(1)=H(1)/CP+T0
        ELSE
          CALL KUBINT(TW,T1,C,IMAX,IINT,EPS,TEMP,CPST,CPMI)
          THETA(1)=CPMI*(TW-T1)/DH01H1
          H(1)=H1+DH01H1*THETA(1)
          T(1)=TW
        END IF
C
      ELSE
        G(1)=0.
        RKY4(1)=1.
        RKY8(1)=0.
        THETA(1)=THETAS
        H(1)=H1+DH01H1*THETA(1)
        T(1)=H(1)/CP+T0
      END IF
C
C
C *** T DIMENSIONSLOS MIT TEMPERATUR DER AUSSENSTROEMUNG **************
C
      TDL(1)=T(1)/T1
C
      CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,T(1),RMUE(1))
C
      DRMUE(1)=RMUE(1)/RMUE1
      F(1)=DUDH(1)*DRMUE(1)
C
      IF(ICPPRV.EQ.2) THEN
C
        IF(ICOM.EQ.0) THEN

C *** VORGABE VON CP=KONSTANT ALS STARTWERT FUER BERECHNUNG MIT
C     VARIABLEM CP BEI ADIABATER WAND *******************************
C
          DO 10 I=1,ITST
            CPM(I)=CP
            CPMS1(I)=CP
            CPMS2(I)=CP
            CPMGS(I)=CP
   10     CONTINUE
        ELSE
C
C *** VORGABE DES WERTES VON CP AN DER WAND BEI RECHNUNG MIT KONSTANTER
C     WANDTEMPERATUR ***
C
          DO 20 I=1,ITST
            CPM(I)=(THETA(1)*DH01H1+H1)/(T(1)-T0)
            CPMS1(I)=CPM(I)
            CPMS2(I)=CPM(I)
            CPMGS(I)=CPM(I)
   20     CONTINUE

        END IF
C
      END IF
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG DER DISKRETEN WERTE VON CP(T)
C *********************************************************************
C
      SUBROUTINE BECPST (CPDR,R,IMAX,CPST)
      implicit real*8 (a-h,o-z)
      DIMENSION CPDR(IMAX),CPST(IMAX)
C
      DO 1000 I=1,IMAX
        CPST(I)=CPDR(I)*R
 1000 CONTINUE
C
      END
C
C *********************************************************************
C     INTERPOLATION VON CP
C *********************************************************************
C
      SUBROUTINE CPINT(T,C,IMAX,TEMP,CPST,CPI)
      implicit real*8 (a-h,o-z)
      DIMENSION C(IMAX-1,3),TEMP(IMAX),CPST(IMAX)
C
C *** SUCHEN DER TEMPERATURSTUETZSTELLE *******************************
C
      DO 100 I=1,IMAX-1
         IF ((TEMP(I).LE.T).AND.(T.LT.TEMP(I+1))) GOTO 200
  100 CONTINUE
  200 CONTINUE
C
      D=T-TEMP(I)
      CPI=((C(I,3)*D+C(I,2))*D+C(I,1))*D+CPST(I)
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG DER GRENZSCHICHTDICKE UND DER
C     GRENZSCHICHTKONSTANTEN
C *********************************************************************
C
      SUBROUTINE DELTI(DELTA,DELT1G,DELT2G,DELT3G,DELTUG,DELTHG,ITST,
     &                 GRENZ)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION T(MYG),RK(MYG),RMUE(MYG),H(MYG)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION DELTA1(MYG),DELTA2(MYG),DELTA3(MYG),DELTAU(MYG),
     &          DELTAH(MYG)
      DIMENSION US1(MYG),US2(MYG),UGS(MYG)
      DIMENSION HS1(MYG),HS2(MYG),HGS(MYG)
      DIMENSION RDLS1(MYG),RDLS2(MYG),RDLGS(MYG)
      DIMENSION DH(MYG),DH2(MYG),DH4(MYG),DH6(MYG),SDH(MYG)
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
      COMMON / CDELT / DELTA1,DELTA2,DELTA3,DELTAU,DELTAH
      COMMON / CURK / US1,US2,UGS
      COMMON / CDIM / T,RK,RMUE,H
      COMMON / CHRK / HS1,HS2,HGS
      COMMON / CRHORK / RDLS1,RDLS2,RDLGS
      COMMON / CST / DH,DH2,DH4,DH6,SDH
C
C *** VORGABE VON ANFANGSBEDINGUNGEN **********************************
C
      DELTA1(1)=0.
      DELTA2(1)=0.
      DELTA3(1)=0.
      DELTAH(1)=0.
      DELTAU(1)=0.
C
      C1=RHODL(ITST)*U(ITST)
C
C *** INTEGRATION DER DIFFERENTIALGLEICHUNGEN DURCH RUNGE KUTTA O(4) **
C
      DO 50 I=2,ITST
        M=I-1
        C2HS1=US1(I)/U(ITST)
        C2HS2=US2(I)/U(ITST)
        C2GS=UGS(I)/U(ITST)
        C2=U(I)/U(ITST)
        PRODS1=RDLS1(I)*US1(I)/C1
        PRODS2=RDLS2(I)*US2(I)/C1
        PRODGS=RDLGS(I)*UGS(I)/C1
        PROD=RHODL(I)*U(I)/C1
        DELTA1(I)=DELTA1(M)+DH6(M)*(1.-PRODS1+2.*(2.-PRODS2-PRODGS)
     &            +1.-PROD)
        DELTA2(I)=DELTA2(M)+DH6(M)*(PRODS1*(1.-US1(I)/U(ITST))+2.*(
     &            PRODS2*(1.-US2(I)/U(ITST))+PRODGS*(1.-UGS(I)/
     &            U(ITST)))+PROD*(1.-U(I)/U(ITST)))
        DELTA3(I)=DELTA3(M)+DH6(M)*(PRODS1*(1.-C2HS1*C2HS1)+2.*
     &            (PRODS2*(1.-C2HS2*C2HS2)+PRODGS*(1.-C2GS*C2GS))+
     &            PROD*(1.-C2*C2))
        DELTAH(I)=DELTAH(M)+DH6(M)*(PRODS1*(HS1(I)/H(ITST)-1.)+2.*
     &            (PRODS2*(HS2(I)/H(ITST)-1.)+PRODGS*(HGS(I)/H(ITST)
     &            -1.))+PROD*(H(I)/H(ITST)-1.))
        DELTAU(I)=DELTAU(M)+DH6(M)*(1.-C2HS1+2.*(2.-C2HS2-C2GS)+1.-
     &            C2)
   50 CONTINUE
      ETA1=0.
      DO 100 I=1,ITST
        IF(U(I).GT.GRENZ) GOTO 500
        ETA1=ETA1+DH(I)
  100 CONTINUE
  500 CONTINUE
      ETA1=ETA1-DH(I-1)
      DELTA=ETA1+(GRENZ-U(I-1))/(U(I)-U(I-1))*DH(I-1)
      DELT1G=DELTA1(I-1)+(DELTA-ETA1)/DH(I-1)*(DELTA1(I)-DELTA1(I-1))
      DELT2G=DELTA2(I-1)+(DELTA-ETA1)/DH(I-1)*(DELTA2(I)-DELTA2(I-1))
      DELT3G=DELTA3(I-1)+(DELTA-ETA1)/DH(I-1)*(DELTA3(I)-DELTA3(I-1))
      DELTHG=DELTAH(I-1)+(DELTA-ETA1)/DH(I-1)*(DELTAH(I)-DELTAH(I-1))
      DELTUG=DELTAU(I-1)+(DELTA-ETA1)/DH(I-1)*(DELTAU(I)-DELTAU(I-1))
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUM AUSDRUCKEN
C *********************************************************************
C
      SUBROUTINE DRUCK(DELTA,DELT1G,DELT2G,DELT3G,DELTUG,DELTHG,WENDEP,
     &                 ITST,KIT,DA,DB,JW)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION GS(MYG),G(MYG),THETA(MYG),F(MYG)
      DIMENSION DH(MYG),DH2(MYG),DH4(MYG),DH6(MYG),SDH(MYG)
      DIMENSION RVS1(MYG),RVS2(MYG),RVGS(MYG),RV(MYG),ETAV(MYG),
     &          ETAV2(MYG)
      DIMENSION VIKOEF(8),WENDEP(MYG)
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
      COMMON / CINT / GS,G,THETA,F
      COMMON / CRV / RVS1,RVS2,RVGS,ETAV,ETAV2,RV
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      COMMON / CST / DH,DH2,DH4,DH6,SDH
C
      WRITE (IW,134) 'GRENZSCHICHTDICKE BEI U=',GRENZ
  134 FORMAT(////,2X,A,F5.3)
C
      WRITE (IW,291) DELTA
  291 FORMAT (//,2X,'GRENZSCHICHTDICKE:',F13.9)
C
      WRITE (IW,292) DELT1G
  292 FORMAT (//,2X,'VERDRAENGUNGSDICKE:',F13.9)
C
      WRITE (IW,293) DELT2G
  293 FORMAT (//,2X,'IMPULSVERLUSTDICKE:',F13.9)
C
      WRITE (IW,294) DELT3G
  294 FORMAT (//,2X,'ENERGIEVERLUSTDICKE:',F13.9)
C
      WRITE (IW,295) DELTHG
  295 FORMAT (//,2X,'ENTHALPIEVERMEHRUNGSDICKE:',F13.9)
C
      WRITE (IW,296) DELTUG
  296 FORMAT (//,2X,'GESCHWINDIGKEITSVERLUSTDICKE:',F13.9,//)
C
      DO 11 I=1,JW
        WRITE (IW,297) I,WENDEP(I)
  297   FORMAT(2X,'VERALLGEMEINERTER WENDEPUNKT ',I2,' :',F13.9)
   11 CONTINUE
C
      WRITE (IW,105) KIT
  105 FORMAT (////,2X,'KONVERGENZ NACH',I7,' ITERATIONEN')
C
      FEHLER=ABS(DA)
      WRITE (IW,'(/,2X,A9,F13.10,//)') 'FEHLER = ',FEHLER
C
      IF (ICOM.EQ.1) THEN
        WRITE(IW,112) G(1)
  112   FORMAT (/,2X,'DIE ABLEITUNG DER DIMENSIONSLOSEN ENTHALPIE G(1) B
     &  ETRAEGT ',F12.6,///)
      ELSE
        WRITE(IW,111) THETA(1)
  111   FORMAT (/,2X,'DIMENSIONSLOSE ENTHALPIE THETA(1): ',F12.6,///)
      END IF
C
      IF(LSTR.EQ.2) THEN
      WRITE(IW,107) SDH(ITST)
  107 FORMAT(2X,'LETZTE BEI DER INTEGRATION VERWENDETE STUETZSTELLE:',
     &       F7.3)
      WRITE(IW,108) DH(ITST)
  108 FORMAT(2X,'BEI EINER MAXIMALEN SCHRITTWEITE VON',E13.5,///)
C
C *** BELEGUNG DER PHYS.FELDER FUER PUNKTE AUSSERHALB DER GRENZSCHICHT *
C
      DO 98 I=ITST+1,(NSTR-NSC+1)
        U(I)=U(I-1)
        TDL(I)=TDL(I-1)
        RHODL(I)=RHODL(I-1)
        DUDH(I)=DUDH(I-1)
        D2UDH(I)=D2UDH(I-1)
        DTDH(I)=DTDH(I-1)
        D2TDH(I)=D2TDH(I-1)
        DRHODH(I)=DRHODH(I-1)
        D2RDH(I)=D2RDH(I-1)
        DRMUE(I)=DRMUE(I-1)
        DMUEDT(I)=DMUEDT(I-1)
        D2MDT(I)=D2MDT(I-1)
        DRD2U(I)=DRD2U(I-1)
        PR(I)=PR(I-1)
        DRK(I)=DRK(I-1)
        DRKDT(I)=DRKDT(I-1)
        D2RKDT(I)=D2RKDT(I-1)
   98 CONTINUE
C
C *** BESETZUNG DER PHYS.FELDER AN D. NACH CHEBY. BERECHNETEN STELLEN
C
      DO 99 I=(NSC+1),(NSTR-NSC+1),NSC
        KKPCH=(I+NSC-1)/NSC
        U(KKPCH)=U(I)
        TDL(KKPCH)=TDL(I)
        RHODL(KKPCH)=RHODL(I)
        DUDH(KKPCH)=DUDH(I)
        D2UDH(KKPCH)=D2UDH(I)
        DTDH(KKPCH)=DTDH(I)
        D2TDH(KKPCH)=D2TDH(I)
        DRHODH(KKPCH)=DRHODH(I)
        D2RDH(KKPCH)=D2RDH(I)
        DRMUE(KKPCH)=DRMUE(I)
        DMUEDT(KKPCH)=DMUEDT(I)
        D2MDT(KKPCH)=D2MDT(I)
        DRD2U(KKPCH)=DRD2U(I)
        PR(KKPCH)=PR(I)
        DRK(KKPCH)=DRK(I)
        DRKDT(KKPCH)=DRKDT(I)
        D2RKDT(KKPCH)=D2RKDT(I)
   99 CONTINUE
      ITST=KKPCH
      END IF
C
      IF (LSTR.EQ.1.AND.ITSTC.LT.1000) THEN
        ITST=INT(ETA/DN)+1
C
        DO 19 I=ITSTC+1,ITST
          U(I)=U(I-1)
          TDL(I)=TDL(I-1)
          RHODL(I)=RHODL(I-1)
          DUDH(I)=DUDH(I-1)
          D2UDH(I)=D2UDH(I-1)
          DTDH(I)=DTDH(I-1)
          D2TDH(I)=D2TDH(I-1)
          DRHODH(I)=DRHODH(I-1)
          D2RDH(I)=D2RDH(I-1)
          DRMUE(I)=DRMUE(I-1)
          DMUEDT(I)=DMUEDT(I-1)
          D2MDT(I)=D2MDT(I-1)
          DRD2U(I)=DRD2U(I-1)
          PR(I)=PR(I-1)
          DRK(I)=DRK(I-1)
          DRKDT(I)=DRKDT(I-1)
          D2RKDT(I)=D2RKDT(I-1)
   19   CONTINUE
C
      END IF
C
C *** AUSDRUCKEN DER FELDER *******************************************
C
  140 FORMAT (4(2X,F14.9))
C
      WRITE (IW,150)
  150 FORMAT (/,2X,'LAENGSGESCHWINDIGKEIT U:',2X)
      WRITE (IW,140) (U(I),I=1,ITST)
C
      WRITE (IW,130)
  130 FORMAT (/,2X,'TEMPERATUR T:',2X)
      WRITE (IW,140) (TDL(I),I=1,ITST)
C
      WRITE (IW,170)
  170 FORMAT (/,2X,'DICHTE RHO',2X)
      WRITE (IW,140) (RHODL(I),I=1,ITST)
C
      WRITE (IW,210)
  210 FORMAT (/,2X,'RV:',2X)
      WRITE (IW,140) (RV(I),I=1,ITST)
C
      WRITE (IW,220)
  220 FORMAT (/,2X,'DUDH:')
      WRITE (IW,140) (DUDH(I),I=1,ITST)
C
      WRITE (IW,230)
  230 FORMAT (/,2X,'D2UDH:')
      WRITE (IW,140) (D2UDH(I),I=1,ITST)
C
      WRITE (IW,240)
  240 FORMAT (/,2X,'DTDH:')
      WRITE (IW,140) (DTDH(I),I=1,ITST)
C
      WRITE (IW,250)
  250 FORMAT (/,2X,'D2TDH:')
      WRITE (IW,140) (D2TDH(I),I=1,ITST)
C
      WRITE (IW,260)
  260 FORMAT (/,2X,'DRHODH:')
      WRITE (IW,140) (DRHODH(I),I=1,ITST)
C
      WRITE (IW,270)
  270 FORMAT (/,2X,'D2RDH:')
      WRITE (IW,140) (D2RDH(I),I=1,ITST)
C
      WRITE (IW,280)
  280 FORMAT (/,2X,'VISKOSITAET DRMUE:')
      WRITE (IW,140) (DRMUE(I),I=1,ITST)
C
      WRITE (IW,290)
  290 FORMAT (/,2X,'DMUEDT:')
      WRITE (IW,140) (DMUEDT(I),I=1,ITST)
C
      WRITE (IW,310)
  310 FORMAT (/,2X,'D2MDT:')
      WRITE (IW,140) (D2MDT(I),I=1,ITST)
C
      WRITE (IW,330)
  330 FORMAT (/,2X,'VERALLGEMEINERTER WENDEPUNKT:')
      WRITE (IW,140) (DRD2U(I),I=1,ITST)
C
      WRITE(IW,333)
  333 FORMAT(/,2X,'PRANDTLZAHL:')
      WRITE(IW,140) (PR(I),I=1,ITST)
C
      WRITE(IW,336)
  336 FORMAT(/,2X,'WAERMELEITFAEHIGKEIT:')
      WRITE(IW,140) (DRK(I),I=1,ITST)
C
      WRITE(IW,358)
  358 FORMAT (/,2X,'DRKDT:')
      WRITE(IW,140)  (DRKDT(I),I=1,ITST)
C
      WRITE(IW,359)
  359 FORMAT(/,2X,'D2RKDT:')
      WRITE (IW,140) (D2RKDT(I),I=1,ITST)
C
      IF (LSTR.EQ.1) THEN
        WRITE (IW,320)
  320   FORMAT (/,2X,'DETA:')
        WRITE (IW,140) (DETA(I),I=1,ITST)
C
        WRITE (IW,340)
  340   FORMAT (/,2X,'D2ETA:')
        WRITE (IW,140) (D2ETA(I),I=1,ITST)
      END IF
C
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM FELDAENDERUNG BEI SCHRITTWEITE GROESSER 0.1
C *********************************************************************
C
      SUBROUTINE FEAEND(ITST1,ITST,ITEIL,DN,DN1)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION RVS1(MYG),RVS2(MYG),RVGS(MYG),RV(MYG),ETAV(MYG),
     &          ETAV2(MYG)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      COMMON / CRV / RVS1,RVS2,RVGS,ETAV,ETAV2,RV
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
C
C *** TRAFO AUF URSRUENGLICHE ANZAHL DER STUETZSTELLEN ****************
C
      DO 10 I=1,ITST1
        J=ITEIL*(I-1)+1
        U(I)=U(J)
        TDL(I)=TDL(J)
        RHODL(I)=RHODL(J)
        FV(I)=FV(J)
        RV(I)=RV(J)
        DUDH(I)=DUDH(J)
        D2UDH(I)=D2UDH(J)
        DTDH(I)=DTDH(J)
        D2TDH(I)=D2TDH(J)
        DRHODH(I)=DRHODH(J)
        D2RDH(I)=D2RDH(J)
        DRMUE(I)=DRMUE(J)
        DMUEDT(I)=DMUEDT(J)
        D2MDT(I)=D2MDT(J)
        DRD2U(I)=DRD2U(J)
        PR(I)=PR(J)
        DRK(I)=DRK(J)
        DRKDT(I)=DRKDT(J)
        D2RKDT(I)=D2RKDT(J)
   10 CONTINUE
C
      DO 20 I=ITST1+1,MYG
        U(I)=0.
        TDL(I)=0.
        RHODL(I)=0.
        FV(I)=0.
        RV(I)=0.
        DUDH(I)=0.
        D2UDH(I)=0.
        DTDH(I)=0.
        D2TDH(I)=0.
        DRHODH(I)=0.
        D2RDH(I)=0.
        DRMUE(I)=0.
        DMUEDT(I)=0.
        D2MDT(I)=0.
        DRD2U(I)=0.
        PR(I)=0.
        DRK(I)=0.
        DRKDT(I)=0.
        D2RKDT(I)=0.
   20 CONTINUE
C
      ITST=ITST1
      DN=DN1
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG VON FV UND DICHTE
C *********************************************************************
C
      SUBROUTINE FVDI(DH6,ITST)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION RDLS1(MYG),RDLS2(MYG),RDLGS(MYG)
      DIMENSION US1(MYG),US2(MYG),UGS(MYG),DH6(MYG)
      DIMENSION TDLS1(MYG),TDLS2(MYG),TDLGS(MYG)
      COMMON / CURK / US1,US2,UGS
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
      COMMON / CTDLRK / TDLS1,TDLS2,TDLGS
      COMMON / CRHORK / RDLS1,RDLS2,RDLGS
C
C *** BERECHNUNG VON FV ***********************************************
C
      FV(1)=0.
      DO 10 I=2,ITST
        M=I-1
        FV(I)=FV(M)+DH6(M)*(U(M)+2.*US1(I)+2.*US2(I)+UGS(I))
   10 CONTINUE
C
C *** BERECHNUNG DER DICHTE ******************************************
C
      DO 30 I=1,ITST
        RHODL(I)=1./TDL(I)
        IF(I.EQ.1) THEN
C
C *** BESETZUNG AN NICHT DEFINIERTER STELLE **************************
C
          RDLS1(I)=1./TDL(I)
          RDLS2(I)=1./TDL(I)
          RDLGS(I)=1./TDL(I)
        ELSE
          RDLS1(I)=1./TDLS1(I)
          RDLS2(I)=1./TDLS2(I)
          RDLGS(I)=1./TDLGS(I)
        END IF
   30 CONTINUE
C
      END
C
C *********************************************************************
C     BERECHNUNG ENTHALPIE, TEMPERATUR UND VISKOITAET BEI KONSTANTER
C     PRANDTLZAHL UND KONSTANTER SPEZ. WAERMEKAPAZITAET
C *********************************************************************
C
      SUBROUTINE HTM0(H,T,TDL,RMUE,DRMUE,F,DUDH,THETA,PR)
      implicit real*8 (a-h,o-z)
      DIMENSION VIKOEF(8)
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
C
      H=H1+DH01H1*THETA
      T=H/CP+T0
      PR=PR0
      TDL=T/T1
      CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,T,RMUE)
      DRMUE=RMUE/RMUE1
      DUDH=F/DRMUE
C
      END
C
C *********************************************************************
C     BERECHNUNG ENTHALPIE, TEMPERATUR UND VISKOITAET BEI KONSTANTER
C     SPEZ. WAERMEKAPAZITAET UND VARIABLER PRANDTLZAHL
C *********************************************************************
C
      SUBROUTINE HTM1(H,T,TDL,RMUE,DRMUE,F,DUDH,THETA,PR)
      implicit real*8 (a-h,o-z)
      DIMENSION VIKOEF(8)
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
C
      H=H1+DH01H1*THETA
      T=H/CP+T0
      TDL=T/T1
      CALL PRV1(T,RMUE,PR)
      DRMUE=RMUE/RMUE1
      DUDH=F/DRMUE
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG VON ENTHALPIE,TEMPERATUR UND
C     VISKOSITAET FUER RUNGE KUTTA VERFAHREN
C *********************************************************************
C
      SUBROUTINE HTM2(H,T,TDL,RMUE,DRMUE,F,DUDH,THETA,PR,CPMI,M)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION TEMP(IMAX),CPDR(IMAX),C(IMAX-1,3),CPST(IMAX),
     &          CPST2(MYG),CPM(MYG)
      DIMENSION VIKOEF(8)
      COMMON / CCP / CPDR,CPST,TEMP,C,CPST2,CPM
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV

      H=H1+DH01H1*THETA
      IF((ICOM.EQ.1).AND.(M.EQ.1)) THEN
        T=TW
      ELSE
        T=H/CPMI+T0
      END IF
      CALL PRV2(T,C,IMAX,TEMP,CPST,C1,C2,C3,C4,VIKOEF,TV1,TV2,EMUE,PR,
     &          RMUE)
      TDL=T/T1
      DRMUE=RMUE/RMUE1
      DUDH=F/DRMUE
      END
C
C *********************************************************************
C     BERECHNUNG DES MITTLEREN WERTES VON CP DURCH
C     ELEMENTARE INTEGRATION DER KUBISCHEN SPLINES
C *********************************************************************
C
      SUBROUTINE KUBINT(TO,TU,C,IMAX,IINT,EPS,TEMP,CPST,CPMI)
      implicit real*8 (a-h,o-z)
      DIMENSION TEMP(IMAX),CPST(IMAX),C(IMAX-1,3)
C
C *** FESTLEGUNG DER OBEREN UND DER UNTEREN TEMPERATUR
C
      IF((TO-TU).GT.EPS) THEN
        T1=TU
        T2=TO
      ELSE IF((TU-TO).GT.EPS) THEN
        T1=TO
        T2=TU
      ELSE
        CALL CPINT(T0,C,IMAX,TEMP,CPST,CPMI)
        GOTO 9999
      END IF
C
C *** AUFSUCHEN DER UNTEREN UND DER OBEREN TEMPERATURSTUETZSTELLE *****
C
      DO 10 I=1,IMAX
        IF((TEMP(I).LE.T1).AND.(TEMP(I+1).GT.T1)) THEN
          N1=I
          GOTO 20
        END IF
   10 CONTINUE
   20 CONTINUE
      DO 50 I=1,IMAX
        IF((TEMP(I).LE.T2).AND.(TEMP(I+1).GT.T2)) THEN
          N2=I
          GOTO 60
        END IF
   50 CONTINUE
   60 CONTINUE
C
C *** BERECHNUNG FUER T1,T2 IM SELBEN TEMPERATURINTERVALL *************
C
      IF(N1.EQ.N2) THEN
        TINT1=T1-TEMP(N1)
        TINT2=T2-TEMP(N2)
        WINTEG=0.25*C(N1,3)*(TINT2*TINT2*TINT2*TINT2-TINT1*TINT1*TINT1
     &  *TINT1)+C(N1,2)/3.*(TIN2*TINT2*TINT2-TINT1*TINT1*TINT1)+0.5*
     &  C(N1,1)*(TINT2*TINT2-TINT1*TINT1)+CPST(N1)*(TINT2-TINT1)
      ELSE
C
C *** BERECHNUNG DURCH ABSCHNITTSWEISE INTEGRATION UEBER TEMPERATUR-
C     STUETZSTELLEN ***
C
        IF(N1.GE.IINT) THEN
          TINT1=50.
        ELSE
          TINT1=10.
        END IF
        TINT2=T1-TEMP(N1)
        WINTEG=0.25*C(N1,3)*(TINT1*TINT1*TINT1*TINT1-TINT2*TINT2*TINT2*
     &         TINT2)+C(N1,2)/3.*(TINT1*TINT1*TINT1-TINT2*TINT2*TINT2)+
     &         0.5*C(N1,1)*(TINT1*TINT1-TINT2*TINT2)+CPST(N1)*(TINT1-
     &         TINT2)
        DO 100 I=N1+1,N2-1
          IF(I.GE.IINT) THEN
            TINT=50.
          ELSE
            TINT=10.
          END IF
          WINTEG=WINTEG+0.25*C(I,3)*TINT*TINT*TINT*TINT+C(I,2)/3.*TINT*
     &          TINT*TINT+0.5*C(I,1)*TINT*TINT+CPST(I)*TINT
  100   CONTINUE
        TINT=T2-TEMP(N2)
        WINTEG=WINTEG+0.25*C(N2,3)*TINT*TINT*TINT*TINT+C(N2,2)/3.*TINT*
     &         TINT*TINT+0.5*C(N2,1)*TINT*TINT+CPST(N2)*TINT
      END IF
C
C *** DIVISION DES INTEGRALS DURCH TEMPERATURDIFFERENZ ERGIBT MITTLERES
C     CP ***
C
      CPMI=WINTEG/ABS(T2-T1)
 9999 CONTINUE
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON GS NACH ETA
C *********************************************************************
C
      SUBROUTINE FUGS(DGSDH,U,TDL)
      implicit real*8 (a-h,o-z)
      DGSDH=0.5*U/TDL
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON F NACH ETA
C *********************************************************************
C
      SUBROUTINE FUF(DFDH,GS,DRMUE,F)
      implicit real*8 (a-h,o-z)
      DFDH=-GS/DRMUE*F
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON THETA NACH ETA
C *********************************************************************
C
      SUBROUTINE FUTET(DTEDH,PR,DRMUE,G)
      implicit real*8 (a-h,o-z)
      DTEDH=PR/DRMUE*G
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON G NACH ETA
C *********************************************************************
C
      SUBROUTINE FUG(DGDH,PR,DRMUE,GS,G,DUDH)
      implicit real*8 (a-h,o-z)
      DGDH=-PR/DRMUE*GS*G-2.*DRMUE*DUDH*DUDH
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y1 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY1(DY1DH,DRMUE,Y5)
      implicit real*8 (a-h,o-z)
      DY1DH=Y5/DRMUE
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y2 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY2(DY2DH,DRMUE,Y6)
      implicit real*8 (a-h,o-z)
      DY2DH=Y6/DRMUE
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y3 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY3(DY3DH,PR,DRMUE,Y7)
      implicit real*8 (a-h,o-z)
      DY3DH=PR/DRMUE*Y7
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y4 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY4(DY4DH,PR,DRMUE,Y8)
      implicit real*8 (a-h,o-z)
      DY4DH=PR/DRMUE*Y8
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y5 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY5(DY5DH,GS,DRMUE,Y5)
      implicit real*8 (a-h,o-z)
      DY5DH=-GS/DRMUE*Y5
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y6 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY6(DY6DH,GS,DRMUE,Y6)
      implicit real*8 (a-h,o-z)
      DY6DH=-GS/DRMUE*Y6
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y7 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY7(DY7DH,DRMUE,F,PR,GS,Y5,Y7)
      implicit real*8 (a-h,o-z)
      DY7DH=-(4.*F*Y5+PR*GS*Y7)/DRMUE
      END
C
C *********************************************************************
C     BERECHNUNG DER ABLEITUNG VON Y8 NACH ETA
C *********************************************************************
C
      SUBROUTINE FUY8(DY8DH,DRMUE,F,PR,GS,Y6,Y8)
      implicit real*8 (a-h,o-z)
      DY8DH=-(4.*F*Y6+PR*GS*Y8)/DRMUE
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG KONSTANTER GROESSEN
C *********************************************************************
C
      SUBROUTINE KONST (C5)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION VIKOEF(8)
      DIMENSION TEMP(IMAX),CPDR(IMAX),C(IMAX-1,3),CPST(IMAX),
     &          CPST2(MYG),CPM(MYG)
      COMMON / CCP / CPDR,CPST,TEMP,C,CPST2,CPM
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
C
C *** BERECHNUNG DER ANSTROEMGESCHWINDIGKEIT *************************
C
      IF ((ICPPRV.EQ.1).OR.(ICPPRV.EQ.0)) THEN
        GAM=CP/CV
      ELSE
        CALL CPINT(T1,C,IMAX,TEMP,CPST,CPI)
        GAM=R/(CPI-R)+1.
      END IF
      CC=SQRT(R*1.E3*T1*GAM)
      U1=RM1*CC
C
C *** DEFINIERTE KONSTANTE ********************************************
C
      TV1=100.
      TV2=130.
      C1=1.458*1.E-6
      C2=0.6325*4.186*1.E-6
      C3=0.693873*1.E-7
      C4=0.222964*4.186*1.E-7
      C5=R/101.325
C
      VIKOEF(1)=-4.479148053679334E1
      VIKOEF(2)=3.195188079744342E2
      VIKOEF(3)=-9.716235566382709E2
      VIKOEF(4)=1.632645086771892E3
      VIKOEF(5)=-1.637375578884298E3
      VIKOEF(6)=9.802775658900685E2
      VIKOEF(7)=-3.234667180557399E2
      VIKOEF(8)=4.58157988617632E1
C
      EMUE=7.659704848E-6
C
C *** BERECHNUNG VON ANSTROEMDICHTE,ANSTROEMVISKOSITAET UND ANSTROEM-
C     ENTHALPIE BEI FESTLEGUNG VON H0 ALS STANDARDREAKTIONSENTHALPIE **
C
      IF((ICPPRV.EQ.1).OR.(ICPPRV.EQ.0)) THEN
        T01=T1+T1*(GAM-1.)/2.*RM1*RM1
        DH01H1=0.5*U1*U1*1.E-3
        H1=CP*(T1-T0)
      ELSE
        DH01H1=0.5*U1*U1*1.E-3
        CALL KUBINT(T1,T0,C,IMAX,IINT,EPS,TEMP,CPST,CPMI)
        H1=CPMI*(T1-T0)
      END IF
C
      RHO1=1./(C5*T1)
C
      CALL VISKOS (C1,C3,VIKOEF,TV1,TV2,EMUE,T1,RMUE1)
C
      CALL RKFT(T1,C2,C4,RK1)
C
C *** AUSDRRUCKEN DER WERTE KONSTANTER GROESSEN ***********************
C
      WRITE(IW,6)
    6 FORMAT(////,4X,'KONSTANTE GROESSEN:',//)
C
      WRITE(IW,11) U1
   11 FORMAT(/,7X,'ANSTROEMGESCHWINDIGKEIT:',F13.5,' M/S')
C
      WRITE(IW,8) CC
    8 FORMAT(/,7X,'SCHALLGESCHWINDIGKEIT:',F13.5,' M/S')
C
      WRITE(IW,4) RHO1
    4 FORMAT(/,7X,'ANSTROEMDICHTE:',F13.5,' KG/M/M/M')
C
      WRITE (IW,5) RMUE1
    5 FORMAT(/,7X,'VISKOSITAET DER AUSSENSTROEMUNG:'
     &       ,F13.7,' KG/M/S')
C
      WRITE (IW,9) RK1
    9 FORMAT(/,7X,'WAERMELEITFAEHIGKEIT DER AUSSENSTROEMUNG:'
     &       ,F13.7,' KJ/M/S')
C
      END
C
C ********************************************************************
C     LESEN DES INPUTS
C ********************************************************************
C
      SUBROUTINE LES(RMSCH,RSCH)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION TEMP(IMAX),CPDR(IMAX),C(IMAX-1,3),CPST(IMAX),
     &          CPST2(MYG),CPM(MYG)
      DIMENSION RSCH(IRSCH),RMSCH(IRSCH)
      dimension inf(100)
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      COMMON / CCP / CPDR,CPST,TEMP,C,CPST2,CPM
C
      WRITE(IW,1)
    1 FORMAT(//,5X,'NUMERISCHE INTEGRATION DER KOMPRESSIBLEN GRENZSCHICH
     &ICHTGLEICHUNGEN MIT EINSCHIESSVERFAHREN UND RUNGE KUTTA O(4)')
      WRITE (IW,90)
   90 FORMAT (/,5X,'RECHENGROESSEN',//)
C
      READ(IR,100) (INF(I),I=1,10)
      READ (IR,110) RM1,T1,CP,CV
      WRITE (IW,120) (INF(I),I=1,10)
      WRITE (IW,130) RM1,T1,CP,CV
C
      READ (IR,100) (INF(I),I=1,10)
      WRITE (IW,120) (INF(I),I=1,10)
      READ (IR,110) DN,ETA,ETAS,PR0
      WRITE (IW,130) DN,ETA,ETAS,PR0
C
      READ (IR,100) (INF(I),I=1,10)
      WRITE (IW,120) (INF(I),I=1,10)
      READ (IR,110) ASTR,YSTR,TW
      WRITE (IW,130) ASTR,YSTR,TW
C
      READ (IR,100) (INF(I),I=1,10)
      WRITE (IW,120) (INF(I),I=1,10)
      READ (IR,110) DUDHS,THETAS,GSCH
      WRITE (IW,130) DUDHS,THETAS,GSCH
C
      READ (IR,100) (INF(I),I=1,10)
      WRITE (IW,120) (INF(I),I=1,10)
      READ (IR,150) ICOM,ICPPRV,NSTR,IAUT,LSTR
      WRITE (IW,160) ICOM,ICPPRV,NSTR,IAUT,LSTR
C
  100 FORMAT (10A7)
  110 FORMAT (5F14.0)
  120 FORMAT (/,2X,5(2A7,1X))
  130 FORMAT (2X,5(F14.7,1X))
  150 FORMAT (5I14)
  160 FORMAT (2X,5(I14,1X))
C
      IF (ICOM.EQ.0) THEN
        WRITE(IW,155) 'BERECHNUNG BEI ADIABATER WAND'
  155   FORMAT(//,5X,A)
      ELSE
        WRITE(IW,156) 'BERECHNUNG ERFOLGT MIT KONSTANTER WANDTEMPERATUR'
  156   FORMAT(//,5X,A)
      END IF
C
      IF(ICPPRV.EQ.0) THEN
        WRITE(IW,158) 'CP UND PRANDTLZAHL KONSTANT'
  158   FORMAT(//,5X,A)
      ELSE IF(ICPPRV.EQ.1) THEN
        WRITE(IW,159) 'CP= KONSTANT,KAPPA=KAPPA(T),MUE=MUE(T)'
  159   FORMAT(//,5X,A,A)
      ELSE
        WRITE(IW,161) 'CP=CP(T),KAPPA=KAPPA(T),MUE=MUE(T),PR=PR(T)'
  161   FORMAT(//,5X,A)
        TEMP(1)=10.
        DO 200 I=2,IINT
          TEMP(I)=TEMP(I-1)+10.
  200   CONTINUE
        DO 300 I=IINT+1,IMAX
          TEMP(I)=TEMP(I-1)+50.
  300   CONTINUE
C
C *** EINLESEN DER DISKRETEN WERTE VON CP/R ***
C
        READ (IR,320) CPDR(1),CPDR(2),CPDR(3),CPDR(4)
  320   FORMAT (/,4F14.0)
        DO 400 I=4,108,4
          READ (IR,350) CPDR(I+1),CPDR(I+2),CPDR(I+3),CPDR(I+4)
  350     FORMAT (4F14.0)
  400   CONTINUE
        WRITE (IW,450) 'EINLESEN DER DISKRETEN WERTE VON'
        WRITE (IW,451) 'CP/R ALS FUNKTION DER TEMPERATUR'
  450   FORMAT (/,2X,A)
  451   FORMAT  (2X,A)
        DO 600 I=1,109,3
          WRITE (IW,500) TEMP(I),CPDR(I),TEMP(I+1),CPDR(I+1),TEMP(I+2),
     A                  CPDR(I+2)
  500     FORMAT (/,2X,'CP/R(',F6.1,')=',F6.4,5X,'CP/R(',F6.1,')=',F6.4,
     &            5X,'CP/R(',F6.1,')=',F6.4)
  600   CONTINUE
      END IF
C
C *** PROGRAMM SUCHT SELBSTSTAENDIG SCHAETZWERTE ***
C
      IF(IAUT.EQ.1) THEN
        CALL SCHATZ(ICOM,DUDHS,THETAS,GSCH,RSCH,IRSCH,RMSCH,RM1)
        WRITE(IW,'(//,2X,A)') 'SCHAETZWERTE FUER ANFANGSBEDINGUNGEN :'
        IF(ICOM.EQ.0) THEN
          WRITE(IW,'(/,2X,A,F14.7)') 'DUDHS  = ',DUDHS
          WRITE(IW,'(/,2X,A,F14.7,/)') 'THETAS = ',THETAS
        ELSE IF(ICOM.EQ.1) THEN
          WRITE(IW,'(/,2X,A,F14.7)') 'DUDHS  = ',DUDHS
          WRITE(IW,'(/,2X,A,F14.7,/)') 'GSCH   = ',GSCH
        END IF
      END IF
C
      END
C
C *********************************************************************
C    UNTERPROGRAMM ZUR BERECHNUNG DER STUETZSTELLEN BEI KONSTANTER
C    SCHRITTWEITE
C *********************************************************************
C
      SUBROUTINE LSTR0 (ETA,DN,ITST,ITEIL,ITST1,DN1)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION DH(MYG),DH2(MYG),DH4(MYG),DH6(MYG),SDH(MYG)
      COMMON / CST / DH,DH2,DH4,DH6,SDH
C
C *** SCHRITTWEITENVERKLEINERUNG BEI ZU GROSSER SCHRITTWEITE **********
C
      ITEIL=1
      DN1=DN
      ITST1=int(ETA/DN)+1
      DO 10 I=1,20
        TEIL=FLOAT(I)
        DN=DN/TEIL
        ITEIL=I
        IF(DN.LE.0.1) GOTO 20
   10 CONTINUE
   20 CONTINUE
C
      ITST=ITEIL*ITST1
C
      IF(ITST.GT.MYG)          STOP 'MYG'
C
      WRITE (IW,2) ITST1,DN1
    2 FORMAT(/,5X,I4,' STUETZSTELLEN BEI KONSTANTER SCHRITTWEITE DH=',
     &       F6.5)
C
      DO 15 I=1,ITST-1
        DH(I)=DN
        DH2(I)=DN/2.
        DH4(I)=DN/4.
        DH6(I)=DN/6.
   15 CONTINUE
C
      SDH(1)=0.
      DO 30 I=2,ITST
        SDH(I)=SDH(I-1)+DN
   30 CONTINUE
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM FUER STRECKUNG MIT ASINH
C *********************************************************************
C
      SUBROUTINE LSTR1 (ITST,ITSTC,DETA,D2ETA)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION DETA(MYG),D2ETA(MYG),SDH(MYG),DH(MYG),DH2(MYG),DH4(MYG),
     &          DH6(MYG)
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      COMMON / CST / DH,DH2,DH4,DH6,SDH
      dimension dummy(myg+1)
C
      ITST=int(ETA/DN)+1
C
      write(*,*)'Warnung: Aufruf von SBCUF mit dummy Feld'
      CALL SBCUF(DN,ETAS,ITST,ASTR,DETA,D2ETA,SDH,dummy)
      WRITE(IW,89)
   89 FORMAT(////,4X,'NICHTAEQUIDISTANTE STUETZSTELLEN DURCH',
     &               ' STRECKUNG MIT ASINH:',/)
      WRITE(IW,91) (SDH(I),I=1,ITST)
   91 FORMAT(/,4X,5F17.9)
      ITSTC=1000
C
      DO 14 I=1,ITST-1
        DH(I)=SDH(I+1)-SDH(I)
        DH2(I)=DH(I)/2.
        DH4(I)=DH(I)/4.
        DH6(I)=DH(I)/6.
        IF (DH(I).GT.DHMAX)   ITSTC=MIN(I,ITSTC)
   14 CONTINUE
C
      IF(ITSTC.LT.1000)       ITST=ITSTC
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG MIT SCHRITTWEITEN NACH
C     CHEBYSHEV'SCHER KOLL.M.
C *********************************************************************
C
      SUBROUTINE LSTR2(PI,ITST,NSC,ITSTC,NTSTR,NCHEB)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION DH(MYG),DH2(MYG),DH4(MYG),DH6(MYG),SDH(MYG)
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      COMMON / CST / DH,DH2,DH4,DH6,SDH
C
      NTSTR=2*NSTR-1
C
C *** TRANSFORMATION AUF BERECHNUNGSBEREICH ***************************
C
      NCHEB=NSTR
      NSC=5
      NSTR=NSC*NSTR
      NTSTR=NSC*NTSTR
C
C *** BERECHNUNG DER KOLLOKATIONSPUNKTE NACH CHEBYSHEV-APPROXIMATION *
C
      DO 25 I=1,NSTR
        DH(I)=COS(PI*(I-1)/NTSTR)
   25 CONTINUE
C
C *** BERECHNUNG DER AKTUELLEN SCHRITTWEITE AN DER STELLE I ***********
C
      DO 27 I=1,NSTR-NSC
        DH(I)=-YSTR*(DLOG(DH(I+1))-DLOG(DH(I)))
   27 CONTINUE
C
C *** BERECHNUNG DER KOLLOKATIONSPUNKTE IM TRANSFORMIERTEN BEREICH ****
C
      SDH(1)=0.
      DO 28 I=2,(NSTR-NSC+1)
        SDH(I)=DH(I-1)+SDH(I-1)
   28 CONTINUE
C
      WRITE (IW,3)
    3 FORMAT(/,5X,'BERECHNUNG ERFOLGT NACH CHEBYSHEV')
      WRITE(IW,9)
    9 FORMAT(////,4X,' STUETZSTELLEN NACH CHEBYSHEV APPROXIMATION:',/)
      WRITE(IW,91) (SDH(I),I=1,(NSTR-NSC+1),NSC)
   91 FORMAT(/,4X,5F17.9)
C
C *** ITST IST DIE ANZAHL DER STUETZSTELLEN AN DENEN DIE PROFILE
C     BERECHNET WERDEN ************************************************
C
      ITST=NSTR-NSC+1
      ITSTC=1000
C
      DO 26 I=1,ITST
        DH2(I)=DH(I)/2.
        DH4(I)=DH(I)/4.
        DH6(I)=DH(I)/6.
        IF(DH(I).GT.DHMAX)     ITSTC=MIN(I,ITSTC)
   26 CONTINUE
C
      IF(ITSTC.LT.1000)      ITST=ITSTC
      WRITE(IW,7)
    7 FORMAT(////,4X,' INTERVALLE NACH CHEBYSHEV APPROXIMATION:',/)
      WRITE(IW,92) (DH(I),I=1,NSTR-NSC)
   92 FORMAT(/,4X,5F17.9)
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG NEUER CPM UND NEUER ANFANGS-
C     BEDINGUNGEN
C *********************************************************************
C
      SUBROUTINE NEUAB(CPMS1,CPMS2,CPMGS,ITEND,ITST,DA,DB,TDLS1,
     &                 TDLS2,TDLGS)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION TEMP(IMAX),CPDR(IMAX),C(IMAX-1,3),CPST(IMAX),
     &          CPST2(MYG)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION RKY1(MYG),RKY2(MYG),RKY3(MYG),RKY4(MYG),RKY5(MYG),
     &          RKY6(MYG),RKY7(MYG),RKY8(MYG)
      DIMENSION GS(MYG),G(MYG),THETA(MYG),F(MYG)
      DIMENSION TDLS1(MYG),TDLS2(MYG),TDLGS(MYG)
      DIMENSION CPMS1(MYG),CPMS2(MYG),CPMGS(MYG),CPM(MYG)
      DIMENSION VIKOEF(8)
      COMMON / CINT / GS,G,THETA,F
C     COMMON / CTDLRK / TDLS1,TDLS2,TDLGS
      COMMON / CY / RKY1,RKY2,RKY3,RKY4,RKY5,RKY6,RKY7,RKY8
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      COMMON / CCP / CPDR,CPST,TEMP,C,CPST2,CPM
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
C
C *** BESTIMMUNG VON NEUEM MITTLEREM CP DURCH TEMPERATURSCHRITTE DES
C     VORIGEN ITERATIONSSCHRITTS **************************************
C
      IF(ICPPRV.EQ.2) THEN
        DO 555 J=1,ITST
          T=TDL(J)*T1
          TS1=TDLS1(J)*T1
          TS2=TDLS2(J)*T1
          TGS=TDLGS(J)*T1
          CALL KUBINT(T,T0,C,IMAX,IINT,EPS,TEMP,CPST,CPM(J))
          IF(J.GE.2) THEN
            CALL KUBINT(TS1,T0,C,IMAX,IINT,EPS,TEMP,CPST,CPMS1(J))
            CALL KUBINT(TS2,T0,C,IMAX,IINT,EPS,TEMP,CPST,CPMS2(J))
            CALL KUBINT(TGS,T0,C,IMAX,IINT,EPS,TEMP,CPST,CPMGS(J))
          END IF
  555   CONTINUE
      END IF
C
C *** BERECHNUNG NEUER ANFANGSBEDINGUNGEN *****************************
C
      ITEND=0
      DB=-(THETA(ITST)*RKY1(ITST)+RKY3(ITST)*(1.-U(ITST)))/
     &       (RKY4(ITST)*RKY1(ITST)-RKY3(ITST)*RKY2(ITST))
      DA=(1.-U(ITST))/RKY1(ITST)-RKY2(ITST)/RKY1(ITST)*DB
      IF((ABS(DA).LE.EPS).AND.(ABS(DB).LE.EPS))              ITEND=1
C
C *** ANFANGSBEDINGUNGEN FUER ADIABATE WAND ***************************
C
      IF(ICOM.EQ.0) THEN
        F(1)=F(1)+DA
        THETA(1)=THETA(1)+DB
C
C *** ANFANGSBEDINGUNGEN FUER KONSTANTE WANDTEMPERATUR ***************
C
      ELSE
        F(1)=F(1)+DA
        IF(ITEND.EQ.0) THEN
          G(1)=G(1)+DB
        END IF
C
        IF((ICPPRV.EQ.1).OR.(ICPPRV.EQ.0)) THEN
          THETA(1)=CP*(TW-T1)/DH01H1
        ELSE
          CALL KUBINT(TW,T1,C,IMAX,IINT,EPS,TEMP,CPST,CPMI)
          THETA(1)=CPMI*(TW-T1)/DH01H1
        END IF
C
      END IF
C
      END
C
C *********************************************************************
C     BERECHNUNG DER PRANDTLZAHL FUER CP KONSTANT,KAPPA
C     UND MUE VARIABEL
C *********************************************************************
C
      SUBROUTINE PRV1(T,RMUE,PRFT)
      implicit real*8 (a-h,o-z)
      DIMENSION VIKOEF(8)
      COMMON / CKON1/ C1,C2,C3,C4,DH01H1,H1,RMUE1,RHO1,RK1,R,VIKOEF
      COMMON / CKON2/ EMUE,TV1,TV2,T0,U1,CC,GRENZ
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
C
      CALL RKFT(T,C2,C4,RK)
      CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,T,RMUE)
C
      PRFT=CP*RMUE/RK
C
      END
C
C *********************************************************************
C     BERECHNUNG DER PRANDTLZAHL FUER KAPPA,CP UND MUE
C     ALS FUNKTIONEN DER TEMPERATUR
C *********************************************************************
C
      SUBROUTINE PRV2(T,C,IMAX,TEMP,CPST,C1,C2,C3,C4,VIKOEF,TV1,
     &                TV2,EMUE,PRV,RMUE)
      implicit real*8 (a-h,o-z)
      DIMENSION TEMP(IMAX),CPST(IMAX),C(IMAX-1,3)
      DIMENSION VIKOEF(8)
C
      CALL RKFT(T,C2,C4,RK)
      CALL VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,T,RMUE)
C
      CALL CPINT(T,C,IMAX,TEMP,CPST,CPI)
      PRV=CPI*RMUE/RK
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG DER WAERMELEITFAEHIGKEIT ODER
C     DEREN ABLEITUNGEN NACH DER TEMPERATUR BEI KONSTANTER
C     PRANDTLZAHL UND KONSTANTER SPEZ. WAERMEKAPAZITAET
C *********************************************************************
C
      SUBROUTINE RKFMUE(RMUE,RK)
      implicit real*8 (a-h,o-z)
      COMMON / CLESR / RM1,T1,CP,CV,DN,ETA,ETAS,PR0,ASTR,YSTR,TW,DUDHS,
     &                 THETAS,GSCH
      COMMON / CLESI / NSTR,LSTR,ICOM,ICPPRV
      RK=CP/PR0*RMUE
      END
C
C *********************************************************************
C     BERECHNUNG VON KAPPA IN ABHAENGIGKEIT VON DER TEMPERATUR
C *********************************************************************
C
      SUBROUTINE RKFT(T,C2,C4,RKVAR)
C
C *** BESTIMMUNG VON KAPPA IN VERSCHIEDENEN TEMPERATURBEREICHEN *******
C
      implicit real*8 (a-h,o-z)
      IF(T.GT.80.) THEN
        RKVAR=C2*SQRT(T)/(1.+245.4/T*10**(-12./T))
      ELSE
        RKVAR=C4*T
      END IF
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM RUNGE KUTTA 1.,2. HALBSCHRITT UND 1.VOLLSCHRITT
C *********************************************************************
C
      SUBROUTINE RK12G(Y,YM,D,A)
      implicit real*8 (a-h,o-z)
      Y=YM+D*A
      END
C
C *********************************************************************
C     UNTERPROGRAMM RUNGE KUTTA SCHRITT Y(N+1)
C *********************************************************************
C
      SUBROUTINE RKYNP1(Y,YM,D,AM,A1,A2,AG)
      implicit real*8 (a-h,o-z)
      Y=YM+D*(AM+2.*(A1+A2)+AG)*1.D0
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR ERZEUGUNG DES TRANSFORMATIONSGITTERS BEI
C     NICHTKONSTANTER SCHRITTWEITE
C *********************************************************************
C
      SUBROUTINE SBCUF(DY,ETAS,MY,B,DETA,D2ETA,YK,eta)
      implicit real*8 (a-h,o-z)
      DIMENSION ETA(MY),DETA(MY),D2ETA(MY),YK(MY)
      MG= int(ETAS/DY)
C
      DO 10 J=1,MG
        ETA(J)=FLOAT(J-1)*DY
        YK(J)=FLOAT(J-1)*DY
        DETA(J)=1.0
        D2ETA(J)=0.0
  10  CONTINUE
C
      DO 20 J=MG+1,MY
        ETA(J)=FLOAT(J-1)*DY
        YY=FLOAT(J-MG)*DY
        YK(J)=FLOAT(MG-1)*DY+1./B*SINH(YY*B)
        DETA(J)=1./SQRT(B*B*YY*YY+1.)
        D2ETA(J)=-B*B*YY/(B*B*YY*YY+1.)/SQRT(B*B*YY*YY+1.)
  20  CONTINUE
C
      END
C
C ********************************************************************
C     SCHATZWERTE FUER ANFANGSBEDINGUNGEN BERECHNEN
C ********************************************************************
C
      SUBROUTINE SCHATZ(ICOM,DUDHS,THETAS,GSCH,RSCH,IRSCH,RMSCH,RM1)
      implicit real*8 (a-h,o-z)
      DIMENSION RSCH(IRSCH),RMSCH(IRSCH)
C
      IF(ICOM.EQ.0) THEN
        THETAS=0.846
      ELSE
        GSCH=0.319
      END IF
C
      DO 10 I=1,IRSCH
        IF(RM1.LT.RMSCH(I)) THEN
          DUDHS=RSCH(I)
          GOTO 20
        END IF
   10 CONTINUE
   20 CONTINUE
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG DER MATRIX DER SPLINE-KOEFFIZIENTEN
C *********************************************************************
C
      SUBROUTINE SPKOEF (TEMP,CPST,IMAX,C)
      implicit real*8 (a-h,o-z)
      DIMENSION TEMP(IMAX),CPST(IMAX),C(IMAX-1,3)
      CALL ICSCCU(TEMP,CPST,IMAX,C,IMAX-1,IER)
      IF(IER.NE.0) THEN
        PRINT*,IER
        STOP 'IER'
      END IF
      END
c ********************************************************************
      subroutine icsccu
c     dummy - fuehrt zum Programmabbruch wenn aufgerufen
c     ist im original eine bibliotheksroutine von IMSL
      stop'ICSCCU'
      end
C
C ********************************************************************
C     BERECHNUNG DER VISKOSITAET
C ********************************************************************
C
      SUBROUTINE VISKOS(C1,C3,VIKOEF,TV1,TV2,EMUE,T,RMUE)
      implicit real*8 (a-h,o-z)
      DIMENSION VIKOEF(8)
C
C *** NORMIERUNG DER TEMPERATUR MIT 110.4K ***
C
      TN=T/110.4
C
      IF(T.GT.TV2) THEN
        RMUE=C1*T*SQRT(T)/(T+110.4)
      ELSE IF (T.GT.TV1) THEN
        RMUE=EMUE*(VIKOEF(1)*TN**7+VIKOEF(2)*TN**6+VIKOEF(3)*TN**5+
     &       VIKOEF(4)*TN**4+VIKOEF(5)*TN**3+VIKOEF(6)*TN*TN+
     &       VIKOEF(7)*TN+VIKOEF(8))
      ELSE
        RMUE=C3*T
      END IF
C
      END
C
C *********************************************************************
C     UNTERPROGRAMM ZUR BERECHNUNG DES VERALLGEMEINERTEN WENDEPUNKTS
C *********************************************************************
C
      SUBROUTINE WENDE(ITST,ICOM,WENDEP,J)
      implicit real*8 (a-h,o-z)
      PARAMETER(MYG=1000,ITMAX=30,EPS=1.E-9,IRSCH=7,DHMAX=0.5,IR=5,
     &          IW=6,IMAX=112,IINT=80)
      DIMENSION U(MYG),TDL(MYG),RHODL(MYG),FV(MYG),DUDH(MYG),DTDH(MYG),
     &          D2UDH(MYG),D2TDH(MYG),DRHODH(MYG),D2RDH(MYG),DETA(MYG),
     &          DRMUE(MYG),DMUEDT(MYG),D2MDT(MYG),D2ETA(MYG),
     &          DRD2U(MYG),PR(MYG),DRK(MYG),DRKDT(MYG),D2RKDT(MYG)
      DIMENSION DH(MYG),DH2(MYG),DH4(MYG),DH6(MYG),SDH(MYG),WENDEP(MYG)
      COMMON / COUT / U,TDL,RHODL,FV,DUDH,D2UDH,DTDH,D2TDH,DRHODH,
     &                D2RDH,DRMUE,DMUEDT,D2MDT,DRD2U,PR,DRK,DRKDT,
     &                D2RKDT,DETA,D2ETA
      COMMON / CST / DH,DH2,DH4,DH6,SDH
C
      DO 100 I=1,ITST
        DRD2U(I)=DRHODH(I)*DUDH(I)+RHODL(I)*D2UDH(I)
  100 CONTINUE
C
C *** BERECHNUNG DER STUETZSTELLEN ************************************
C
      J=1
      DO 200 I=2,ITST-1
        IF(((DRD2U(I).GT.EPS).AND.(DRD2U(I+1).LT.EPS)).OR.((DRD2U(I).
     &    LT.EPS).AND.(DRD2U(I+1).GT.EPS))) THEN
C
C *** BESTIMMUNG DES WENDEPUNKTS DURCH LINEARE INTERPOLATION **********
C
          WENDEP(J)=SDH(I)+DRD2U(I)/(DRD2U(I)-DRD2U(I+1))*(SDH(I+1)-
     &              SDH(I))
          J=J+1
C
        END IF
  200 CONTINUE
      J=J-1
C
      END
