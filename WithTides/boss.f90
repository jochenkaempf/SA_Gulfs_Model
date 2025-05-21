
!-------------------------------------------------------------------------------
! General Sigma Coordinate Model  - GSCM Version 1.0
!-------------------------------------------------------------------------------
! By:
!   Jochen Kaempf
!     School of Chemistry, Physics and Earth Sciences, Flinders University Adelaide
! 
!-------------------------------------------------------------------------------
PROGRAM RUNA
  !-------------------------------------
  USE param
  USE functions
  USE cohini
  USE cohrun
  INTEGER(4) :: NTTT, NTTOT, NTOUT,NCOUT,NTRA
  CHARACTER(3) :: OFILE
  INTEGER(4) :: IOUT,JOUT,KOUT, k,j,i,n
  REAL :: depp, zposreal(ntrac), dista(ntrac)
  INTEGER(4) :: YEAROLD
  LOGICAL :: PARTI
! LAGRANGIAN MODEL JK
   REAL :: DX2,DY2,epsx,epsy,term1,term2,term3,term4
  INTEGER(4) :: II,JJ,WEST,SOUTH

! 
 ! ********** INITILIZE OCEAN MODEL *********

idayold = 0

PARTI = .TRUE.
DO n = 1,ntrac
  trac_start(n) = .true.
END DO

metmo = 0

CALL INITCOH

 ! time step is delt, simulation time:

NTTOT = 8*365*24*3600/INT(delt)

 ! snapshot outputs

NTOUT = 365*24*3600/INT(delt)/12

 ! float outputs

NTRA = 24*3600/INT(delt)

time = 0.

IF(IOUTT>0) CALL ANALYS
IF(IOUTF.EQ.1) CALL INTEGR0 

IOUT = 0
JOUT = 0
KOUT = 0

CALL IDATES(IDATE,IMONT,IDAY,IHOUR,IMIN) 

MONTAV = IMONT
IGLOB = 0
NMEAN = 0
DAYAV = IDAY
ICOUNT = 0

Yearold  = 2000

! set initial heat fluxes to zero
CALL zero33(outSW,1,nr,nc)
CALL zero33(outLW,1,nr,nc)
CALL zero33(outSEN,1,nr,nc)
CALL zero33(outLAT,1,nr,nc)
CALL zero33(outEVA,1,nr,nc)
CALL zero33(outPRE,1,nr,nc)

QUNIT = 987

OPEN(QUNIT,file = 'dat/qflux.dat',form='formatted',status='unknown')

DO n = 1,nreg

QTNIT1(n) = 550+n
QTNIT2(n) = 650+n
QTNIT3(n) = 750+n
QTNIT4(n) = 850+n

END DO

!========================================================================
IF(PARTI)THEN
! --particle module
OPEN(QTNIT1(1),file = 'dat/1floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(1),file = 'dat/1floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(1),file = 'dat/1floatz.dat',form='formatted',status='unknown',recl=1000000)

OPEN(QTNIT1(2),file = 'dat/2floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(2),file = 'dat/2floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(2),file = 'dat/2floatz.dat',form='formatted',status='unknown',recl=1000000)

OPEN(QTNIT1(3),file = 'dat/3floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(3),file = 'dat/3floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(3),file = 'dat/3floatz.dat',form='formatted',status='unknown',recl=1000000)

OPEN(QTNIT1(4),file = 'dat/4floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(4),file = 'dat/4floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(4),file = 'dat/4floatz.dat',form='formatted',status='unknown',recl=1000000)

OPEN(QTNIT1(5),file = 'dat/5floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(5),file = 'dat/5floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(5),file = 'dat/5floatz.dat',form='formatted',status='unknown',recl=1000000)

OPEN(QTNIT1(6),file = 'dat/6floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(6),file = 'dat/6floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(6),file = 'dat/6floatz.dat',form='formatted',status='unknown',recl=1000000)

OPEN(QTNIT1(7),file = 'dat/7floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(7),file = 'dat/7floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(7),file = 'dat/7floatz.dat',form='formatted',status='unknown',recl=1000000)

OPEN(QTNIT1(8),file = 'dat/8floatx.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT2(8),file = 'dat/8floaty.dat',form='formatted',status='unknown',recl=1000000)
OPEN(QTNIT3(8),file = 'dat/8floatz.dat',form='formatted',status='unknown',recl=1000000)

END IF

!========================================================================
! read warm-up file

  open(unit=97,file='dat/input.dat',form='unformatted',status='unknown')
  read(97) t,s,ro,rpress,conc1,conc2,age
  read(97) u2,v2,w2,w2phys,zeta2
  read(97) ud2,vd2
  close(unit=97)

DO i = 1,nc
DO j = 1,nr
DO k = 1,nz
 age(k,j,i) = 0.0
 DO n = 1,nreg
  conc(n,k,j,i) = 0.0
 END DO
END DO
END DO
END DO

 ! ********** TIME ITERATION *********

  DO NTTT = 1, NTTOT

time = time + delt

! write(*,*)'TIME ',time/(24.*3600.)

! call dynamics module

  CALL RUNCOH(NTTT)

!===================


! particle output
     IF ((MOD(NTTT,NTRA)==0).AND.PARTI) THEN 

DO n = 1,ntrac

  II = IT(n)
  JJ = JT(n)
  KK = KT(n)
  xpos(n) = REAL(II-0.5)*gx2(1,1)+XT(n)
  ypos(n) = REAL(JJ-0.5)*gy2(1)+YT(n)
  DX2 = 0.5*GX2(JJ,II)
  DY2 = 0.5*GY2(JJ)
  DZ2 = 0.5*GZSC(KK,JJ,II)

   WEST = 0
  epsx = XT(n)/GX2(JJ,II)
  IF(XT(n)<0.0)THEN
     WEST = -1
     epsx = 1.0+XT(n)/GX2(JJ,II-1)
  END IF

  SOUTH = 0
  epsy = YT(n)/GY2(JJ)
  IF(YT(n)<0.0)THEN
     SOUTH = -1
     epsy = 1.0+YT(n)/GY2(JJ-1)
  END IF

  term1 = h2atc(JJ+SOUTH,II+WEST)*(1.0-epsx)*(1.0-epsy)
  term2 = h2atc(JJ+SOUTH,II+WEST+1)*epsx*(1.0-epsy)
  term3 = h2atc(JJ+SOUTH+1,II+WEST)*(1.0-epsx)*epsy
  term4 = h2atc(JJ+SOUTH+1,II+WEST+1)*epsx*epsy
  depp = term1+term2+term3+term4
  zpos(n) = (1.0-(gz0(KK,JJ,II)-DZ2+ZT(n)) )*depp
  zposreal(n) = zpos(n)
  posH(n) = depp
END DO

     WRITE(QTNIT1(1),'(1000F12.6)')(xpos(n)/1000.0,n=1,1000)
     WRITE(QTNIT2(1),'(1000F12.6)')(ypos(n)/1000.0,n=1,1000)
     WRITE(QTNIT3(1),'(1000F14.8)')(zposreal(n),n=1,1000)

     WRITE(QTNIT1(2),'(1000F12.6)')(xpos(n)/1000.0,n=1001,2000)
     WRITE(QTNIT2(2),'(1000F12.6)')(ypos(n)/1000.0,n=1001,2000)
     WRITE(QTNIT3(2),'(1000F14.8)')(zposreal(n),n=1001,2000)

     WRITE(QTNIT1(3),'(1000F12.6)')(xpos(n)/1000.0,n=2001,3000)
     WRITE(QTNIT2(3),'(1000F12.6)')(ypos(n)/1000.0,n=2001,3000)
     WRITE(QTNIT3(3),'(1000F14.8)')(zposreal(n),n=2001,3000)

     WRITE(QTNIT1(4),'(1000F12.6)')(xpos(n)/1000.0,n=3001,4000)
     WRITE(QTNIT2(4),'(1000F12.6)')(ypos(n)/1000.0,n=3001,4000)
     WRITE(QTNIT3(4),'(1000F14.8)')(zposreal(n),n=3001,4000)

     WRITE(QTNIT1(5),'(1000F12.6)')(xpos(n)/1000.0,n=4001,5000)
     WRITE(QTNIT2(5),'(1000F12.6)')(ypos(n)/1000.0,n=4001,5000)
     WRITE(QTNIT3(5),'(1000F14.8)')(zposreal(n),n=4001,5000)

     WRITE(QTNIT1(6),'(1000F12.6)')(xpos(n)/1000.0,n=5001,6000)
     WRITE(QTNIT2(6),'(1000F12.6)')(ypos(n)/1000.0,n=5001,6000)
     WRITE(QTNIT3(6),'(1000F14.8)')(zposreal(n),n=5001,6000)

     WRITE(QTNIT1(7),'(1000F12.6)')(xpos(n)/1000.0,n=6001,7000)
     WRITE(QTNIT2(7),'(1000F12.6)')(ypos(n)/1000.0,n=6001,7000)
     WRITE(QTNIT3(7),'(1000F14.8)')(zposreal(n),n=6001,7000)

     WRITE(QTNIT1(8),'(1000F12.6)')(xpos(n)/1000.0,n=7001,8000)
     WRITE(QTNIT2(8),'(1000F12.6)')(ypos(n)/1000.0,n=7001,8000)
     WRITE(QTNIT3(8),'(1000F14.8)')(zposreal(n),n=7001,8000)

     END IF

! output of spatial distribution arrays
!     IF (MOD(NTTT,NTOUT)==0) THEN 

IF(IYEAR/=YEAROLD)THEN

         YEAROLD=IYEAR         

         IOUT = IOUT + 1
         jout = iout/10
         kout = iout-jout*10
         OFILE = '.00'
         IF (iout<10) THEN
         WRITE(OFILE(3:3),'(I1.1)')iout
         ELSE
         WRITE(OFILE(2:2),'(I1.1)')jout
         WRITE(OFILE(3:3),'(I1.1)')kout
         ENDIF
         WRITE(*,*)'OUTPUT TO FILE =',OFILE,'AT TIME',time/(24.*3600.)
 
         open(unit=97,file='dat/ageb'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(200F16.8)')(age(1,j,i)/(24.*3600),i=1,nc)
         END DO
         close(unit=97)

         open(unit=97,file='dat/ages'//OFILE,form='formatted',status='unknown')
         DO j = 1,nr
         write(97,FMT='(200F16.8)')(age(nz,j,i)/(24.*3600),i=1,nc)
         END DO
         close(unit=97)

  open(unit=97,file='dat/output'//OFILE,form='unformatted',status='unknown')
  write(97) t,s,ro,rpress
  write(97) u2,v2,w2,w2phys,zeta2
  write(97) ud2,vd2
  close(unit=97)

     ENDIF
!END IF
  
  END DO

close(unit=33)


END PROGRAM RUNA
!-------------------------------------------------------------------------------

