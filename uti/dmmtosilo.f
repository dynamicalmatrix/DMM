      Program DMMTOSILO
C     This program elaborates the data calculated by micro3D for
C     producing the required silo/visit graphs.
C
C.......................................................................
C
C     UNITS:
C     L --> cm
C     M --> g
C     T --> s
C
C     MAGN. FIELD --> Oe
C.......................................................................
C
C     REQUIRED EXTERNAL LIBRARIES: silo
C     (https://wci.llnl.gov/simulation/computer-codes/silo)
C
C.......................................................................

      use string_utility

      implicit NONE
 
      include "/usr/local/include/silo.inc"

C     .. Parameters ..

C     ------------------------------------------------------------------
      real ZERO, ONE
      complex CZERO, CI
      parameter (ZERO=0.E0, ONE=1.E0, CZERO=(0.E0,0.E0), CI=(0.E0,1.E0))
C     MAXNMAT: Maximum number of different materials
      integer MAXNMAT
      parameter (MAXNMAT=10)
C     ------------------------------------------------------------------

C     .. Variables ..

C     ------------------------------------------------------------------
C     The complex Eigenvalues
      real Eigvalr
C     ------------------------------------------------------------------
C     The Eigenvectors
      real, allocatable :: Eigvect(:,:)
C     ------------------------------------------------------------------
C     Orientation of the equilibrium magnetization 
      real, allocatable :: Theta0(:), Phi0(:)
C     ------------------------------------------------------------------
C     Arrays for idendifying the active (magnetic) elements in the part
C     and keeping track of their indices
      integer(8), allocatable :: Active(:), Iactive(:)
c     integer Active(XCMAX*YCMAX*ZCMAX), Iactive(XCMAX*YCMAX*ZCMAX)
C     Total number of active cells and of variables
      integer(8) Cactive
      integer Nactive
C     ------------------------------------------------------------------
      integer Sel, Pqual, unicelx, unicely
      real H0, H0angle, H0polar
      real, allocatable :: Exchange(:,:), Ms(:), Gamma(:)
      real, allocatable :: K2p(:), Anisang(:,:)
      real PolAngle, IncAngle, ScatAngle
      real, allocatable :: Freq(:), CrossSect(:), NewCrossSect(:)
      real Thickover, fmin, fmax
      complex, allocatable :: Eps(:), KMO(:)
      logical Array, unformatted
      logical, allocatable ::  check(:)
      real  KBloch(2)
C     ------------------------------------------------------------------
C     Misc variables
      logical Err
      integer dbid, ierr, iierr, optlistid, optlists(6)
      integer dims(3), dX, dY
      integer jj, n, l, ll, Isterr, Line, Lgood, lx, ly
      integer ii, i, j, p, base, far, infilelength
      integer ROW, COLUMN, LAYER
      integer Ccount, Xborder, Yborder
      integer(8) NumMat, Goodmodes
      integer Degeneracy, Savedmodes
      integer(4) Xcount, Ycount, Zcount
      integer argomenti(3), ipid, istat, iretpid
      integer types(6), lnames(6), ldefs(6)
      integer, allocatable :: Iindex(:)
      integer(8), allocatable :: MatMap(:)
      integer, allocatable :: matnos(:)
      integer, allocatable :: matlist(:)
      real, allocatable :: x(:), y(:), z(:)
      real, allocatable :: MxR(:),MxI(:),MyR(:),MyI(:),MzR(:),MzI(:)
      real PI, xr
      real, allocatable :: Fvar(:), Xvar(:)
      real, allocatable :: XpsF(:), YpsF(:), ZpsF(:)
      real, allocatable :: XpsH(:), YpsH(:), ZpsH(:)
      real Xcenter, Ycenter
      real Temp(9), Energy(1)
c     real Cpmin, Cpmax
      real Minfreq, Maxfreq, rad
      complex Ctemp(3), Ctemp2(3)
      character*80 Filename, OldFilename, Infilename, maininfilename
      character*30 comandi(3)
      character*6 names(6), defs(6)
      character*12 forms
      character*12, allocatable :: MatName(:)
      complex Mx, My, Mz, expfactor
      complex, allocatable :: Pv(:,:), Ep(:,:), Em(:,:), A(:,:,:)
      complex, allocatable :: PP(:,:)
      character*11 fnformat1, fnformat2
      data types/6*DB_VARTYPE_SCALAR/
      data lnames/6*6/
      data ldefs/6*4/
      data optlists/6*DB_F77NULL/
C     ------------------------------------------------------------------

C     .. commons ..
C     ------------------------------------------------------------------
      real Cellsize, ZCellsize
      COMMON/PGEOM/Cellsize, ZCellsize
C     ------------------------------------------------------------------
C     .. on-line functions ..
C
C     We assume the following reference frame for the part:
C           .
C        y / \
C           |
C           | Ycount
C           | . |
C           | . |
C           | . |
C           |
C           | 5 |
C           | 4 |
C           | 3 |  cell-index
C           | 2 |
C           | 1 | 12345...
C           | r \-------------
C           | \ c 12345... Xcount
C           -----------------------------> x
C
C     The z-axis is perpendicular to xy and points outward.
C     The position of an element is: 
C     (x,y,z)=(cellsize c, cellsize r, cellsize p).
C     In the first layer (p=1) we label each cell of the grid with an
C     index which runs from 1 to Xcount in the first row, from 1+Xcount 
C     to 2*Xcount in the second row, and so on up to Xcount*Ycount.
C     On the next layer(s) the count continues from Xcount*Ycount+1
C     with the same rule.
C     Therefore: cell-index=c+Xcount*(r-1)+Xcount*Ycount*(p-1). Reversing 
C     this equation we find the row and the column indices as a function 
C     of cell-index:
      LAYER(i)=(i-1)/(Xcount*Ycount)+1
      ROW(i)=MOD(i-1,Xcount*Ycount)/Xcount+1
      COLUMN(i)=MOD(MOD(i-1,Xcount*Ycount),Xcount)+1
C     ------------------------------------------------------------------
C     Note: our reference frame coincides with that of oommf. 
C     ------------------------------------------------------------------
      RAD(xr)=(xr*PI)/180

C     .. Executable Statements ..

      PI=4.*ATAN(ONE)
C     ------------------------------------------------------------------
C     Get arguments
      CALL GET_COMMAND_ARGUMENT (1,maininfilename,infilelength,ierr)
      if (ierr.ne.0) then
        print *, 'Wrong progam invocation. Syntax: '//
     2  'dmmtosilo input-file-name'
        stop
      endif
C     ------------------------------------------------------------------
C     Data input 
      open(UNIT=2,FILE=maininfilename,STATUS='OLD',ERR=22,IOSTAT=Isterr)
      Line=1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11) Infilename, unformatted
C     read some data from Infilename
      if (index(Infilename,'.dmm').EQ.0) then
        print *, 'The main input file must have ''.dmm'' extension.'
        stop
      endif
      if (unformatted) then
        forms='UNFORMATTED'
      else
        forms='FORMATTED'
      endif
      open(UNIT=1,FILE=Infilename,STATUS='UNKNOWN',FORM=forms,
     2 ERR=23,IOSTAT=Isterr)
      if (unformatted) then
        read(1,ERR=53) NumMat
      else
        read(1,*,ERR=53) NumMat
      endif
      if (NumMat.gt.MAXNMAT) then
        print *, 'dmmtosilo: the material number exceed its maximum'
        print *, 'Increase MAXNUMMAT to ', NumMat, 'and recompile'
        stop
      endif
      allocate (MatName(NumMat))
      allocate (matnos(NumMat+1))
      allocate (Ms(NumMat))
      allocate (Gamma(NumMat))
      allocate (K2p(NumMat))
      allocate (Eps(NumMat))
      allocate (KMO(NumMat))
      allocate (Anisang(NumMat,2))
      allocate (Exchange(NumMat,NumMat))
      allocate (check(NumMat+2))
      allocate (Pv(NumMat,0:3))
      allocate (Ep(NumMat,0:3))
      allocate (Em(NumMat,0:3))
      allocate (A(NumMat,0:3,2))
      allocate (PP(NumMat,0:3))
      if (unformatted) then
        read(1,ERR=53) Cellsize,ZCellsize,
     2  Xcount,Ycount,Zcount,MatName,Ms,Gamma,Exchange,H0,H0angle,
     3  H0polar,K2p,Anisang,Cactive,fmin,fmax,Xcenter,Ycenter,
     4  Array,KBloch
        if (unicelx.ne.1.or.unicely.ne.1) then
          print *, 'Warning, unformatted and array repetition '
     2    //'not working yet'
        endif
      else
        read(1,*,ERR=53) Cellsize,ZCellsize,
     2  Xcount,Ycount,Zcount,MatName,Ms,Gamma,Exchange,H0,H0angle,
     3  H0polar,K2p,Anisang,Cactive,fmin,fmax,Xcenter,Ycenter,
     4  Array,KBloch
      endif
C     end read some data from Infilename
C     continue reading data input
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11) Sel, Pqual
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11) Minfreq, Maxfreq
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11)
      Line=Line+1
      read(2,*,ERR=11) unicelx, unicely
      close(2)
C     ------------------------------------------------------------------
C     Conversions
      PolAngle=RAD(PolAngle)
      IncAngle=RAD(IncAngle)
      ScatAngle=RAD(ScatAngle)
      Minfreq=Minfreq*1.E09*2*PI
      Maxfreq=Maxfreq*1.E09*2*PI
C     ------------------------------------------------------------------
C     Sanity tests
      Err=.FALSE.
      if(Sel.lt.1.or.Sel.gt.3) then
        print *, 'Datain: ',
     2  '"Sel" out of range [1-3]: what do you want to do?'
        Err=.TRUE.
      endif
      if((Pqual.lt.1.or.Pqual.gt.2).and.Sel.ge.2) then
        print *, 'Datain: ',
     2  '"Pqual" out of range [0-2]: what do you want to do?'
        Err=.TRUE.
      endif
      if(Sel.ge.4) then
        if(IncAngle.lt.-PI.or.IncAngle.gt.PI) then
          print *, 'Datain: the Incidence Angle seems unreasonable!'
          Err=.TRUE.
        endif
        if(ScatAngle.lt.-PI.or.ScatAngle.gt.PI) then
          print *, 'Datain: the ScatAngle Angle seems unreasonable!'
          Err=.TRUE.
        endif
        if (ScatAngle.ne.-IncAngle) then
          print *, 'Warning: you are NOT in a backscattering geometry.'
        endif
        if(Thickover.lt.0) then
          print *, 'Datain: Overlayer thickness is negative!'
          Err=.TRUE.
        endif
        if(Thickover.gt.100.E-04) then
          print *, 'Datain: ',
     2    'Overlayer thickness is larger than 100 micron!'
          Err=.TRUE.
        endif
      endif
      if(MinFreq.lt.0.or.MinFreq.gt.1.E12) then
        print *, 'Datain: the requested min frequency in plots seems '
        print *, 'unreasonable!'
        Err=.TRUE.
      endif
      if(MaxFreq.lt.0.or.MaxFreq.gt.1.E13.or.MaxFreq.lt.Minfreq) then
        print *, 'Datain: the requested max frequency in plots seems '
        print *, 'unreasonable!'
        Err=.TRUE.
      endif
      if (unicelx.lt.1.or.unicelx.gt.10.or.unicely.lt.1.or.
     2 unicely.gt.10) then
        print *, 'Datain: the requested representation of periodic '
     2  //'array seems unreasonable!'
        Err=.TRUE.
      endif
      if (Err) then
        print *
        print *, 'Stopping due to input error(s)...'
        stop
      endif
C     continue reading from the main input file and do conversions and tests
      Ccount=Xcount*Ycount*Zcount
      Nactive=2*Cactive
      dX=0
      dY=0
      if(Minfreq.lt.fmin) Minfreq=fmin
      if(Maxfreq.gt.fmax) Maxfreq=fmax
      do i=1,NumMat+1
        matnos(i)=i-1
      enddo
C     ------------------------------------------------------------------
      if (unformatted) then
        read(1) Goodmodes
      else
        read(1,*) Goodmodes
      endif
      print *, 'goodmodes', goodmodes
      allocate (Active(0:Ccount))
      allocate (IActive(Cactive))
      allocate (MatMap(Cactive))
      allocate (Theta0(Cactive))
      allocate (Phi0(Cactive))
      if (unformatted) then
        do i=1, Ccount
          read(1) Active(i)
        enddo
        Active(0)=0
        do i=1, Cactive
          read(1) Iactive(i), MatMap(i), Theta0(i), Phi0(i)
        enddo
      else
        do i=1, Ccount
          read(1,*) Active(i)
        enddo
        Active(0)=0
        do i=1, Cactive
          read(1,*) Iactive(i), MatMap(i), Theta0(i), Phi0(i)
        enddo
      endif
C     ------------------------------------------------------------------
      allocate (Eigvect(Nactive,2))
      allocate (Freq(MAX(Nactive,1000)))
      allocate (CrossSect(Nactive))
      allocate (NewCrossSect(MAX(Nactive,1000)))
      allocate (Iindex(Nactive))
      allocate (Xvar(Nactive))
      allocate (XpsF(Nactive))
      allocate (YpsF(Nactive))
      allocate (ZpsF(Nactive))
      allocate (XpsH(Nactive))
      allocate (YpsH(Nactive))
      allocate (ZpsH(Nactive))
      allocate (Fvar(Nactive))
      allocate (MxR((Xcount+dX)*(Ycount+dY)*Zcount*unicelx*unicely))
      allocate (MxI((Xcount+dX)*(Ycount+dY)*Zcount*unicelx*unicely))
      allocate (MyR((Xcount+dX)*(Ycount+dY)*Zcount*unicelx*unicely))
      allocate (MyI((Xcount+dX)*(Ycount+dY)*Zcount*unicelx*unicely))
      allocate (MzR((Xcount+dX)*(Ycount+dY)*Zcount*unicelx*unicely))
      allocate (MzI((Xcount+dX)*(Ycount+dY)*Zcount*unicelx*unicely))
      allocate (matlist((Xcount+dX)*(Ycount+dY)*Zcount*unicelx*unicely))
      allocate (x((Xcount+dX)*unicelx+1))
      allocate (y((Ycount+dY)*unicely+1))
      allocate (z(Zcount+1))
C     No errors; prepare for printing the results
      do i=1, Xcount+1
        x(i)=(i-1)
      enddo
      do i=1,Ycount+1
        y(i)=(i-1)
      enddo
      do i=1, Zcount+1
        z(i)=(i-1)
      enddo
C     ------------------------------------------------------------------
C     Preparation for profiles calculation
C     ------------------------------------------------------------------
C     ifort
!     call PXFACCESS('./output',8,0,Info)
!     if (Info.eq.0) then
C       exists, remove it
!       comandi(1)='rm'
!       argomenti(1)=0
!       comandi(2)='-rf '
!       argomenti(2)=0
!       comandi(3)='./output'
!       argomenti(3)=0
!       CALL PXFFORK(ipid,Info)
!       if (Info .ne. 0) then
!          print *,'FAILED: PXFFORK call with error = ',Info
!       else
!         if (ipid .eq. 0) then
!            print *,'CHILD1: execing rm'
!            call PXFEXECVP ('rm',0,comandi,argomenti,3,Info)
!            print *,'FAILED: PXFEXEC call with error = ',Info
!            stop
!         else
!            print *,'PARENT: waiting for CHILD1'
!            CALL PXFWAIT(istat,iretpid,Info)
!            print *,'PARENT: CHILD1 finished'
!         endif
!       endif
!     endif
c     (in any case) does not exist, create it
!     call PXFMKDIR ('./output',8,488,Info)
!     comandi(1)='cp'
!     argomenti(1)=0
!     comandi(2)=Maininfilename
!     argomenti(2)=0
!     comandi(3)='output/in.m3s'
!     argomenti(3)=0
!     CALL PXFFORK(ipid,Info)
!     if (Info .ne. 0) then
!        print *,'FAILED: PXFFORK call with error = ',Info
!     else
!       if (ipid .eq. 0) then
!          print *,'CHILD1: execing cp'
!          call PXFEXECVP ('cp',0,comandi,argomenti,3,Info)
!          print *,'FAILED: PXFEXEC call with error = ',Info
!          stop
!       else
!          print *,'PARENT: waiting for CHILD1'
!          CALL PXFWAIT(istat,iretpid,Info)
!          print *,'PARENT: CHILD1 finished'
!       endif
!     endif
C     ------------------------------------------------------------------
C     Gfortran
      call system('rm -rf output')
      call system('mkdir output')
      call system('cp '//Maininfilename//' output')
C     ------------------------------------------------------------------
      open(UNIT=2,FILE='output/mat.tab')
      write (2,*) 0, 'non-magnetic'
      do i=1, NumMat
        write (2,*) i, MatName(i)
      enddo
      close(2)
C     ------------------------------------------------------------------
C     Print the ground state
      iierr = dbcreate('output/static.silo', 18, DB_CLOBBER, 
     2      DB_LOCAL, "informazioni", 12, DB_PDB, dbid)
      ierr = dbmkoptlist(9, optlistid)
      ierr = dbaddcopt  (optlistid, DBOPT_XLABEL, "X Axis", 6)
      ierr = dbaddcopt  (optlistid, DBOPT_YLABEL, "Y Axis", 6)
      ierr = dbaddcopt  (optlistid, DBOPT_ZLABEL, "Z Axis", 6)
      ierr = dbaddcopt  (optlistid, DBOPT_XUNITS, "cells", 5)
      ierr = dbaddcopt  (optlistid, DBOPT_YUNITS, "cells", 5)
      ierr = dbaddcopt  (optlistid, DBOPT_ZUNITS, "cells", 5)
      dims(1)=Xcount+1
      dims(2)=Ycount+1
      dims(3)=Zcount+1
      iierr = dbputqm (dbid, "mesh3D", 6, "xcoords", 7, "ycoords", 7,
     2 "zcoords", 7, x, y, z, dims, 3,
     3 DB_FLOAT, DB_COLLINEAR, optlistid, ierr)
      do p=1,Zcount
        do i=1,Ycount
          do n=1,Xcount
            ii= n+Xcount*(i-1)+Xcount*Ycount*(p-1)
            if(Active(ii).eq.0) then
              MxR(ii)=ZERO
              MyR(ii)=ZERO
              MzR(ii)=ZERO
            else
              do l=1,Cactive
                ll=Iactive(l)
                if(ii.eq.ll) then
                  MxR(ii)=SIN(Phi0(l))*COS(Theta0(l))
                  MyR(ii)=SIN(Phi0(l))*SIN(Theta0(l))
                  MzR(ii)=COS(Phi0(l))
                  goto 41
                endif
              enddo
              print *, 'Micro3D: this should never happen!'
              stop
  41        endif
          enddo
        enddo
      enddo
      dims(1)=Xcount
      dims(2)=Ycount
      dims(3)=Zcount
      iierr=dbputqv1(dbid, "Msx", 3, "mesh3D", 6, MxR, dims, 3, 
     2 DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
      iierr=dbputqv1(dbid, "Msy", 3, "mesh3D", 6, MyR, dims, 3, 
     2 DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
      iierr=dbputqv1(dbid, "Msz", 3, "mesh3D", 6, MzR, dims, 3, 
     2 DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
      iierr = dbputdefvars(dbid, "defvars", 7, 1, "Ms", 2, 
     2 DB_VARTYPE_VECTOR, "{Msx,Msy,Msz}", 13, DB_F77NULL, ierr)
      iierr=dbclose(dbid)
      Xborder=(Xcount+dX)*unicelx
      Yborder=(Ycount+dY)*unicely
      do i=1, (Xcount+dX)*unicelx+1
        x(i)=(i-1)
      enddo
      do i=1,(Ycount+dY)*unicely+1
        y(i)=(i-1)
      enddo
      do i=1, Zcount+1
        z(i)=(i-1)
      enddo
C     ------------------------------------------------------------------
C     Perform the cycle on eigenvalues/eigenvector
      Degeneracy=1
      Lgood=0
      do j=1,Goodmodes
C       Read the eigenvalues
        if (unformatted) then
          read(1,END=20) Eigvalr
          read(1) Energy(1)
        else
          read(1,*,END=20) Eigvalr
          read(1,*) Energy(1)
        endif
        if (Energy(1).eq.0) then
          Energy(1)=ONE
        endif
        print *, 'Frequency: ', Eigvalr/2/PI/1.E09, ' GHz'
C       Read the eigenvectors
        if (unformatted) then
          do i=1,Nactive
            read(1) Eigvect(i,1), Eigvect(i,1+1)
          enddo
        else
          do i=1,Nactive
            read(1,*) Eigvect(i,1), Eigvect(i,1+1)
          enddo
        endif
        if (Eigvalr.gt.Maxfreq.or.
     2  Eigvalr.lt.Minfreq.or.Sel.le.1) goto 31
        Lgood=Lgood+1
C       ----------------------------------------------------------------
C       Various calculation (profiles, power spectrum)
C       ----------------------------------------------------------------
C       Profiles
        do p=1,Zcount
C         Calculate each layer separately
          do i=1,Ycount+dY
            do n=1,Xcount+dX
              if(i.le.Ycount.and.n.le.Xcount) then
                jj= n+Xcount*(i-1)+Xcount*Ycount*(p-1)
              else
                jj=0
              endif
              base=(p-1)*Xborder*Yborder+(i-1)*Xborder+n
              if(Active(jj).eq.0) then
c               (n,i,p) --> (p-1)*Xcount*Ycount+(i-1)*Xcount+n
                MxR(base)=ZERO
                MxI(base)=ZERO
                MyR(base)=ZERO
                MyI(base)=ZERO
                MzR(base)=ZERO
                MzI(base)=ZERO
                matlist(base)=0
              else
                do l=1,Cactive
                  ll=Iactive(l)
                  if(jj.eq.ll) then
                    matlist(base)=MatMap(l)
C                   Perpendicular (z) component:
                    MzR(base)=-Eigvect(2*l,1)*SIN(Phi0(l))
                    MzI(base)=-Eigvect(2*l,1+1)*SIN(Phi0(l))
C                   Parallel (x) component:
                    MxR(base)=-SIN(Phi0(l))*SIN(Theta0(l))*
     2              Eigvect(2*l-1,1)+COS(Phi0(l))*
     3              COS(Theta0(l))*Eigvect(2*l,1)
                    MxI(base)=-SIN(Phi0(l))*SIN(Theta0(l))*
     2              Eigvect(2*l-1,1+1)+COS(Phi0(l))*
     3              COS(Theta0(l))*Eigvect(2*l,1+1)
C                   Parallel (y, along the applied field) component:
                    MyR(base)=+SIN(Phi0(l))*COS(Theta0(l))*
     2              Eigvect(2*l-1,1)+COS(Phi0(l))*
     3              SIN(Theta0(l))*Eigvect(2*l,1)
                    MyI(base)=+SIN(Phi0(l))*COS(Theta0(l))*
     2              Eigvect(2*l-1,1+1)+COS(Phi0(l))*
     3              SIN(Theta0(l))*Eigvect(2*l,1+1)
                    goto 42
                  endif
                enddo
                print *, 'Micro3D: this should never happen!'
                stop
  42          endif
c             print *, p, i, n, base, matlist(base)
            enddo
          enddo
        enddo
C       Periodic array representation: define M in other supercells
        do lx=1, unicelx
          do ly=1, unicely
            if (lx.gt.1.or.ly.gt.1) then
              do p=1,Zcount
                do i=1,Ycount+dY
                  do n=1,Xcount+dX
                    base=(p-1)*Xborder*Yborder+(i-1)*Xborder+n
                    far=(p-1)*Xborder*Yborder+(i-1+(ly-1)*(Ycount+dY))*
     2              Xborder+n+(lx-1)*(Xcount+dX)
                    expfactor=exp(2*CI*PI*(KBloch(1)*(lx-1)+
     2              KBloch(2)*(ly-1)))
                    Mx=CMPLX(MxR(base),MxI(base))*expfactor
                    My=CMPLX(MyR(base),MyI(base))*expfactor
                    Mz=CMPLX(MzR(base),MzI(base))*expfactor
                    MxR(far)=REAL(Mx)
                    MxI(far)=IMAG(Mx)
                    MyR(far)=REAL(My)
                    MyI(far)=IMAG(My)
                    MzR(far)=REAL(Mz)
                    MzI(far)=IMAG(Mz)
                    matlist(far)=matlist(base)
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
C       ----------------------------------------------------------------
        Temp(1)=Eigvalr/2/PI/1.E09
        if (Temp(1).lt.10) then
          write(Filename,'(F5.3,5H.silo)') Temp(1)
          line=10
        else if (Temp(1).lt.100) then
          write(Filename,'(F6.3,5H.silo)') Temp(1)
          line=11
        else if (Temp(1).lt.1000) then
          write(Filename,'(F7.3,5H.silo)') Temp(1)
          line=12
        else if (Temp(1).lt.10000) then
          write(Filename,'(F8.3,5H.silo)') Temp(1)
        endif
        if (OldFilename.EQ.Filename) then
          Degeneracy=Degeneracy+1
          if (Temp(1).lt.10) then
            fnformat1='(F5.3,1H.,'
            line=12
          else if (Temp(1).lt.100) then
            fnformat1='(F6.3,1H.,'
            line=13
          else if (Temp(1).lt.1000) then
            fnformat1='(F7.3,1H.,'
            line=14
          else if (Temp(1).lt.10000) then
            fnformat1='(F8.3,1H.,'
            line=15
          endif
          if (Degeneracy.lt.10) then
            fnformat2='I1,5H.silo)'
          else if (Degeneracy.lt.100) then
            fnformat2='I2,5H.silo)'
            line=line+1
          else if (Degeneracy.lt.1000) then
            fnformat2='I3,5H.silo)'
            line=line+2
          endif
          write(Filename,fnformat1//fnformat2) Temp(1), Degeneracy
        else
          Degeneracy=1
          OldFilename=Filename
        endif
        iierr = dbcreate('output/'//Filename, 7+line, DB_CLOBBER, 
     2        DB_LOCAL, "informazioni", 12, DB_PDB, dbid)
        ierr = dbmkoptlist(9, optlistid)
        ierr = dbaddcopt  (optlistid, DBOPT_XLABEL, "X Axis", 6)
        ierr = dbaddcopt  (optlistid, DBOPT_YLABEL, "Y Axis", 6)
        ierr = dbaddcopt  (optlistid, DBOPT_ZLABEL, "Z Axis", 6)
        ierr = dbaddcopt  (optlistid, DBOPT_XUNITS, "cells", 5)
        ierr = dbaddcopt  (optlistid, DBOPT_YUNITS, "cells", 5)
        ierr = dbaddcopt  (optlistid, DBOPT_ZUNITS, "cells", 5)
C       nodes:
        dims(1)=(Xcount+dX)*unicelx+1
        dims(2)=(Ycount+dY)*unicely+1
        dims(3)=Zcount+1
        iierr = dbputqm (dbid, "mesh3D", 6, "xcoords", 7, "ycoords", 7,
     2  "zcoords", 7, x, y, z, dims, 3,
     3  DB_FLOAT, DB_COLLINEAR, optlistid, ierr)
C     ZONES:
        dims(1)=(Xcount+dX)*unicelx
        dims(2)=(Ycount+dY)*unicely
        dims(3)=Zcount
        iierr=dbputqv1(dbid, "MxR", 3, "mesh3D", 6, mxr, dims, 3, 
     2  DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
        iierr=dbputqv1(dbid, "MxI", 3, "mesh3D", 6, mxi, dims, 3, 
     2  DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
        iierr=dbputqv1(dbid, "MyR", 3, "mesh3D", 6, myr, dims, 3, 
     2  DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
        iierr=dbputqv1(dbid, "MyI", 3, "mesh3D", 6, myi, dims, 3, 
     2  DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
        iierr=dbputqv1(dbid, "MzR", 3, "mesh3D", 6, mzr, dims, 3, 
     2  DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
        iierr=dbputqv1(dbid, "MzI", 3, "mesh3D", 6, mzi, dims, 3, 
     2  DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, optlistid, ierr)
        iierr=dbputqv1(dbid, "mask", 4, "mesh3D", 6, matlist, dims, 3, 
     2  DB_F77NULL, 0, DB_INT, DB_ZONECENT, optlistid, ierr)
        iierr=dbputmat(dbid, "Material", 8, "mesh3D", 6, 
     2  int(NumMat+1,4), 
     3  matnos, matlist, dims, 3, DB_F77NULL, DB_F77NULL, DB_F77NULL, 
     4  DB_F77NULL, 0, DB_FLOAT, optlistid, ierr)
        names(1)="invMxR"
        names(1)="invMxR"
        names(2)="invMyR"
        names(3)="invMzR"
        names(4)="invMxI"
        names(5)="invMyI"
        names(6)="invMzI"
        defs(1)="-MxR  "
        defs(2)="-MyR  "
        defs(3)="-MzR  "
        defs(4)="-MxI  "
        defs(5)="-MyI  "
        defs(6)="-MzI  "
        iierr = dbset2dstrlen(6)
        iierr = dbputdefvars(dbid, "defvars", 7, 6, names, lnames, 
     2  types, defs, ldefs, optlists, ierr)
        iierr=dbclose(dbid)
C       ----------------------------------------------------------------
C       Power spectrum
        if (Sel.le.2) goto 31
        do i=1,3
          Ctemp(i)=CZERO
          Ctemp2(i)=CZERO
        enddo
        do i=1,Cactive
          ii=Iactive(i)
C         --------------------------------------------------------------
C         Perpendicular (z) power spectrum:
          Ctemp(3)=Ctemp(3)-CMPLX(Eigvect(2*i,1),Eigvect(2*i,1+1))*
     2    SIN(Phi0(i))
C         Parallel (x) power spectrum:
          Ctemp(1)=Ctemp(1)-SIN(Phi0(i))*SIN(Theta0(i))*
     2    CMPLX(Eigvect(2*i-1,1),Eigvect(2*i-1,1+1))+COS(Phi0(i))*
     3    COS(Theta0(i))*CMPLX(Eigvect(2*i,1),Eigvect(2*i,1+1))
C         Parallel (y) power spectrum:
          Ctemp(2)=Ctemp(2)+SIN(Phi0(i))*COS(Theta0(i))*
     2    CMPLX(Eigvect(2*i-1,1),Eigvect(2*i-1,1+1))+COS(Phi0(i))*
     3    SIN(Theta0(i))*CMPLX(Eigvect(2*i,1),Eigvect(2*i,1+1))
C         --------------------------------------------------------------
          if (COLUMN(int(Iactive(i))).lt.Xcenter) then
C           Perpendicular (z) power spectrum:
            Ctemp2(3)=Ctemp2(3)-CMPLX(Eigvect(2*i,1),Eigvect(2*i,1+1))*
     2      SIN(Phi0(i))
C           Parallel (x) power spectrum:
            Ctemp2(1)=Ctemp2(1)-SIN(Phi0(i))*SIN(Theta0(i))*
     2      CMPLX(Eigvect(2*i-1,1),Eigvect(2*i-1,1+1))+COS(Phi0(i))*
     3      COS(Theta0(i))*CMPLX(Eigvect(2*i,1),Eigvect(2*i,1+1))
C           Parallel (y) power spectrum:
            Ctemp2(2)=Ctemp2(2)+SIN(Phi0(i))*COS(Theta0(i))*
     2      CMPLX(Eigvect(2*i-1,1),Eigvect(2*i-1,1+1))+COS(Phi0(i))*
     3      SIN(Theta0(i))*CMPLX(Eigvect(2*i,1),Eigvect(2*i,1+1))
          endif
        enddo
        Xvar(Lgood)= Eigvalr/2/PI/1.E09
        XpsF(Lgood)= ABS(Ctemp(1)/Cactive)
        YpsF(Lgood)= ABS(Ctemp(2)/Cactive)
        ZpsF(Lgood)= ABS(Ctemp(3)/Cactive)
        XpsH(Lgood)= ABS(Ctemp2(1)/Cactive)
        YpsH(Lgood)= ABS(Ctemp2(2)/Cactive)
        ZpsH(Lgood)= ABS(Ctemp2(3)/Cactive)
C       To get the normalized power spectrum change to:
C       ABS(Ctemp/Cactive)/SQRT(Energy(1))
C       ABS(Ctemp2/Cactive)/SQRT(Energy(1))
  31  enddo
  20  Savedmodes=Lgood-1
      close(1)
      deallocate (MxR)
      deallocate (Mxi)
      deallocate (MyR)
      deallocate (MyI)
      deallocate (MzR)
      deallocate (MzI)
      deallocate (matlist)
      deallocate (x)
      deallocate (y)
      deallocate (z)
C     ------------------------------------------------------------------
C     Postprocessing
C     ------------------------------------------------------------------
C     Power spectrum
      if (Sel.le.2) stop
      iierr = dbcreate('output/powsp.silo', 17, DB_CLOBBER, 
     2      DB_LOCAL, "informazioni", 12, DB_PDB, dbid)
      ierr = dbmkoptlist(4, optlistid)
      ierr = dbaddcopt  (optlistid, DBOPT_XLABEL, "Frequency", 9)
      ierr = dbaddcopt  (optlistid, DBOPT_YLABEL, "Power spectrum", 14)
      ierr = dbaddcopt  (optlistid, DBOPT_XUNITS, "GHz", 3)
      ierr = dbaddcopt  (optlistid, DBOPT_YUNITS, "arb. units", 10)
      iierr=dbputcurve(dbid, "FulldotX", 8, Xvar, XpsF, DB_FLOAT,
     2 Lgood, optlistid, ierr)
      iierr=dbputcurve(dbid, "FulldotY", 8, Xvar, YpsF, DB_FLOAT,
     2 Lgood, optlistid, ierr)
      iierr=dbputcurve(dbid, "FulldotZ", 8, Xvar, ZpsF, DB_FLOAT,
     2 Lgood, optlistid, ierr)
      iierr=dbputcurve(dbid, "HalfdotX", 8, Xvar, XpsH, DB_FLOAT,
     2 Lgood, optlistid, ierr)
      iierr=dbputcurve(dbid, "HalfdotY", 8, Xvar, YpsH, DB_FLOAT,
     2 Lgood, optlistid, ierr)
      iierr=dbputcurve(dbid, "HalfdotZ", 8, Xvar, ZpsH, DB_FLOAT,
     2 Lgood, optlistid, ierr)
      iierr=dbclose(dbid)
      deallocate(Xvar)
      deallocate(XpsF)
      deallocate(YpsF)
      deallocate(ZpsF)
      deallocate(XpsH)
      deallocate(YpsH)
      deallocate(ZpsH)
      stop
  11  print *, 'dmmtosilo: in.m3s: error reading line ', Line
      stop
  22  print *, 'dmmtosilo: cannot open in.m3s, status:', Isterr
      stop
  23  print *, 'dmmtosilo: cannot open main input file (dmm)',
     2         ', status:', Isterr
      stop
  53  print *
      print *, 'dmmtosilo: error reading file: ', Infilename
      stop
      end
