      program Micro3d
C     This program calculates the spin wave frequency and dynamic
C     magnetization profile in a dot by a micromagnetic calculation,
C     using a 3D model.
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
C     REQUIRED LIBRARIES: LAPACK/MKL (Intel)
C
C.......................................................................

      use input_parameters

      implicit NONE

C     .. Parameters ..
      real ZERO, ONE, EPS
      complex CZERO, CI
C     EPS:   Roundoff error affecting the demagnetization tensor
      parameter (EPS=1.E-06)
      parameter (ZERO=0.E0, ONE=1.E0, CZERO=(0.E0,0.E0), CI=(0.E0,1.E0))
C     MAXR: Max. number of R to be considered in the sums for interparticle
C           coupling
C     Reps: Relative error tolerated on the sums over R in the interparticle
C           coupling
      integer MAXR
      real REPS
      parameter (REPS=0.002)
      data MAXR/100/
C     ------------------------------------------------------------------
      real Fabstol
      integer BORDER, BORDERZ
C     BORDER e BORDERZ: width and height of the zone for the calculation of
C     the dynamic field entering the normalization (also the
C     visualization zone for meshtv)
C     fmin: lowest frequency of modes to be calculated and plotted
C     Fabstol: the absolute error tolerance for the frequencies
      parameter (BORDER=1, BORDERZ=1, Fabstol=0.001*1.E09*2*3.14159)
C     .. Variables ..
C     Total number of active cells and of variables
      integer Nactive
C     ------------------------------------------------------------------
      integer LWORK
      data LWORK/1/
      real, allocatable :: Rabs(:), Rvec(:,:), TRvec(:,:)
      real, allocatable :: Eigvalr(:), Eigvali(:), Eigvalc(:)
      integer, allocatable :: Rindex(:)
      real, allocatable :: Eigvect(:,:)
      real, allocatable :: Theta0(:), Phi0(:)
      real, allocatable :: M(:,:), Dmdtheta(:,:), Dmdphi(:,:)
      real, allocatable :: D2mdtheta2(:,:), D2mdphi2(:,:) 
      real, allocatable :: D2mdthetadphi(:,:)
      real, allocatable :: Work(:)
      real, allocatable :: RWork(:)
      complex, allocatable :: CWork(:)
      complex, allocatable :: CB(:,:), CA(:,:)
      complex, allocatable :: CEigvect(:,:)
      integer, allocatable :: IWork(:), IFail(:)
C     ------------------------------------------------------------------
C     Input variables
C     Sel=0 workspace query,
C        =1 full calculation with approximate normalization
C        =2 full calculation with exact normalization
C     Verbosity=0 Print error messages and results only
C              =1 Print warning messages
C              =2 Print information messages
C              =3 Debug: print everything
C     The external field has intensity H0 and is directed at an angle
C     H0angle with respect to the x direction (enters the Zeeman
C     term of the derivatives and the BLS cross section).
C     The exchange interaction is calculated taking into account:
C     Fmax is the highest frequency of modes to be calculated and plotted
      logical unformatted
      integer Array, Sel, Verbosity, Cactive
      integer(4) xcount,ycount,zcount
      integer NumMat
      real(4) Cellsize, ZCellsize
      real H0, H0angle, H0polar, fmin, fmax, KBloch(2)
      character*80 infilename, Filename
      real, allocatable :: mm(:,:), Gamma(:)
      real, allocatable :: Exchange(:,:), Ms(:)
      real, allocatable :: K2p(:), Anisang(:,:)
      character*12, allocatable :: MatName(:)
      integer, allocatable :: Active(:), Iactive(:)
      integer, allocatable :: MatMap(:)
C     ------------------------------------------------------------------
C     Misc variables
      real Partwidth, Partheight
      real, allocatable :: anisvv(:,:)
      integer jj, n, l, Isterr, ierr
      integer Info, ii, k, i, j, p, infilelength
      integer ROW, COLUMN, LAYER
      integer Ccount, mati, matj
      integer jc, Goodmodes
      integer ix, iy, iz
      integer ir, sx, sy
      real Partthickness
      real PI
      real Time, second
      real Xcenter, Ycenter
      real Temp(13), Energy
      real DISTANCE
c     real DISTANCE, PDISTANCE
c     real Tdemag
      complex Eder(4,4)
      complex CTemp(13)
      real, allocatable :: demag(:,:,:,:,:)
      complex, allocatable :: edemag(:,:,:,:,:)
      real Efield, Ezee, Eexc, Edip, Eani
      real delta
      real Td
      complex Edelta
      integer XDISTANCE, YDISTANCE, ZDISTANCE
      complex Phase, Cemp
      INTERFACE 
        function tdemag(ex,ey,ez,erow,ecolumn,cellsize,zcellsize)
        !DEC$ ATTRIBUTES C, DECORATE, REFERENCE :: tdemag
        !DEC$ ATTRIBUTES REFERENCE :: ex, ey, ez, erow, ecolumn
        !DEC$ ATTRIBUTES REFERENCE :: cellsize, zcellsize
        real(4) tdemag, cellsize, zcellsize
        integer(4) ex, ey, ez, erow, ecolumn
        END FUNCTION
      END INTERFACE 
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
C     This is the x-component of the vector pointing from i to j in cell-units:
      XDISTANCE(i,j)=COLUMN(j)-COLUMN(i)
C     This is the y-component of the vector pointing from i to j in cell-units:
      YDISTANCE(i,j)=ROW(j)-ROW(i)
C     This is the z-component of the vector pointing from i to j in cell-units:
      ZDISTANCE(i,j)=LAYER(j)-LAYER(i)
C     This is the surface distance between cells i and j:
c     PDISTANCE(i,j)=SQRT(XDISTANCE(i,j)**2+YDISTANCE(i,j)**2+ZERO)
C     This is the distance between cells i and j:
      DISTANCE(i,j)=SQRT(XDISTANCE(i,j)**2+YDISTANCE(i,j)**2+
     2 ZDISTANCE(i,j)**2+ZERO)
C     ------------------------------------------------------------------
C     Note: our reference frame coincides with that of oommf.
C     ------------------------------------------------------------------

C     .. Executable Statements ..
      
C     First timemark
      Time = second()

      PI=4.*ATAN(ONE)
C     ------------------------------------------------------------------
C     Get arguments
      CALL GET_COMMAND_ARGUMENT (1,infilename,infilelength,ierr)
      if (ierr.ne.0) then
        print *, 'Wrong progam invocation. Syntax: '//
     2  'micro3d input-file-name'
        stop
      endif
C     ------------------------------------------------------------------
C     Data input
      call DATAIN(infilename,Sel,Verbosity,Cellsize,ZCellsize,
     2 Xcount,Ycount,Zcount,NumMat,MatName,Ms,Gamma,Exchange,H0,H0angle,
     3 H0polar,K2p,Anisang,mm,MatMap,Cactive,Active,Iactive,fmin,fmax,
     4 Filename,Array,KBloch,unformatted)
C     ------------------------------------------------------------------
      if (fmin.lt.1.E06) fmin=1.E06
c     fmin: 1.E06 = 1.E-03/2/PI GHz
      Ccount=Xcount*Ycount*Zcount
      print *, 'Number of active cells: ', Cactive
      allocate (anisvv(NumMat,3))
      Partwidth=Cellsize*Xcount
      Partheight=Cellsize*Ycount
      Partthickness=ZCellsize*Zcount
      Nactive=2*Cactive
C     ------------------------------------------------------------------
      if (Sel.GE.2) then
C       Open the derivatives scratch file
        open(UNIT=2,STATUS='SCRATCH',
     2  ERR=25,IOSTAT=Isterr,FORM='UNFORMATTED')
      endif
C     Dimension of the problem
      if(Verbosity.ge.2) then
        print *, 'Micro3d: array of ',
     2  Xcount,' x ',Ycount,' x ', Zcount,' cells'
        print *, 'Dimensions (cm):', Partwidth,Partheight,Partthickness
      endif
C     ------------------------------------------------------------------
      allocate (CB(Nactive,Nactive))
      allocate (CA(Nactive,Nactive))
      allocate (Eigvalc(Nactive))
      allocate (Theta0(Cactive))
      allocate (Phi0(Cactive))
      allocate (M(Cactive,3))
      allocate (Dmdtheta(Cactive,3))
      allocate (Dmdphi(Cactive,3))
      allocate (D2mdtheta2(Cactive,3))
      allocate (D2mdphi2(Cactive,3))
      allocate (D2mdthetadphi(Cactive,3))
      allocate (demag(-Xcount-BORDER+1:Xcount+BORDER-1,-Ycount-BORDER+1:
     2 Ycount+BORDER-1,-Zcount-BORDERZ+1:Zcount+BORDERZ-1,3,3))
      allocate (edemag(-Xcount-BORDER+1:Xcount+BORDER-1,
     2 -Ycount-BORDER+1:Ycount+BORDER-1,
     2 -Zcount-BORDERZ+1:Zcount+BORDERZ-1,3,3))
C     ------------------------------------------------------------------
      Xcenter=Xcount/2.+0.5
      Ycenter=Ycount/2.+0.5
C     ------------------------------------------------------------------
      Temp(1)=PI
      do k=1,Cactive  
C       Convert the oommf magnetizations to microdyn angles:
        if (mm(k,1).eq.ZERO.AND.mm(k,2).gt.ZERO) Theta0(k)=PI/2.
        if (mm(k,1).eq.ZERO.AND.mm(k,2).lt.ZERO) Theta0(k)=(3./2.)*PI
        if (mm(k,1).gt.ZERO) Theta0(k)=ATAN(mm(k,2)/mm(k,1))
        if (mm(k,1).lt.ZERO) Theta0(k)=ATAN(mm(k,2)/mm(k,1))+PI
        Phi0(k)=ACOS(mm(k,3))
        if(Phi0(k).lt.Temp(1)) Temp(1)=Phi0(k)
      enddo
      if (Temp(1)*180./PI.lt.5) then
C       At least one of the polar angles is less than 5 degrees
        print *, 
     2  'Micro3d: caution, magnetization almost in-plane somewhere'
        print *, 'Micro3d: smallest value among the polar angles:', 
     2  Temp(1)*180./PI, ' degree.'
      endif
C     ------------------------------------------------------------------
C     Generate the first smallest R (in number of MAXR)
      if (Array.eq.0) MAXR=1
      allocate (Rabs(MAXR*MAXR+2*MAXR+1))
      allocate (Rindex(MAXR*MAXR+2*MAXR+1))
      allocate (Rvec(MAXR+1,2))
      allocate (TRvec(MAXR*MAXR+2*MAXR+1,2))
      Rindex(1)=1
      Rvec(1,1)=0
      Rvec(1,2)=0
      if (Array.eq.1) then
        do i=1, MAXR/2
          Rvec(2*i,1)=i
          Rvec(2*i,2)=0
          Rvec(2*i+1,1)=-i
          Rvec(2*i+1,2)=0
        enddo
      else if (Array.eq.2) then
        do j=1, MAXR/2
          Rvec(2*j,1)=0
          Rvec(2*j,2)=j
          Rvec(2*j+1,1)=0
          Rvec(2*j+1,2)=-j
        enddo
      else if (Array.eq.3) then
C       Very brute force algorithm...
        do j=-MAXR/2, MAXR/2
          do i=-MAXR/2, MAXR/2
            ii=(MAXR/2+j)*(MAXR+1)+i+MAXR/2+1
            Rabs(ii)=(i*Partwidth)**2+(j*Partheight)**2
            TRvec(ii,1)=i
            TRvec(ii,2)=j
            Rindex(ii)=ii
          enddo
        enddo
        Call SORT(ii,MAXR*MAXR+2*MAXR+1,Rabs,Rindex)
        do i=1, MAXR
          Rvec(i,1)=TRvec(Rindex(i),1)
          Rvec(i,2)=TRvec(Rindex(i),2)
        enddo
      endif
      deallocate (Rabs)
      deallocate (Rindex)
      deallocate (TRvec)
C     Results: Rvec: components of the first MAXR smallest R
C              sorted for ascending R
c     do i=1, MAXR
c       print *, i, Rvec(i,1), Rvec(i,2)
c     enddo
C     ------------------------------------------------------------------
C     Assign the magnetization vectors calculated at equilbrium
      do i=1,Cactive
        M(i,1)=SIN(Phi0(i))*COS(Theta0(i))
        M(i,2)=SIN(Phi0(i))*SIN(Theta0(i))
        M(i,3)=COS(Phi0(i))
      enddo
C     ------------------------------------------------------------------
C     Assign the magnetization derivatives, calculated at equilbrium
      do i=1,Cactive
        Dmdtheta(i,1)=-SIN(Phi0(i))*SIN(Theta0(i))
        Dmdtheta(i,2)=SIN(Phi0(i))*COS(Theta0(i))
        Dmdtheta(i,3)=ZERO
        Dmdphi(i,1)=COS(Phi0(i))*COS(Theta0(i))
        Dmdphi(i,2)=COS(Phi0(i))*SIN(Theta0(i))
        Dmdphi(i,3)=-SIN(Phi0(i))
        D2mdtheta2(i,1)=-SIN(Phi0(i))*COS(Theta0(i))
        D2mdtheta2(i,2)=-SIN(Phi0(i))*SIN(Theta0(i))
        D2mdtheta2(i,3)=ZERO
        D2mdphi2(i,1)=-SIN(Phi0(i))*COS(Theta0(i))
        D2mdphi2(i,2)=-SIN(Phi0(i))*SIN(Theta0(i))
        D2mdphi2(i,3)=-COS(Phi0(i))
        D2mdthetadphi(i,1)=-COS(Phi0(i))*SIN(Theta0(i))
        D2mdthetadphi(i,2)=COS(Phi0(i))*COS(Theta0(i))
        D2mdthetadphi(i,3)=ZERO
      enddo
C     ------------------------------------------------------------------
C     Anisotropy easy axis cartesian coordinates:
      do i=1,NumMat
        anisvv(i,1)=SIN(Anisang(i,2))*COS(Anisang(i,1))
        anisvv(i,2)=SIN(Anisang(i,2))*SIN(Anisang(i,1))
        anisvv(i,3)=COS(Anisang(i,2))       
      enddo
C     ------------------------------------------------------------------
C     Allocate the work array (at least temporarily)
      allocate (CWork(LWORK))
      allocate (RWork(7*Nactive))
      allocate (IWork(5*Nactive))
      allocate (IFail(Nactive))
C     ------------------------------------------------------------------
C     Find the optimal LWORK
      allocate (Eigvalr(Nactive))
      allocate (CEigvect(Nactive,Nactive))
      call CHEGVX(1,'V','V','U',Nactive,CA,Nactive,CB,Nactive,
     2          1.E-05/fmax,1.E-05/fmin,0,0,Temp(1),
     3          Jc,Eigvalr,CEigvect,Nactive,CWork,-1,
     4          RWork,IWork,IFail,Info)
C     The following is for finding all the eigenvalues/eigenvectors
c     call CHEGV(1,'V','U',Nactive,CA,Nactive,CB,Nactive,
c    2          Eigvalr,CWork,-1,
c    3          RWork,Info)
c     info=0
      deallocate (Eigvalr)
      deallocate (CEigvect)
      if (Info.NE.0) then
        print *, 'Workspace query failed.'
        print *,'Cannot determine LWORK (and allocate the work array)'
        print *, 'Stopping...'
        stop
      else
        print *, 'Optimal LWORK:', REAL( CWork(1))
        if (Sel.eq.0) stop
        LWORK=INT(CWork(1))
        deallocate (CWork)
        allocate (CWork(LWORK))
      endif
C     ------------------------------------------------------------------
C     The (modified) demagnetization tensor is always real and symmetric 
C     w.r. to the inversion of the cell index, while the terms with 
C     exponentials (edemag) change by complex conjugation with inversion.
      ii=-Xcount-BORDER+1
C     Evaluate (once) the demagnetizing tensor
      do i=ii,Xcount+BORDER-1
        do j=-Ycount-BORDER+1,Ycount+BORDER-1
          do p=-Zcount-BORDERZ+1,Zcount+BORDERZ-1
            do k=1,3
              do l=1,3
C               the R=0 term must be always set
                demag(i,j,p,k,l)=Tdemag(int(i,4),int(j,4),int(p,4),
     2          int(k,4),int(l,4),Cellsize,ZCellsize)
                edemag(i,j,p,k,l)=demag(i,j,p,k,l)
C               Now, if MAXR>3, add other R's, in blocks of four
                do n=1, MAXR/4
                  delta=ZERO
                  edelta=CZERO
                  do ir=-2+4*n, 1+4*n
C                   Shift along x (number of cells):
                    sx=NINT(Partwidth/Cellsize)*Rvec(ir,1)
C                   Shift along y (number of cells):
                    sy=NINT(Partheight/Cellsize)*Rvec(ir,2)
c                   print *, 'preliminaries', i, j, p, k, l, ir
c                   print *, 'shift', sx, sy
c                   print *, 'Rvec', Rvec(ir,1), Rvec(ir,2)
                    Td=Tdemag(int(sx+i,4),int(sy+j,4),int(p,4),int(k,4),
     2              int(l,4),Cellsize,ZCellsize)
                    delta=delta+Td
                    edelta=edelta+EXP(2*PI*CI*(Kbloch(1)*Rvec(ir,1)+
     2              Kbloch(2)*Rvec(ir,2)))*Td
                  enddo
                  demag(i,j,p,k,l)=demag(i,j,p,k,l)+delta
                  edemag(i,j,p,k,l)=edemag(i,j,p,k,l)+edelta
C                 Is the required convergence obtained?
c                 print *, 'convergence', delta, demag(i,j,p,k,l)
                  if (ABS(delta).le.ABS(REPS*demag(i,j,p,k,l)).or.
     2            ABS(delta).le.EPS) goto 32
                enddo
                if (Array.ne.0) then
                  print *, 'micro3d: warning, problem at', i,j,p,k,l,n
                  print *, 'Too many steps for evaluating the '//
     2            'interdot coupling'
                  print *, demag(i,j,p,k,l)
                endif
  32          enddo
c             print *, 'max n:', n
              do l=1,k-1
                demag(i,j,p,k,l)=demag(i,j,p,l,k)
                edemag(i,j,p,k,l)=edemag(i,j,p,l,k)
  99          enddo
            enddo
          enddo
        enddo
      enddo
C     Here we exploit the symmetries of the demagnetization tensor 
C     mentioned above
      do i=-Xcount-BORDER+1-ii,-1
        do j=-Ycount-BORDER+1,Ycount+BORDER-1
          do p=-Zcount-BORDERZ+1,Zcount+BORDERZ-1
            do k=1,3
              do l=1,3
                demag(i,j,p,k,l)=demag(-i,-j,-p,k,l)
                edemag(i,j,p,k,l)=CONJG(edemag(-i,-j,-p,k,l))
              enddo
            enddo
          enddo
        enddo
      enddo
C     ------------------------------------------------------------------
      if(Verbosity.ge.2) 
     2 print *, 'Initialization time:', second()-Time, ' s'
C     Second timemark
      Time = second()
C     ------------------------------------------------------------------
C     The energy derivatives will be saved in the scratch file,
C     with the following scheme:
C     Fixed a submatrix corresponding to the couple of cells i,j (indexes
C     not shown below), the theta/phi derivatives are labeled as:
C     Eder(k,n)= derivative of E w.r. to theta/phi, with the following
C     ruleset:
C     k=1 phi theta
C     k=2 phi phi
C     k=3 theta theta
C     k=4 theta phi
C     The values of n correspond to:
C     n=1 Zeeman
C     n=2 Exchange
C     n=3 Dipolar
C     n=4 Anisotropy
C     The values of Eder are saved in the file as a function of cell
C     indexes i, j, with the following scheme:
C     repeat for i=1, Nactive
C       1. Diagonal term of the block j=i
C       2. Super-diagonal terms, cycle over j
C     end repeat i
C     so that the file is organized as follows:
C     (1,1), (1,2), (1,3), (1,4) ... (1,Nactive),
C     (2,2), (2,3), (2,4), ... (2,Nactive),
C     (3,3), (3,4), ... (3,Nactive),
C     ...
C     (Nactive,Nactive).
C     ------------------------------------------------------------------
C     Assign The Matrix
      do i=1,Nactive
        do j=1,Nactive
          CB(i,j)=CZERO
        enddo
      enddo
      do i=1,Cactive
C       i: block row index (and diagonal element block column index)
        ii=Iactive(i)
        mati=MatMap(i)
C       Diagonal terms:
        do j=1,6 
          Temp(j)=ZERO
        enddo
        do j=7,9
          CTemp(j)=CZERO
        enddo
        do j=1,Cactive
          jj=Iactive(j)
          matj=MatMap(j)
C         We need these sums for the dipolar term (including self interaction)
          do l=4,6
            Temp(l)=ZERO
            do n=1,3
              Temp(l)=Temp(l)+
     2	      demag(XDISTANCE(ii,jj),YDISTANCE(ii,jj),ZDISTANCE(ii,jj),
     3        (l-3),n)*D2mdtheta2(i,n)*Ms(matj)
            enddo
            CTemp(7)=CTemp(7)+M(j,(l-3))*Temp(l) 
          enddo
          do l=4,6
            Temp(l)=ZERO
            do n=1,3
              Temp(l)=Temp(l)+
     2        demag(XDISTANCE(ii,jj),YDISTANCE(ii,jj),ZDISTANCE(ii,jj),
     3        (l-3),n)*D2mdphi2(i,n)*Ms(matj)
            enddo
            CTemp(8)=CTemp(8)+M(j,(l-3))*Temp(l) 
          enddo
          do l=4,6
            Temp(l)=ZERO
            do n=1,3
              Temp(l)=Temp(l)+
     2        demag(XDISTANCE(ii,jj),YDISTANCE(ii,jj),ZDISTANCE(ii,jj),
     3        (l-3),n)*D2mdthetadphi(i,n)*Ms(matj)
            enddo
            CTemp(9)=CTemp(9)+M(j,(l-3))*Temp(l) 
          enddo
C         No self-interaction (exchange)
C         and first nearest neighbours
          if(DISTANCE(ii,jj).eq.1) then
C           We need these sums for the exchange term
            if(ZDISTANCE(ii,jj).eq.0) then
C             Couple in the same plane:
              do n=1,3
                Temp(1)=Temp(1)+D2mdtheta2(i,n)*M(j,n)/Cellsize**2
                Temp(2)=Temp(2)+D2mdphi2(i,n)*M(j,n)/Cellsize**2
                Temp(3)=Temp(3)+D2mdthetadphi(i,n)*M(j,n)/Cellsize**2
              enddo
            else
C             Couple in the same column:
              do n=1,3
                Temp(1)=Temp(1)+D2mdtheta2(i,n)*M(j,n)/ZCellsize**2
                Temp(2)=Temp(2)+D2mdphi2(i,n)*M(j,n)/ZCellsize**2
                Temp(3)=Temp(3)+D2mdthetadphi(i,n)*M(j,n)/ZCellsize**2
              enddo
            endif
          endif
          if (LAYER(ii).eq.LAYER(jj).and.
C         Same layer, different supercells
     2    (((Array.eq.1.or.Array.eq.3).and.COLUMN(ii).eq.1.and.
     3    COLUMN(jj).eq.Xcount.and.ROW(ii).eq.ROW(jj)).or.
C         We are dealing with an array periodic along x, i-cell is on the
C         first column, j-cell is on the last column, on the same row.
     4    ((Array.eq.1.or.Array.eq.3).and.COLUMN(ii).eq.Xcount.
     5    and.COLUMN(jj).eq.1.and.ROW(ii).eq.ROW(jj)).or.
C         We are dealing with an array periodic along x, i-cell is on the
C         last column, j-cell is on the first column, on the same row.
     6    ((Array.eq.2.or.Array.eq.3).and.ROW(ii).eq.1.and.
     7    ROW(jj).eq.Ycount.and.COLUMN(ii).eq.COLUMN(jj)).or.
C         We are dealing with an array periodic along y, i-cell is on the
C         first row, j-cell is on the last row, on the same column.
     8    ((Array.eq.2.or.Array.eq.3).and.ROW(ii).eq.Ycount.and.
     9    ROW(jj).eq.1.and.COLUMN(ii).eq.COLUMN(jj)))) then
C         We are dealing with an array periodic along y, i-cell is on the
C         last row, j-cell is on the first row, on the same column.
            do n=1,3
              Temp(1)=Temp(1)+D2mdtheta2(i,n)*M(j,n)/Cellsize**2
              Temp(2)=Temp(2)+D2mdphi2(i,n)*M(j,n)/Cellsize**2
              Temp(3)=Temp(3)+D2mdthetadphi(i,n)*M(j,n)/Cellsize**2
            enddo
          endif
  20    enddo
C       ----------------------------------------------------------------
C       phi i         theta j     terms with i=j
        do l=4,6
          CTemp(l)=CZERO
          do n=1,3
            CTemp(l)=CTemp(l)+edemag(0,0,0,(l-3),n)
     3      *Dmdtheta(i,n)
          enddo
          CTemp(9)=CTemp(9)+Dmdphi(i,(l-3))*CTemp(l)*Ms(mati)
        enddo
C       anisotropy ----  phi theta ----
        do n=10,13
          Temp(n)=ZERO
        enddo
        do n=1,3
          Temp(10)=Temp(10)+Dmdphi(i,n)*anisvv(mati,n)
          Temp(11)=Temp(11)+Dmdtheta(i,n)*anisvv(mati,n)
          Temp(12)=Temp(12)+M(i,n)*anisvv(mati,n)
          Temp(13)=Temp(13)+D2mdthetadphi(i,n)*anisvv(mati,n)
        enddo
C       -------------------------------------------------
        Eder(1,1)=-Ms(mati)*H0*(sin(H0polar)*cos(H0angle)*
     2  D2mdthetadphi(i,1)+sin(H0polar)*sin(H0angle)*D2mdthetadphi(i,2)
     3  +cos(H0polar)*D2mdthetadphi(i,3))
        Eder(1,2)=-2*Exchange(mati,mati)*Temp(3)
        Eder(1,3)=Ms(mati)*CTemp(9)
        Eder(1,4)=-2*K2p(mati)*(Temp(10)*Temp(11)+Temp(12)*Temp(13))
        CB(2*i,2*i-1)=Eder(1,1)+Eder(1,2)+Eder(1,3)+Eder(1,4)
C       phi i         phi j       terms with i=j
        do l=4,6
          CTemp(l)=CZERO
          do n=1,3
            CTemp(l)=CTemp(l)+edemag(0,0,0,(l-3),n)
     3      *Dmdphi(i,n)
          enddo
          CTemp(8)=CTemp(8)+ Dmdphi(i,(l-3))*CTemp(l)*Ms(mati)
        enddo
C       anisotropy ----  phi phi ----------------------
        do n=10,12
          Temp(n)=ZERO
        enddo 
        do n=1,3
          Temp(10)=Temp(10)+Dmdphi(i,n)*anisvv(mati,n)
          Temp(11)=Temp(11)+M(i,n)*anisvv(mati,n)
          Temp(12)=Temp(12)+D2mdphi2(i,n)*anisvv(mati,n)
        enddo
C       -----------------------------------------------------------
        Eder(2,1)=-Ms(mati)*H0*(D2mdphi2(i,1)*sin(H0polar)*cos(H0angle)+
     2  D2mdphi2(i,2)*sin(H0polar)*sin(H0angle) 
     3  +D2mdphi2(i,3)*cos(H0polar))
        Eder(2,2)=-2*Exchange(mati,mati)*Temp(2)
        Eder(2,3)=Ms(mati)*CTemp(8)
        Eder(2,4)=-2*K2p(mati)*(Temp(10)**2+Temp(11)*Temp(12))
        CB(2*i,2*i)=Eder(2,1)+Eder(2,2)+Eder(2,3)+Eder(2,4)
C       theta i       theta j      terms with i=j
        do l=4,6
          CTemp(l)=CZERO
          do n=1,3
            CTemp(l)=CTemp(l)+edemag(0,0,0,(l-3),n)
     3      *Dmdtheta(i,n)
          enddo
          CTemp(7)=CTemp(7)+Dmdtheta(i,(l-3))*CTemp(l)*Ms(mati)
        enddo
C       anisotropy ---- theta theta ------------------
        do n=10,12
          Temp(n)=ZERO
        enddo
        do n=1,3
          Temp(10)=Temp(10)+Dmdtheta(i,n)*anisvv(mati,n)
          Temp(11)=Temp(11)+M(i,n)*anisvv(mati,n)
          Temp(12)=Temp(12)+D2mdtheta2(i,n)*anisvv(mati,n)
        enddo
C       ----------------------------------------------------------------
        Eder(3,1)=-Ms(mati)*H0*(D2mdtheta2(i,1)*sin(H0polar)*
     2  cos(H0angle)+D2mdtheta2(i,2)*sin(H0polar)*sin(H0angle) 
     2  +D2mdtheta2(i,3)*cos(H0polar))
        Eder(3,2)=-2*Exchange(mati,mati)*Temp(1)
        Eder(3,3)=Ms(mati)*CTemp(7)
        Eder(3,4)=-2*K2p(mati)*(Temp(10)**2+Temp(11)*Temp(12))
        CB(2*i-1,2*i-1)=Eder(3,1)+Eder(3,2)+Eder(3,3)+Eder(3,4)
C       theta i       phi j        terms with i=j
        Eder(4,1)=Eder(1,1)
        Eder(4,2)=Eder(1,2)
        Eder(4,3)=Eder(1,3)
        Eder(4,4)=Eder(1,4)
C       When i=j (diagonal terms) the partial derivatives of E are real
        CB(2*i-1,2*i)= CB(2*i,2*i-1) 
        if (Sel.GE.2) then
C         Write the Eder diagonal element i:
          write(2) Eder
        endif
C       ----------------------------------------------------------------
C       Off-diagonal terms
C       Super-diagonal elements:
        do j=i+1,Cactive
C         j: block column index
          jj=Iactive(j)
          matj=MatMap(j)
          Eder(1,1)=CZERO
          Eder(2,1)=CZERO
          Eder(3,1)=CZERO
          Eder(4,1)=CZERO
          Eder(1,2)=CZERO
          Eder(2,2)=CZERO
          Eder(3,2)=CZERO
          Eder(4,2)=CZERO
C         first nearest neighbours
          if(DISTANCE(ii,jj).eq.1) then
            if(ZDISTANCE(ii,jj).eq.0) then
C             Couple in the same plane:
              do n=1,3
C               phi i         theta j
                Eder(1,2)=Eder(1,2)-2*Exchange(mati,matj)/Cellsize**2*
     2          Dmdphi(i,n)*Dmdtheta(j,n)
C               phi i         phi j
                Eder(2,2)=Eder(2,2)-2*Exchange(mati,matj)/
     2          Cellsize**2*Dmdphi(i,n)*Dmdphi(j,n)
C               theta i       theta j
                Eder(3,2)=Eder(3,2)-2*Exchange(mati,matj)/
     2          Cellsize**2*Dmdtheta(i,n)*Dmdtheta(j,n)
C               theta i       phi j
                Eder(4,2)=Eder(4,2)-2*Exchange(mati,matj)/
     2          Cellsize**2*Dmdtheta(i,n)*Dmdphi(j,n)
              enddo
            else
C             Couple in the same column:
              do n=1,3
C               phi i         theta j
                Eder(1,2)=Eder(1,2)-2*Exchange(mati,matj)/
     2          ZCellsize**2*Dmdphi(i,n)*Dmdtheta(j,n)
C               phi i         phi j
                Eder(2,2)=Eder(2,2)-2*Exchange(mati,matj)/
     2          ZCellsize**2*Dmdphi(i,n)*Dmdphi(j,n)
C               theta i       theta j
                Eder(3,2)=Eder(3,2)-2*Exchange(mati,matj)/
     2          ZCellsize**2*Dmdtheta(i,n)*Dmdtheta(j,n)
C               theta i       phi j
                Eder(4,2)=Eder(4,2)-2*Exchange(mati,matj)/
     2          ZCellsize**2*Dmdtheta(i,n)*Dmdphi(j,n)
              enddo
            endif
          endif
          if ((Array.eq.1.or.Array.eq.3).and.COLUMN(ii).eq.1.and.
     2    COLUMN(jj).eq.Xcount.and.ROW(ii).eq.ROW(jj).and.
     3    LAYER(ii).eq.LAYER(jj)) then
C           We are dealing with an array periodic along x, i-cell is on the
C           first column, j-cell is on the last column, on the same row
C           and layer.
            do n=1,3
              Eder(1,2)=Eder(1,2)-2*Exchange(mati,matj)/Cellsize**2*
     2        Dmdphi(i,n)*Dmdtheta(j,n)
     3        *EXP(-2*PI*CI*Kbloch(1))
              Eder(2,2)=Eder(2,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdphi(i,n)*Dmdphi(j,n)
     3        *EXP(-2*PI*CI*Kbloch(1))
              Eder(3,2)=Eder(3,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdtheta(j,n)
     3        *EXP(-2*PI*CI*Kbloch(1))
              Eder(4,2)=Eder(4,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdphi(j,n)
     3        *EXP(-2*PI*CI*Kbloch(1))
            enddo
          endif
          if ((Array.eq.1.or.Array.eq.3).and.COLUMN(ii).eq.Xcount.
     2    and.COLUMN(jj).eq.1.and.ROW(ii).eq.ROW(jj).and.
     3    LAYER(ii).eq.LAYER(jj)) then
C           We are dealing with an array periodic along x, i-cell is on the
C           last column, j-cell is on the first column, on the same row
C           and layer.
            do n=1,3
              Eder(1,2)=Eder(1,2)-2*Exchange(mati,matj)/Cellsize**2*
     2        Dmdphi(i,n)*Dmdtheta(j,n)
     3        *EXP(2*PI*CI*Kbloch(1))
              Eder(2,2)=Eder(2,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdphi(i,n)*Dmdphi(j,n)
     3        *EXP(2*PI*CI*Kbloch(1))
              Eder(3,2)=Eder(3,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdtheta(j,n)
     3        *EXP(2*PI*CI*Kbloch(1))
              Eder(4,2)=Eder(4,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdphi(j,n)
     3        *EXP(2*PI*CI*Kbloch(1))
            enddo
          endif
          if ((Array.eq.2.or.Array.eq.3).and.ROW(ii).eq.1.and.
     2    ROW(jj).eq.Ycount.and.COLUMN(ii).eq.COLUMN(jj).and.
     3    LAYER(ii).eq.LAYER(jj)) then
C           We are dealing with an array periodic along y, i-cell is on the
C           first row, j-cell is on the last row, on the same column
C           and layer.
            do n=1,3
              Eder(1,2)=Eder(1,2)-2*Exchange(mati,matj)/Cellsize**2*
     2        Dmdphi(i,n)*Dmdtheta(j,n)*
     3        EXP(-2*PI*CI*Kbloch(2))
              Eder(2,2)=Eder(2,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdphi(i,n)*Dmdphi(j,n)*
     3        EXP(-2*PI*CI*Kbloch(2))
              Eder(3,2)=Eder(3,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdtheta(j,n)*
     3        EXP(-2*PI*CI*Kbloch(2))
              Eder(4,2)=Eder(4,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdphi(j,n)*
     3        EXP(-2*PI*CI*Kbloch(2))
            enddo
          endif
          if ((Array.eq.2.or.Array.eq.3).and.ROW(ii).eq.Ycount.and.
     2    ROW(jj).eq.1.and.COLUMN(ii).eq.COLUMN(jj).and.
     3    LAYER(ii).eq.LAYER(jj)) then
C           We are dealing with an array periodic along y, i-cell is on the
C           last row, j-cell is on the first row, on the same column
C           and layer.
            do n=1,3
              Eder(1,2)=Eder(1,2)-2*Exchange(mati,matj)/Cellsize**2*
     2        Dmdphi(i,n)*Dmdtheta(j,n)*
     3        EXP(2*PI*CI*Kbloch(2))
              Eder(2,2)=Eder(2,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdphi(i,n)*Dmdphi(j,n)*
     3        EXP(2*PI*CI*Kbloch(2))
              Eder(3,2)=Eder(3,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdtheta(j,n)*
     3        EXP(2*PI*CI*Kbloch(2))
              Eder(4,2)=Eder(4,2)-2*Exchange(mati,matj)/
     2        Cellsize**2*Dmdtheta(i,n)*Dmdphi(j,n)*
     3        EXP(2*PI*CI*Kbloch(2))
            enddo
          endif
          CB(2*i,2*j-1)=Eder(1,2)
          CB(2*i,2*j)=Eder(2,2)
          CB(2*i-1,2*j-1)=Eder(3,2)
          CB(2*i-1,2*j)=Eder(4,2)
C         phi i         theta j
          Eder(1,3)=CZERO
          Eder(1,4)=CZERO
          do n=1,3
            CTemp(n)=CZERO
            do l=1,3    
              CTemp(n)=CTemp(n)+edemag(XDISTANCE(ii,jj),
     2        YDISTANCE(ii,jj),ZDISTANCE(ii,jj),n,l)*Dmdtheta(j,l)
            enddo
            Eder(1,3)=Eder(1,3)+Ms(mati)*Ms(matj)*Dmdphi(i,n)*CTemp(n)
          enddo
          CB(2*i,2*j-1)=CB(2*i,2*j-1)+Eder(1,3)
C         phi i         phi j
          Eder(2,3)=CZERO
          Eder(2,4)=CZERO
          do n=1,3
            CTemp(n)=CZERO
            do l=1,3    
              CTemp(n)=CTemp(n)+edemag(XDISTANCE(ii,jj),
     2        YDISTANCE(ii,jj),ZDISTANCE(ii,jj),n,l)*Dmdphi(j,l)
            enddo
            Eder(2,3)=Eder(2,3)+Ms(mati)*Ms(matj)*Dmdphi(i,n)*CTemp(n)
          enddo
          CB(2*i,2*j)=CB(2*i,2*j)+Eder(2,3)
C         theta i       theta j
          Eder(3,3)=CZERO
          Eder(3,4)=CZERO
          do n=1,3
            CTemp(n)=CZERO
            do l=1,3    
	      CTemp(n)=CTemp(n)+edemag(XDISTANCE(ii,jj),
     2        YDISTANCE(ii,jj),ZDISTANCE(ii,jj),n,l)*Dmdtheta(j,l)
            enddo
            Eder(3,3)=Eder(3,3)+Ms(mati)*Ms(matj)*Dmdtheta(i,n)*CTemp(n)
	  enddo
          CB(2*i-1,2*j-1)=CB(2*i-1,2*j-1)+Eder(3,3)
C	  theta i       phi j
          Eder(4,3)=CZERO
          Eder(4,4)=CZERO
          do n=1,3
            CTemp(n)=CZERO
            do l=1,3    
	      CTemp(n)=CTemp(n)+edemag(XDISTANCE(ii,jj),
     2        YDISTANCE(ii,jj),ZDISTANCE(ii,jj),n,l)*Dmdphi(j,l)
            enddo
            Eder(4,3)=Eder(4,3)+Ms(mati)*Ms(matj)*Dmdtheta(i,n)*CTemp(n)
	  enddo
          CB(2*i-1,2*j)=CB(2*i-1,2*j)+Eder(4,3)
          if (Sel.GE.2) then
C           Write the Eder element i,j:
            write(2) Eder
          endif
        enddo
C       Sub-diagonal elements not needed
      enddo
      if(Verbosity.ge.2) 
     2 print *, 'Time for assigning The Matrix:', second()-Time, ' s'
C     Third timemark
      Time = second()
C     ------------------------------------------------------------------
C     Define the Hermitian matrix A
      do i=1,Nactive
        mati=MatMap((i+1)/2)
        do j=i,Nactive
          if (j-i.EQ.1.and.2*(j/2).EQ.j) then
            CA(i,j)=CI*SIN(Phi0((i+1)/2))*Ms(mati)/(Gamma(mati)/1.E05)
          else
            CA(i,j)=CZERO
          endif
        enddo
      enddo
      allocate (Eigvalr(Nactive))
      allocate (CEigvect(Nactive,Nactive))
C     ------------------------------------------------------------------
C     Set the absolute error tolerance for the eigenvalues
      Temp(1)=1.E-05*Fabstol/fmax/fmax
C     ------------------------------------------------------------------
C     The following is for the highest accuracy
c     temp(1)=2*SLAMCH('S')
C     ------------------------------------------------------------------
C     The following is for finding all the eigenvalues/eigenvectors
      call CHEGV(1,'V','U',Nactive,CA,Nactive,CB,Nactive,
     2           Eigvalr,CWork,LWork,RWork,Info)
      Jc=Nactive
      do i=1, Nactive
        do j=1, Nactive
          CEigvect(i,j)=CA(i,j)
        enddo
      enddo
C     ------------------------------------------------------------------
  55  if (Info.NE.0) then
        if (Info.LT.0) then
          print *, 'The argument ',-info,' has an illegal value'
          stop
        else
          if (Info.GT.Nactive) then
            print *, 'Micro3d: the given static magnetization does not '
     2      //'correspond to an energy minimum!'
            print *
            print *, 'Micro3d: calculation of eigenfrequencies failed.'
            stop
          else
            print *, 'Micro3d: several eigenfrequencies failed to '
     2      //'converge. They are:'
            do i=1,Info
              print *, '  Mode', IFail(i)
              print *, '  Frequency', 
     2        1.E+05/Eigvalr(IFail(i))/2/PI/1.E09
            enddo
            print *, 'Micro3d: we continue using the latest '//
     2      ' approximation to the eigenvalue and eigenvector'
          endif
        endif
      endif
C     Goodmodes is the number of good modes
      Goodmodes=Jc 
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      if(Verbosity.ge.2) 
     2 print *, 'Time for solving the eigenproblem:',second()-Time,' s'
C     Fourth timemark
      Time = second()
C     ------------------------------------------------------------------
C     Print the results
      if (Verbosity.GE.2) print *, 'Micro3d: optimal LWORK:', 
     2 CWork(1), '; actual value:', LWORK
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      if (ALLOCATED(Work)) deallocate (Work)
      if (ALLOCATED(CWork)) deallocate (CWork)
      if (ALLOCATED(RWork)) deallocate (RWork)
      if (ALLOCATED(CB)) deallocate (CB)
      if (ALLOCATED(Eigvalc)) deallocate (Eigvalc)
      if (ALLOCATED(Eigvali)) deallocate (Eigvali)
      if (ALLOCATED(Eigvect)) deallocate (Eigvect)
C     In the following we will use the first eigenvalue (positive frequency)
C     ------------------------------------------------------------------
C     write the relevant data in the output file
      if (unformatted) then
        write(1) SIZE(Ms)
        write(1) Cellsize,ZCellsize,
     2  Xcount,Ycount,Zcount,MatName,Ms,Gamma,Exchange,H0,H0angle,
     3  H0polar,K2p,Anisang,Cactive,fmin,fmax,Xcenter,Ycenter,
     4  (Array.ne.0),KBloch
        write(1) Goodmodes
        do i=1, Ccount
          write(1) Active(i)
        enddo
        do i=1, Cactive
          write(1) Iactive(i), MatMap(i), Theta0(i), Phi0(i)
        enddo
      else
        write(1,*) SIZE(Ms)
        write(1,*) Cellsize,ZCellsize,
     2  Xcount,Ycount,Zcount,MatName,Ms,Gamma,Exchange,H0,H0angle,
     3  H0polar,K2p,Anisang,Cactive,fmin,fmax,Xcenter,Ycenter,
     4  (Array.ne.0),KBloch
        write(1,*) Goodmodes
        do i=1, Ccount
          write(1,*) Active(i)
        enddo
        do i=1, Cactive
          write(1,*) Iactive(i), MatMap(i), Theta0(i), Phi0(i)
        enddo
      endif

c       print *, Cellsize,ZCellsize,
c    2  Xcount,Ycount,Zcount,MatName,Ms,Gamma,Exchange,H0,H0angle,
c    3  H0polar,K2p,Anisang,Cactive,fmin,fmax,Xcenter,Ycenter,
c    4  (Array.ne.0),KBloch

      do j=1,Goodmodes
        if (1.E05/Eigvalr(j).ge.fmin.and.1.E05/Eigvalr(j).lt.fmax) 
     2  then
          if (unformatted) then
            write(1) 1.E05/Eigvalr(j)
          else
            write(1,*) 1.E05/Eigvalr(j)
          endif
          Temp(1)=ZERO
          Temp(2)=ZERO
          Phase=CZERO
          do i=1,Cactive
C           Temp(1) stores the square modulus of the dynamic magnetization
C           Phase stores the complex phase of dPhi and dTheta, obtained
C           by successive approximations.
            Temp(1)=Temp(1)+ABS(CEigvect(2*i-1,j))**2
     2      *SIN(Phi0(i))**2+ ABS(CEigvect(2*i,j))**2
            Cemp=CEigvect(2*i-1,j)* EXP(-CI*Phase)
            Temp(2)=ATAN(-REAL(Cemp)/IMAG(Cemp))
            if(Temp(2).lt.-PI/2) Temp(2)=Temp(2)+PI/2
            if(Temp(2).gt.PI/2) Temp(2)=Temp(2)-PI/2
            Phase=Phase+Temp(2)/(2*i-1)
            Cemp=CEigvect(2*i,j)*EXP(-CI*Phase)
            Temp(2)=ATAN(IMAG(Cemp)/REAL(Cemp))
            if(Temp(2).lt.-PI/2) Temp(2)=Temp(2)+PI/2
            if(Temp(2).gt.PI/2) Temp(2)=Temp(2)-PI/2
            Phase=Phase+Temp(2)/(2*i)
          enddo
          Temp(1)=Temp(1)/Cactive
          do i=1,Cactive
            Cemp=CEigvect(2*i-1,j)* EXP(-CI*Phase)/SQRT(Temp(1))
            CEigvect(2*i-1,j)=Cemp
            Cemp=CEigvect(2*i,j)* EXP(-CI*Phase)/SQRT(Temp(1))
            CEigvect(2*i,j)=Cemp
          enddo
C         --------------------------------------------------------------
C         True normalization (actually we calculate the energy of a
C         mode, letting the following code to use it)
          if (Sel.LE.1) then
C           skip the true normalization:
            Energy=ONE
            GOTO 60
          endif
          rewind(2)
C         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
C         Energy of the magnetic field
C         We calculate the energy of the magnetic field in the space
C         by summing over a volume extending BORDERxBORDERxBORDERZ cells
C         (linearly) beyond the magnetic particle:
          Efield=ZERO
          do ix=1-BORDER, Xcount+BORDER
            do iy=1-BORDER, Ycount+BORDER
              do iz=1-BORDERZ, Zcount+BORDERZ
C               ... and now let's sum the contributions of the active cells:
                do n=1,3
                  Temp(n)=ZERO
                  do i=1,Cactive
                    ii=Iactive(i)
                    mati=MatMap(i)
                    Temp(n)=Temp(n)+Ms(mati)*(
     2              demag(ix-COLUMN(ii),iy-ROW(ii),iz-LAYER(ii),n,1)*
     3              ( -SIN(Phi0(i))*SIN(Theta0(i))*
     4              CEigvect(2*i-1,j)+
     5              COS(Phi0(i))*COS(Theta0(i))*
     6              CEigvect(2*i,j) )+
     7              demag(ix-COLUMN(ii),iy-ROW(ii),iz-LAYER(ii),n,2)*
     8              ( SIN(Phi0(i))*COS(Theta0(i))*
     9              CEigvect(2*i-1,j) +
     *              COS(Phi0(i))*SIN(Theta0(i))*
     1              CEigvect(2*i,j) )-
     2              demag(ix-COLUMN(ii),iy-ROW(ii),iz-LAYER(ii),n,3)*
     3              SIN(Phi0(i))*CEigvect(2*i,j))
                  enddo
                enddo
                Efield=Efield+Temp(1)**2+Temp(2)**2+Temp(3)**2
              enddo
            enddo
          enddo
          Efield=Efield/8/PI
C         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
C         All the other energy contributions: Zeeman, exchange, dipolar,
C         anisotropy
          Ezee=ZERO
          Eexc=ZERO
          Edip=ZERO
          Eani=ZERO
          do i=1,Cactive
C           i: block row index (and diagonal element block index)
C           Eigvect(2*i,j:j+1)=d phi; Eigvect(2*i-1,j:j+1)=d theta
C           Read the Eder diagonal element i:
            read(2) Eder
            Temp(1)=REAL(CEigvect(2*i,j))
            Temp(2)=IMAG(CEigvect(2*i,j))
            Temp(3)=REAL(CEigvect(2*i-1,j))
            Temp(4)=IMAG(CEigvect(2*i-1,j))
            Ezee=Ezee+Eder(1,1)*(Temp(1)*Temp(3)+
     2      Temp(2)*Temp(4))+
     3      0.5*Eder(2,1)*(Temp(1)*Temp(1)+
     4      Temp(2)*Temp(2))+
     5      0.5*Eder(3,1)*(Temp(3)*Temp(3)+
     6      Temp(4)*Temp(4))
            Eexc=Eexc+Eder(1,2)*(Temp(1)*Temp(3)+
     2      Temp(2)*Temp(4))+
     3      0.5*Eder(2,2)*(Temp(1)*Temp(1)+
     4      Temp(2)*Temp(2))+
     5      0.5*Eder(3,2)*(Temp(3)*Temp(3)+
     6      Temp(4)*Temp(4))
            Edip=Edip+Eder(1,3)*(Temp(1)*Temp(3)+
     2      Temp(2)*Temp(4))+
     3      0.5*Eder(2,3)*(Temp(1)*Temp(1)+
     4      Temp(2)*Temp(2))+
     5      0.5*Eder(3,3)*(Temp(3)*Temp(3)+
     6      Temp(4)*Temp(4))
            Eani=Eani+Eder(1,4)*(Temp(1)*Temp(3)+
     2      Temp(2)*Temp(4))+
     3      0.5*Eder(2,4)*(Temp(1)*Temp(1)+
     4      Temp(2)*Temp(2))+
     5      0.5*Eder(3,4)*(Temp(3)*Temp(3)+
     6      Temp(4)*Temp(4))
            do n=i+1,Cactive
C             n: block column index
              Temp(5)=REAL(CEigvect(2*n,j))
              Temp(6)=IMAG(CEigvect(2*n,j))
              Temp(7)=REAL(CEigvect(2*n-1,j))
              Temp(8)=IMAG(CEigvect(2*n-1,j))
C             Read the Eder element i,n:
              read(2) Eder
              Ezee=Ezee+Eder(1,1)*(Temp(1)*Temp(7)+
     2        Temp(2)*Temp(8))+
     3        Eder(2,1)*(Temp(1)*Temp(5)+
     4        Temp(2)*Temp(6))+
     5        Eder(3,1)*(Temp(3)*Temp(7)+
     6        Temp(4)*Temp(8))+
     7        Eder(4,1)*(Temp(3)*Temp(5)+
     8        Temp(4)*Temp(6))
              Eexc=Eexc+Eder(1,2)*(Temp(1)*Temp(7)+
     2        Temp(2)*Temp(8))+
     3        Eder(2,2)*(Temp(1)*Temp(5)+
     4        Temp(2)*Temp(6))+
     5        Eder(3,2)*(Temp(3)*Temp(7)+
     6        Temp(4)*Temp(8))+
     7        Eder(4,2)*(Temp(3)*Temp(5)+
     8        Temp(4)*Temp(6))
              Edip=Edip+Eder(1,3)*(Temp(1)*Temp(7)+
     2        Temp(2)*Temp(8))+
     3        Eder(2,3)*(Temp(1)*Temp(5)+
     4        Temp(2)*Temp(6))+
     5        Eder(3,3)*(Temp(3)*Temp(7)+
     6        Temp(4)*Temp(8))+
     7        Eder(4,3)*(Temp(3)*Temp(5)+
     8        Temp(4)*Temp(6))
              Eani=Eani+Eder(1,4)*(Temp(1)*Temp(7)+
     2        Temp(2)*Temp(8))+
     3        Eder(2,4)*(Temp(1)*Temp(5)+
     4        Temp(2)*Temp(6))+
     5        Eder(3,4)*(Temp(3)*Temp(7)+
     6        Temp(4)*Temp(8))+
     7        Eder(4,4)*(Temp(3)*Temp(5)+
     8        Temp(4)*Temp(6))
  40        enddo
          enddo
c         Try to not use the dipolar term:
          Energy=Efield+Ezee+Eexc+Eani
c         Energy=Efield+Ezee+Eexc+Edip+Eani
          if(Verbosity.ge.2) then
            print *, 'Energy contributions, mode', 
     2      1.E05/Eigvalr(j)/2/PI/1.E09, ' GHz:'
            print *, 'Field:     ', Efield
            print *, 'Zeeman:    ', Ezee
            print *, 'Exchange:  ', Eexc
            print *, 'Dipolar:   ', Edip
            print *, 'Anisotropy:', Eani
          endif
  60      if (unformatted) then
            write(1) Energy
          else
            write(1,*) Energy
          endif
C         --------------------------------------------------------------
C         Write the eigenvectors
          if (unformatted) then
            do i=1,Nactive
              write(1) REAL(CEigvect(i,j)), IMAG(CEigvect(i,j))
            enddo
          else
            do i=1,Nactive
              write(1,*) REAL(CEigvect(i,j)), IMAG(CEigvect(i,j))
            enddo
          endif
        endif
      enddo
      close(1)
      close(2)
      if(Verbosity.ge.2)
     2 print *, 'Time for writing the results:',second()-Time,' s'
      stop
  25  print *, 'micro3d: error opening scratch file, status: ', Isterr
      stop
      end
