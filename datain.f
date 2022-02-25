      module input_parameters
      contains

      subroutine DATAIN(infilename,Sel,Verbosity,Cellsize,ZCellsize,
     2 Xcount,Ycount,Zcount,NumMat,MatName,Ms,Gamma,Exchange,H0,H0angle,
     3 H0polar,K2p,Anisang,mm,MatMap,Cactive,Active,Iactive,fmin,fmax,
     4 Filename,Array,KBloch,unformatted)
C.......................................................................
C
C     Read magnetic constants of the problem from infilename
C     and dot geometry and spin mesh from a .omf file
C.......................................................................
C
C     UNITS:
C     L --> cm
C     M --> g
C     T --> s
C
C     MAGN. FIELD --> Oe
C
C     Note:
C       The gyromagnetic ratio is Gamma=G*E/(2*M*C)
C       E/(2*M*C)=8.7939532E06 s sC (g cm)^-1   (cgs)
C.......................................................................
C
C     Input variables
C     Sel=0 workspace query,
C        =1 full calculation.
C     Verbosity=0 Print error messages and results only
C              =1 Print warning messages
C              =2 Print information messages
C              =3 Debug: print everything
C     Array=0 Single dot
C           1 1D periodic structure, repeated along x
C           2 1D periodic structure, repeated along y
C           3 2D periodic structure, repeated along x and y

      use string_utility

      implicit none

      logical unformatted
      integer Array, Sel, Verbosity, Cactive
      integer(4) xcount,ycount,zcount
      integer NumMat
      real(4) Cellsize, ZCellsize
      real H0, H0angle, H0polar, fmin, fmax, KBloch(2)
      character*80 infilename, Filename
      real, allocatable, intent(out) :: mm(:,:), Gamma(:)
      real, allocatable, intent(out) :: Exchange(:,:), Ms(:)
      real, allocatable, intent(out) :: K2p(:), Anisang(:,:)
      character*12, allocatable, intent(out) :: MatName(:)
      integer, allocatable, intent(out) :: Active(:), Iactive(:)
      integer, allocatable, intent(out) :: MatMap(:)
C     ------------------------------------------------------------------
C     Internal variables:
      logical Err
      integer Line, Isterr, i, ii, j, k, m
      integer(4) StartLine
C     MAXNMAT: Maximum number of different materials
      integer MAXNMAT
      parameter (MAXNMAT=10)
      real ZERO, PI, avtab(MAXNMAT), absvec
      real*4 MscaleFactor
      real TMs, TK2p, TAnisang1, TAnisang2, TExchange, DMS, TGamma
      character*80 OmfFilename, MatTableFN
      character*12 TMatName, TMatName2, forms
      logical, allocatable :: check(:), mcheck(:,:)
      real, allocatable :: B(:,:,:,:)
C     ------------------------------------------------------------------
      INTERFACE 
        subroutine inquiryomf(xcount,ycount,zcount,startline,cellsize,
     2  zcellsize,mscalefactor,omffilename)
        !DEC$ ATTRIBUTES C, DECORATE, REFERENCE :: inquiryomf 
        !DEC$ ATTRIBUTES REFERENCE :: omffilename
        CHARACTER*(*) omffilename 
        !DEC$ ATTRIBUTES REFERENCE :: xcount,ycount,zcount,startLine
        !DEC$ ATTRIBUTES REFERENCE :: cellsize,zcellsize,mscalefactor
        real(4) cellsize, zcellsize, mscalefactor
        integer(4) xcount,ycount,zcount,startline
        END SUBROUTINE 
      END INTERFACE 
C     ------------------------------------------------------------------
C     DMS: difference of Ms (A/m) for the two corresponding materials to be 
C          considered different
      parameter (ZERO=0.E0, DMS=1.E00)
      PI=4.*ATAN(1.E0)
      Err=.FALSE.
      open(UNIT=1,FILE=infilename,STATUS='OLD',ERR=22,IOSTAT=Isterr)
      Line=1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11) Filename, unformatted
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11) Sel, Verbosity
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11) H0, H0angle, H0polar
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11) OmfFilename, MatTableFN
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11) fmin, fmax
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11)
      Line=Line+1
      read(1,*,ERR=11) Array, KBloch
      close(1)
C     ------------------------------------------------------------------
C     Sanity tests
      if(Sel.lt.0.or.Sel.GT.2) then
        print *, 'Datain: ',
     2  '"Sel" out of range [0-1]: what do you want to do?'
        Err=.TRUE.
      endif
      if(Verbosity.lt.0.or.Verbosity.GT.3) then
        print *, 'Datain: ',
     2  '"Verbosity" out of range [0-3]: what do you want to print?'
        Err=.TRUE.
      endif
      if(H0.lt.0) then
        print *, 'Datain: the applied field is negative!'
        Err=.TRUE.
      endif
      if(H0angle.lt.-360.or.H0angle.gt.360) then
        print *, 'Datain: strange value for the external field  in-plane
     * angle...'
        Err=.TRUE.
      endif
      if(H0polar.lt.-360.or.H0polar.gt.360) then
        print *, 'Datain: strange value for the external field  polar an
     *gle...'
        Err=.TRUE.
      endif
      if (Array.lt.0.OR.Array.gt.3) then
        print *, 'Datain: Array is not valid: it should be either '
        print *, '0, 1, 2 or 3'
        Err=.TRUE.
      endif
C     ------------------------------------------------------------------
      if (Err) then
        print *
        print *, 'Stopping due to input error(s)...'
        stop
      endif
C     ------------------------------------------------------------------
C     Find the dimension of the problem, scale factor and starting line
C     in the OMF input file
      OmfFilename=TRIM(ADJUSTL(OmfFilename))//char(0)
c     print *, 'kind xcount', KIND(Xcount)
      call inquiryomf(xcount,ycount,zcount,startline,cellsize,
     2 zcellsize,mscalefactor,omffilename)
c     print *, 'xcount,ycount,zcount,startline,cellsize,
c    2 zcellsize,mscalefactor,omffilename'
c     print *, xcount,ycount,zcount,startline,cellsize,
c    2 zcellsize,mscalefactor,omffilename
C     Change units
      MscaleFactor=MscaleFactor/1.E03
C     Assign the magnetization
C     ------------------------------------------------------------------
      open (UNIT=3,FILE=OmfFilename,STATUS='OLD',ERR=13,
     2 IOSTAT=Isterr)
      Line=1
      do j=1,StartLine
        Line=Line+1
        read (3,*,ERR=12)
      enddo
      allocate (B(Xcount,Ycount,Zcount,3))
C     ----------------------------------------------------------------
C     Read static magnetization and count the active elements
      Cactive=0
      do k=1,Zcount
        do j=1,Ycount
          do i=1,Xcount
            Line=Line+1
            read (3,*,ERR=12) B(i,j,k,1),B(i,j,k,2),B(i,j,k,3)
c           print *, i, j, k, B(i,j,k,1),B(i,j,k,2),B(i,j,k,3)
            if (B(i,j,k,1)**2+B(i,j,k,2)**2+B(i,j,k,3)**2.ge.1.e-5)
     2      Cactive=Cactive+1
          enddo
        enddo
      enddo  
      close(3)
C     ----------------------------------------------------------------
      allocate (mm(Cactive,3))
      allocate (Active(Xcount*Ycount*Zcount))
      allocate (Iactive(Cactive))
      allocate (MatMap(Cactive))
C     ----------------------------------------------------------------
C     Identify the active elements, store their indices
C     and find out the material map
      ii=0
      NumMat=0
      do i=1,MAXNMAT
        avtab(i)=ZERO
      enddo
      do k=1,Zcount
        do j=1,Ycount
          do i=1,Xcount
            absvec=sqrt(B(i,j,k,1)**2+B(i,j,k,2)**2+B(i,j,k,3)**2)*
     2      Mscalefactor
            if (absvec.lt.0.01) then
C             Empy cell
              Active(i+Xcount*(j-1)+Xcount*Ycount*(k-1))=0     
            else
C             Active cell
              Active(i+Xcount*(j-1)+Xcount*Ycount*(k-1))=1     
              ii=ii+1
              Iactive(ii)=i+Xcount*(j-1)+Xcount*Ycount*(k-1)
C             Identify the constituent material
              do m=1,NumMat
                if (abs(absvec-avtab(m)).LE.DMS) then
C                 Old material found
                  MatMap(ii)=m
                  goto 10
                endif
              enddo
C             New material found
              NumMat=NumMat+1
              if (NumMat.gt.MAXNMAT) THEN
                print *, 'Datain: maximum number of materials '//
     2          'exceeded. Please change MAXNMAT and recompile.'
                stop
              endif
              avtab(NumMat)=absvec
              MatMap(ii)=NumMat
  10          continue
C             Assign the static normalized magnetization
              do m=1,3
                mm(ii,m)=B(MOD(MOD(Iactive(ii)-1,Xcount*Ycount),Xcount)
     2          +1,MOD(Iactive(ii)-1,Xcount*Ycount)/Xcount+1,
     3          (Iactive(ii)-1)/(Xcount*Ycount)+1,m)/
     4          (absvec/Mscalefactor)
              enddo
            endif
          enddo
        enddo
      enddo  
      deallocate (B)
      allocate (MatName(NumMat))
      allocate (Ms(NumMat))
      allocate (Gamma(NumMat))
      allocate (K2p(NumMat))
      allocate (Anisang(NumMat,2))
      allocate (check(NumMat))
      print *, 'Number of materials:', NumMat
      do j=1,NumMat
        check(j)=.false.
      enddo
      open (UNIT=3,FILE=MatTableFN,STATUS='OLD',ERR=14,IOSTAT=Isterr)
C     Read the material specifications
      read (3,*)
      do i=1, MAXNMAT
        read (3,*,ERR=30) TMatName, TMs, TGamma, TK2p, TAnisang1, 
     2  TAnisang2
        do j=1,NumMat
          if (abs(TMS-avtab(j)).LT.DMS) then
            if(check(j)) then
              print *, 'Datain: duplicated material in ', MatTableFN
              stop
            endif
            MatName(j)=STRLOWCASE(ADJUSTL(TMatName))
            Ms(j)=TMs
            Gamma(j)=TGamma
            K2p(j)=TK2p
            Anisang(j,1)=TAnisang1*PI/180.
            Anisang(j,2)=TAnisang2*PI/180.
            check(j)=.true.
          endif
        enddo
      enddo
C     Too many materials in the first part of the table
      print *, 'Datain: too many materials in ', MatTableFN
      stop
C     End of first part of the table reached
  30  continue
      print *, 'MATERIALS:'
      do j=1,NumMat
        if (.not.check(j)) then
          print *, 'Datain: missing required material in ', MatTableFN
          stop
        endif
        print *, j, MatName(j), Ms(j), Gamma(j), K2p(j)
      enddo
      deallocate (check)
      allocate (mcheck(NumMat,NumMat))
      allocate (Exchange(NumMat,NumMat))
      do j=1,NumMat
        do k=1,NumMat
          mcheck(j,k)=.false.
        enddo
      enddo
C     Look for the required intermaterial exchange constants
      do i=1, MAXNMAT*(MAXNMAT-1)/2
        read (3,*,END=31) TMatName, TMatName2, TExchange
        TmatName=STRLOWCASE(ADJUSTL(TMatName))
        TmatName2=STRLOWCASE(ADJUSTL(TMatName2))
        do j=1,NumMat
          do k=j,NumMat
            if((TMatName.eq.MatName(j).and.TMatName2.eq.MatName(k)).or.
     2      (TMatName.eq.MatName(k).and.TMatName2.eq.MatName(j))) then
              if(mcheck(j,k)) then
                print *, 'Datain: duplicated exchange in ', MatTableFN
                stop
              endif
              Exchange(j,k)=Texchange
              Exchange(k,j)=Texchange
              mcheck(j,k)=.true.
            endif

          enddo
        enddo
      enddo
C     Too many materials in the second part of the table
      print *, 'Datain: too many exchanges in ', MatTableFN
      stop
  31  continue
      do j=1,NumMat
        do k=j,NumMat
          if (.not.mcheck(j,k)) then
            print *, 'Datain: missing required exchange in ', MatTableFN
            print *, 'Exchange required: ', MatName(j), MatName(k)
            stop
          endif
c         print *, 'coppia', j, k, Matname(j), Matname(k)
c         print *, 'scambio:', Exchange(j,k)
c         print *
        enddo
      enddo
      deallocate (mcheck)
C     ------------------------------------------------------------------
C     Conversions
      H0angle=H0angle*PI/180.
      H0polar=H0polar*PI/180.
      fmax=fmax*1.E09*2*PI
      fmin=fmin*1.E09*2*PI
      Cellsize=Cellsize*100.
      ZCellsize=ZCellsize*100.
C     ----------------------------------------------------------------
      if (Sel.ne.0) then
C       Open the output file
        if (unformatted) then
          forms='UNFORMATTED'
        else
          forms='FORMATTED'
        endif
        open(UNIT=1,FILE=Filename,STATUS='UNKNOWN',FORM=forms,
     2 ERR=23,IOSTAT=Isterr)
      endif
      return
  11  print *, 'Datain:', infilename, ': error reading line ', Line
      print *, 'or, perhaps, you are reading an old version of input'//
     2 ' file!'
      stop
  12  print *, 'Datain: error reading ',OmfFilename,' line ', Line
      stop
  13  print *, 'Datain: error opening ',OmfFilename, ' ',Isterr
      stop
  14  print *, 'Datain: error opening ',MatTableFN, ' ',Isterr
      stop
  22  print *, 'Datain: cannot open input file, status:', Isterr
      stop
  23  print *, 'Datain: error opening output file, status: ', Isterr
      stop
      end subroutine
      end module
