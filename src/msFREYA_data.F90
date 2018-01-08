! Copyright (c) 2016, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory
! Written by Ramona Voga <vogt2@llnl.gov>, 
!            Jørgen Randrup <jrandrup@lbl.gov>, 
!            Christian Hagmann <hagmann1@llnl.gov>, 
!            Jérôme Verbeke <verbeke2@llnl.gov>.
! LLNL-CODE-701993.
! Title: Fission Reaction Event Yield Algorithm (FREYA), Version: 2.0
! OCEC-16-151
! All rights reserved.
! 
! This file is part of FREYA, Version:  2.0. For details, see <http://nuclear.llnl.gov/simulations>. 
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
! following conditions are met:
! 
! o Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
! 
! o Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer 
! (as noted below) in the documentation and/or other materials provided with the distribution.
! 
! o Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived 
! from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
! DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
! IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! Additional BSD Notice
! 
! 1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was 
! produced at Lawrence Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
! 
! 2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes 
! any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness 
! of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned 
! rights.
! 
! 3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer 
! or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States 
! Government or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed herein do not necessarily 
! state or reflect those of the United States Government or Lawrence Livermore National Security, LLC, and shall not be used 
! for advertising or product endorsement purposes.
!
      module freyaparameters
!     parameters
!     ----------
      implicit none
      save

      integer mxK            ! Number of cases
      integer Mxth
      parameter (Mxth=3)     ! Up to Mxth prefission neutrons
      integer MxEGn
      parameter (MxEGn=1001) ! Max entries in table of Gn/(Gn+Gf)
      integer Nk
      parameter (Nk=1)       ! Number of ejectile types (photons:0)
                             ! Ejectile types (in addition to photons)
      integer mxE
      parameter (mxE=300)    ! Maximum number of energy bins (of width=dE)
      integer Mxg
      parameter (Mxg=2)      ! Maximum number of asymmetric fission modes
      integer mMax
      parameter (mMax=50)    ! Max n+g mult in each chain: pre & post1&2
                             ! Max ejectile multiplicity (incl photons)
      end module freyaparameters

      module freyaG
!     comG common block
!     -----------------
      use freyaparameters
      implicit none
      save

      double precision Gnf
      double precision dEnf
      double precision,dimension(:,:,:),allocatable :: Pnf! Pnf=Gn/(Gn+Gf) [not needed here]
      end module freyaG

      module freyaEv
!     comEv common block
!     ------------------
      implicit none
      save

      integer totEvents, kEv                  ! Event number (just for monitoring)
      end module freyaEv

      module freyaout
!     COMout common block
!     -------------------
      implicit none
      save

      integer L2,L5,L6,L7,L8,L9               ! Output channels
      data L2,L5,L6,L7,L8,L9/2,5,6,7,8,9/
      end module freyaout

      module freyaB
!     comB common block
!     -----------------
!     Thresholds for evaporation and fission
      use freyaparameters
      implicit none
      save

      double precision,dimension(:,:),allocatable :: Sn   ! Neutron separation energy
      double precision,dimension(:,:),allocatable :: Bf   ! Fission barrier
      end module freyaB

      module freyaMth
!     cMth common block
!     -----------------
      implicit none
      save

      integer,dimension(:),allocatable :: Mthk! Max number of pre-fission neutrons
      end module freyaMth

      module freyaK
!     comK common block
!     -----------------
      implicit none
      save

      integer,dimension(:),allocatable :: iZk                  ! Z of up to mK species
      integer,dimension(:),allocatable :: iAk                  ! A of up to mK species
      integer,dimension(:),allocatable :: ntkek                ! number of TKE distributions
      double precision :: EnMax                                            ! maximum incident neutron energy
      double precision,dimension(:),allocatable :: Emax                    ! Maximum excitation
      character(len=100),dimension(:),allocatable :: pAfk      ! name of the p(Af) distribution file
      character(len=100),dimension(:),allocatable ::  preprob  ! names of the files giving the probability 
                                                               ! for pre-equillibrium emission
      character(len=100),dimension(:),allocatable :: prespec   ! names of the files giving the pre-
                                                               ! equillibrium emission neutron spectra
      character(len=100),dimension(:,:),allocatable :: tkek    ! names of the total kinetic energy files
      character(len=30),dimension(:),allocatable :: label      ! Label for process
      character(len=10),dimension(:),allocatable :: reactlab   ! Label for reaction ( sf, (n,f), ...)
      integer,dimension(:),allocatable :: react                ! reaction index:
                                                               ! 0:sf
                                                               ! 1:(n,f)
      end module freyaK

      module freyaZ
!     comZ common block
!     -----------------
      implicit none
      save

      integer mdZ
      data mdZ/5/                     ! mdZ*Zdis is max Zf dev from <Zf>
      double precision,dimension(:),allocatable :: Zdis
      end module freyaZ

      module freyaC
!     comC common block
!     -----------------
      implicit none
      save

      double precision esq,c13,c23,c53
      end module freyaC

      module freyacmass
!     cmass common block
!     ------------------
      use freyaparameters
      implicit none
      save

      double precision wA,Dn,DH
      data wA/931.49386/                ! [MeV/c2] Audi & Wapstra: NPA595 (1995) 409
      data Dn/8.0713228/                ! [MeV/c2] Audi & Wapstra, Table A
      data DH/7.2889694/                ! [MeV/c2] Audi & Wapstra, Table A

      end module freyacmass

      module freyaseed
!     seed common block
!     -------------------
      implicit none
      save

      integer iseed
      data iseed/123457/       ! Seed for RAN0
      end module freyaseed

      module freyaconsts
!     consts common block
!     -------------------
      implicit none
      save

      double precision pi,twopi,hbarc,r0
      data pi/3.141592653/
      data twopi/6.283185308/
      data hbarc/197.33/
      data r0/1.2/             ! Nuclear radius constant (Lysekil)
      end module freyaconsts

      module freyakmass
!     kmass common block
!     ------------------
      implicit none
      save

      double precision,dimension(:),allocatable :: wk ! Mass of k-type ejectiles
      double precision,dimension(:),allocatable :: Dk ! Mass defect of k-type ejectiles
      end module freyakmass

      module freyaCD
!     comCD common block
!     ------------------
      implicit none
      save

      double precision,dimension(:,:),allocatable :: cNZ
      double precision,dimension(:,:),allocatable :: dNZ
      character,dimension(:,:),allocatable :: char
      end module freyaCD

      module freyaWg
!     comWg common block
!     ------------------
      implicit none
      save

      double precision dEf                            ! Width of energy bins in compound nucleus
      data dEf/1.0/
      double precision,dimension(:,:,:,:),allocatable :: CnNE
      double precision,dimension(:,:,:,:),allocatable :: FnNE
      double precision,dimension(:,:,:,:),allocatable :: AnNE
      double precision,dimension(:,:,:,:),allocatable :: WnNE
      end module freyaWg

      module freyacROT
!     comcROT common block
!     ------------------
      implicit none
      save

      double precision,dimension(:),allocatable :: ROT ! Fragment moments of inertia
      end module freyacROT


      module freyacSission
!     cSission common block
!     ---------------------
      implicit none
      save

      integer maxN1, maxZ1
      parameter (maxN1=130,maxZ1=70)        ! Max N & Z of light fission fragment
      double precision,dimension(:,:,:,:),allocatable :: Qsf          ! Q-value for spontaneous fission (gs)
!     double precision,dimension(:,:,:,:),allocatable :: eps1gs       ! gs deformation of light fragment
!     double precision,dimension(:,:,:,:),allocatable :: eps2gs       ! gs deformation of heavy fragment
      double precision,dimension(:,:,:,:),allocatable :: eps1sc       ! scission def of light fragment
      double precision,dimension(:,:,:,:),allocatable :: eps2sc       ! scission def of heavy fragment
      double precision,dimension(:,:,:,:),allocatable :: Esciss       ! Interaction energy at scission
      double precision,dimension(:,:,:,:,:),allocatable :: Escissk    ! Interaction energy at scission
      double precision,dimension(:,:,:,:),allocatable :: Q1sciss      ! Def erg of light fragm at scission
      double precision,dimension(:,:,:,:),allocatable :: Q2sciss      ! Def erg of heavy fragm at scission
#if defined(WRITEL6) || defined(WRITEL8)
      double precision,dimension(:,:,:,:),allocatable :: Qsciss       ! Energy available at scission
      double precision,dimension(:,:,:,:,:),allocatable :: Qscissk    ! Energy available at scission
      integer,dimension(:,:),allocatable :: nc            ! Number of channels open for Nth
#endif
      ! double precision,dimension(:,:,:,:),allocatable :: s12sciss   ! Center sep at scission: s12=c1+c2+d
      ! double precision,dimension(:,:,:,:,:),allocatable :: s12scissk! Center sep at scission: s12=c1+c2+d
      end module freyacSission

      module freyaSAMPLE
!     SAMPLE common block
!     -------------------
      implicit none
      save

      LOGICAL,dimension(:),allocatable :: sampleAf    ! Get Af by direct sampling of P(Af)?
      end module freyaSAMPLE

      module freyacYAf
!     cYAf common block
!     -----------------
      implicit none
      save

      double precision,dimension(:,:),allocatable :: YAfk
      double precision,dimension(:,:),allocatable :: YYAfk
      integer,dimension(:),allocatable :: minAfk
      integer,dimension(:),allocatable :: maxAfk
      end module freyacYAf

      module freyacLOOK
!     cLOOK common block
!     ------------------
      implicit none
      save

      LOGICAL LOOK ! /.FALSE./               ! controls test printouts
      end module freyacLOOK

      module freyaQmin
!     comQmin common block
!     --------------------
      implicit none
      save

      double precision Qmin                           ! Min Q-value & ejectile kin erg (numerics)
                                          ! Q-value cutoff where photons dominate over n;
                                          ! can be adjusted but seems to be close to 0.
      data Qmin/0.010/                    ! Minimum accepted neutron kinetic energy
      end module freyaQmin

      module freyaPath
!     ----------------
      implicit none
      save

      integer maxDATAPATH
      logical :: pathset=.FALSE.          ! set to true after call to msfreya_setpath
      parameter (maxDATAPATH=2000)        ! max length of environment variable
                                          ! FREYADATAPATH
      character(len=maxDATAPATH) freyadir ! path to freya data directory

      end module freyaPath

      module freyaparams
!     params common block
!     -------------------
      implicit none
      save

      double precision Emin   ! Minimum allowed neutron kinetic energy
      parameter(Emin=0.010)
      double precision gmin                                       ! the detection threshold below which photons are unrecorded.
                                                      ! softer photons are simulated but not recorded.
      double precision,dimension(:),allocatable :: alevelk        ! alevel
      double precision,dimension(:),allocatable :: xepsk          ! xeps
      double precision,dimension(:),allocatable :: dTKE0k         ! value of dTKE (if no file)
      character(len=100),dimension(:),allocatable :: dTKEfk! names of the files
                                                      ! containing dTKE vs energy
      logical,dimension(:),allocatable :: dTKEe       ! existence of the dTKE files
      integer,dimension(:),allocatable :: iKtodTKEfk  ! mapping between iK index and 
                                                      ! index of (Ein,dTKE) array
      integer,dimension(:),allocatable :: nEin        ! number of pairs (Ein,dTKE)
                                                      ! in the dTKE file
      double precision,dimension(:,:),allocatable :: dTKEk        ! dTKE versus 
      double precision,dimension(:,:),allocatable :: dTKE_Ek      ! energies
      double precision,dimension(:),allocatable :: ck             ! factor to calculate fluctuations of eps1
      double precision,dimension(:),allocatable :: cTSk           ! factor to calculate spin temperature
      double precision,dimension(:),allocatable :: epsmink        ! maximum heat left after statistial photons
      double precision,dimension(:),allocatable :: gmink          ! detection threshold below which photons are unrecorded
      end module freyaparams

      module freyaerrors
!     common block for errors
!     -----------------------
      implicit none
      save

!     static variable storing whether initerrors has already been called
      logical :: freyaerrorsinitialized=.FALSE. 

!
!      Table of error codes
!      --------------------
!      code    meaning
!         1    data file react.dat could not be found
!         2    error code out of bounds
!         3    case-specific error code out of bounds
!         4    data file SEPn.dat could not be found
!         5    data file fisbar.dat could not be found
!         6    data file alevel.dat could not be found
!         7    data file MassMNMS.dat could not be found
!         8    data file MassAudi.dat could not be found
!         9    data file inputparameters.dat could not be found
!        10    data file acorrection.dat could not be found
!        11    data file Zdis.dat could not be found
!        12    data file gaussfit.dat could not be found
!        13    environment variable FREYADATAPATH too long or
!              path passed to msfreya_setpath too long
!        14    data file Decays/nA.dat could not be found
!        15    data file Decays/iAZ.dat could not be found
!        16    data file Decays/nAZ.dat could not be found
!        17    data file Decays/nl.dat could not be found
!        18    data file Decays/lN.dat could not be found
!        19    data file Decays/Ek.dat could not be found
!        20    data file Decays/Fij.dat could not be found
!        21    data file DefTable.dat could not be found
!
      integer maxErrorCodes
      parameter (maxErrorCodes=21)
      character(len=256),dimension(:),allocatable :: errors! list of errors
      integer,dimension(:),allocatable :: severities       ! severity of errors
      integer,dimension(:),allocatable :: errorctr         ! error counter
!
!     Table of case-specific error codes
!     ----------------------------------
!     code    meaning
!        1    incoming neutron energy Einc out of bounds
!        2    excitation energy eps0 out of bounds
!        3    data file *-Af*.dat could not be found
!        4    data file *-Ktot*.dat could not be found
!        5    data file *.xs could not be found
!        6    data file *.PreEq could not be found
!        7    fission fragment mass mistmatch in data file *-Ktot*.dat
!        8    heat exceeds the maximum allowed
!
      integer maxCaseErrorCodes
      parameter (maxCaseErrorCodes=8)
      character(len=256),dimension(:),allocatable :: caseerrors! list of errors
      integer,dimension(:),allocatable :: caseseverities       ! severity of errors
      integer,dimension(:),allocatable :: caseerrorctr         ! error counter
!
      integer nerrors
      integer ncaseerrors

      logical errorflag
      logical casespec
      integer lasterrorcode

      character(len=256) error                        ! error message
      integer ctr                                     ! error counter
      end module freyaerrors

      module freyarng
!     common block for random number generator
!     ----------------------------------------
      use iso_c_binding, only: C_BOOL

      implicit none
      save

      logical (kind=c_bool) :: user_rng=.FALSE.

      end module freyarng

      module freyaSource
!     common block for Source
!     -----------------------

      implicit none
      save

      integer iZ00,iA00,&
              iZ0,iA0
      double precision Einput,En0

      end module freyaSource

      module freyaRIPL
!     common block for RIPL
!     ---------------------

      implicit none
      save

      integer mZ, mxi, Nmx, mxl, mxd
      parameter (mZ=   80)    ! Maximum CHARGE number accommodated
      parameter (mxi=  50)    ! Max allowed Amax(Z)-Amin(Z)+1 for any Z.
      parameter (Nmx=2000)    ! Maximum Number of Nuclei included
      parameter (mxl=50000)   ! Maximum combined number of levels
      parameter (mxd=300000)  ! Maximum combined number of decays
      
      integer nRIPL
      integer, dimension(:), allocatable ::&
         nA,&   ! Number of isotopes for each Z
         iAN,&  ! Baryon number of nucleus #N
         iZN,&  ! Charge number of nucleus #N
         nl,&   ! Number of energy levels in nucleus #N
         lN,&   ! lN(N) points to the last level of nucleus #N,
                ! so lN(N) = sum{n=1 to N} nl(n); lN(0):=0.
                ! and lN(N-1) is one before the 1st level of #N.
         nkN,&  ! Number of decays from a level:
                ! nk(N,l)=nkN(lN(N-1)+l), l=1,..,nl(N)
         nkNl,& ! Last decay from a level (nkNl(0):=0), so
                ! nkNl(#l)=sum{all decays thru level l} nkN(#l).
         nf     ! Combined table of daughter level indices:
                ! nd(N,l,k)=nf(nkNl(lN(N-1)+l-1)+k)
      integer, dimension(:,:), allocatable ::&
         iAZ,&  ! The included A values for specified Z
         nAZ    ! index #N(Z,A) using i=Z & j=A-iAZ(i,1)+1
      double precision, dimension(:), allocatable ::&
         ElN,&  ! The energy of a level:
                ! E(N,l)=ElN(lN(N-1)+l), l=1,..,nl(N)
         Fij,&  ! Combined table of the decay branchings:
                ! Fij(N,l,k)=Fij(nkNl(lN(N-1)+l-1)+k) 
         tlN    ! Half-life of a level (organized as ElN above)

      end module freyaRIPL

      module freyacount
!     common block for count
!     -----------------------

      implicit none
      save

      double precision, dimension(:), allocatable ::&
        countN,&   ! Count number of various kinds of photons
        countE     ! Count energy of various kinds of photons

      end module freyacount

      module freyaTdist
!     common block for Tdist
!     ----------------------

      implicit none
      save

      integer mT
      parameter (mT=40) ! Number of temperature bins

      integer nT
      double precision dT
      double precision, dimension(:), allocatable ::&
        Ti,&            ! temperature bins
        T1,&
        T2,&
        T3

      end module freyaTdist

      module freyalocal
!     common block for local (in subroutine DecayS)
!      ---------------------------------------------

      implicit none
      save

      double precision h0,h1,h2,d0,d1,d2,E0,E1,E2
      integer, dimension(:), allocatable ::&
        isomer  ! # isomers encountered for this nucleus

      end module freyalocal

      module freyaaveE
!     common block for averages in subroutine DecayS
!     ----------------------------------------------

      implicit none
      save

      double precision ave0,ave1,ave2,ave3,Q0s,Q1s,Q2s,Q3s,Q4s,Q5s,Q6s

      end module freyaaveE

      module freyaDecayS
!     common block for subroutine DecayS
!     ----------------------------------

      implicit none
      save

      integer Ncos
      parameter (Ncos=20)     ! Number of bins for ang dist

      integer NQf,mp,minZ,maxZ
      integer, dimension(:), allocatable :: minA,maxA
      double precision, dimension(:), allocatable :: bin,ben,aves
      double precision, dimension(:,:), allocatable :: vars
      double precision ave,aveQf,Sst1,Sst2,Est1,Est2,&
           Syr1,Syr2,Eyr1,Eyr2

      end module freyaDecayS

      module freyags
!     common block for subroutines gsC & gsM
!     --------------------------------------

      implicit none
      save

      integer mxN,mxZ,mxA,jm
      parameter (mxN=210,mxZ=120,mxA=mxN+mxZ,jm=68)

      end module freyags

      module freyaGD
!     common block for subroutines DecayS & GDRw
!     ------------------------------------------

      implicit none
      save

      double precision GDwidth
      data GDwidth/5.0/

      end module freyaGD

      module freyaID
!     common block for particle IDs of fission fragments
!     --------------------------------------------------
      use freyaparameters

      implicit none
      save

      integer id0,id1,id2

      dimension id0(mMax), & ! Ejectile type of prefission neutrons
                id1(mMax), & ! Ejectile type of 1st fission fragment
                id2(mMax)    ! Ejectile type of 2nd fission fragment

      end module freyaID
