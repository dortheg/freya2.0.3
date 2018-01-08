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
!=======================================================================
!               S P E C I F Y   T H E   S Y S T E M :
!=======================================================================
! The fissioning system can be prepared in two different ways:
! 1) Specification of the initial compound nucleus Z, A, E*, S:
!       Z:      Charge number of the initial compound nucleus;
!       A:      Mass number of the initial compound nucleus;
!       E*:     Total excitation (in MeV) of the initial compound
!       nucleus
!               [E* is the combined rotational and statistical energy].
!       S:      Angular momentum (in hbar) of the initial compound
!       nucleus
!               [the orientation of the spin vector S may be
!               specified];
! 2) Neutron-induced: (Z, A-1) and En (MeV) are specified:
!       It is assumed that the incoming neutron hits the nucleus from a
!       user-specified direction (or random direction if specified 
!       direction is (0,0,0) ) with a random impact parameter and is then
!       fully absorbed (so we ignore pre-equilibrium emission for now),
!       yielding the linear & angular momenta as well as the total excitation.

      SUBROUTINE msfreya_event_c(iK,Einc,eps0,PP0 &
      ,iZ1,iA1,PP1,iZ2,iA2,PP2,m,p,id,ndir) &
      bind (C, name="msfreya_event_c_")
!
!      This subroutine serves as an "overloaded" function for
!      msfreya_event where this function's parameters are type integer 
!      (kind=c_int) and double precision (kind=c_double) instead of the type 
!      integer and double precision that function msfreya_event requires,
!
      use freyaparameters
      use iso_c_binding, only: C_INT, C_DOUBLE

      implicit none

      integer (kind=c_int), value :: iK
      real (kind=c_double), value :: Einc,eps0
      real (kind=c_double), dimension(0:4) :: PP0    ! Four-momentum of the fissioning nucleus
      real (kind=c_double), dimension(0:4) :: PP1    ! Four-momentum of the 1st evaporating nucleus
      real (kind=c_double), dimension(0:4) :: PP2    ! Four-momentum of the 2nd evaporating nucleus
      real (kind=c_double), dimension(1:3) :: ndir   ! normalized direction of the incident neutron
                                                    ! (0,0,0) for random
      integer (kind=c_int) :: iZ1,iA1,iZ2,iA2
      integer (kind=c_int) :: m
      real (kind=c_double), dimension(4,3*mMax) :: p
      integer (kind=c_int), dimension(3*mMax) :: id

      integer :: iK_f
      double precision :: Einc_f,eps0_f
      double precision, dimension(0:4) :: PP0_f    ! Four-momentum of the fissioning nucleus
      double precision, dimension(0:4) :: PP1_f    ! Four-momentum of the 1st evaporating nucleus
      double precision, dimension(0:4) :: PP2_f    ! Four-momentum of the 2nd evaporating nucleus
      double precision, dimension(1:3) :: ndir_f   ! normalized direction of the incident neutron
                                       ! (0,0,0) for random
      integer :: iZ1_f,iA1_f,iZ2_f,iA2_f
      integer :: m_f
      double precision, dimension(4,3*mMax) :: p_f
      integer, dimension(3*mMax) :: id_f

      iK_f=int(iK)
      Einc_f=dble(Einc)
      eps0_f=dble(eps0)
      PP0_f=dble(PP0)
      ndir_f=dble(ndir)

      call msfreya_event(iK_f,Einc_f,eps0_f,PP0_f,iZ1_f,iA1_f,PP1_f,&
                         iZ2_f,iA2_f,PP2_f,m_f,p_f,id_f,ndir_f)
      iA1=int(iA1_f,kind=c_int)
      iZ1=int(iZ1_f,kind=c_int)
      iA2=int(iA2_f,kind=c_int)
      iZ2=int(iZ2_f,kind=c_int)
      m=int(m_f,kind=c_int)
      p=real(p_f,kind=c_double)
      id=int(id_f,kind=c_int)
      PP0=real(PP0_f,kind=c_double)
      PP1=real(PP1_f,kind=c_double)
      PP2=real(PP2_f,kind=c_double)

      RETURN
      END

!************************************************************************

      SUBROUTINE msfreya_event(iK,Einc,eps0,PP0 &
      ,iZ1,iA1,PP1,iZ2,iA2,PP2,m,p,id,ndir)
!
! called from msFREYA_main to generate one complete fission event:
!
!OBS:  Prior preparation by msFREYA_setup is required!

! -----------------------------------------------------------------------
! INPUT:
!      iK         Index identifying the isotope considered (iK=1,...,mxK);
!      Einc       neutron energy (if neutron-induced fission)
!      eps0       its EXCITATION energy: eps0=En+Sn for n-induced, zero for sf.
!      PP0(0:4)   its exc energy, momentum, and kinetic energy.
!      OBS:       It is assumed that the nucleus is initially at rest
!      ndir       normalized direction of incident neutron. (0,0,0) for random
! OUTPUT:
!      iZ1,iA1    Charge & mass number of product 1;
!      PP1(0:4)   its exc energy, momentum, and kinetic energy.
!      iZ2,iA2    Charge & mass number of product 2;
!      PP2(0:4)   its exc energy, momentum, and kinetic energy.
!      n          Number of emitted neutrons;
!      p(1:4,n)   their momentum and kinetic energy.
! -----------------------------------------------------------------------
!OBS:   The four-vector convention follows standard FREYA:
!      0                rest mass
!      1-3              momentum
!      4                total energy
! -----------------------------------------------------------------------
! OBS:  Although the excitation energy eps0 is in principle PP0(0)-W0,
!       we transfer it as an explicit parameter for accuracy reasons:
!       since O(eps0)=1 MeV but O(W0)=200 GeV the subtraction reduces
!       the accuracy by five orders of magnitude!
! -----------------------------------------------------------------------

!========================================================================
      use freyaout
      use freyaEv
      use freyaB
      use freyaK
      use freyaMth
      use freyaparams
      use freyaZ
      use freyaC
      use freyacmass
      use freyaconsts
      use freyaseed
      use freyakmass
      use freyaWg
      use freyacSission
      use freyaSAMPLE
      use freyacYAf
      use freyacLOOK
      use freyaerrors
      use freyacROT
      use freyaID
      use freyainterfaces, only: msfreya_reseterrorflag_c

      use iso_c_binding, only: C_INT, C_DOUBLE, C_BOOL

      implicit none

      integer :: iK
      double precision :: Einc,eps0
      double precision :: S00                      ! Angular momentum
      double precision, dimension(0:4) :: PP0      ! Four-momentum of the fissioning nucleus
      double precision, dimension(0:4) :: PP1      ! Four-momentum of the 1st evaporating nucleus
      double precision, dimension(0:4) :: PP2      ! Four-momentum of the 2nd evaporating nucleus
      double precision, dimension(1:3) :: ndir     ! normalized direction of the incident neutron
                                       ! (0,0,0) for random
      integer :: iZ1,iA1,iZ2,iA2
      integer :: m
      double precision, dimension(4,3*mMax) :: p
      integer, dimension(3*mMax) :: id ! Ejectile type     

      double precision, dimension(3) :: SS0        ! Compound spin vector
      double precision, dimension(3) :: SS1
      double precision, dimension(3) :: SS2

      logical exitonerror, errorflagset

!     Ejectiles types (0-6): photon, n, p, d=2H ,t=3H, tau=3He, alpha=4He

      integer minA,maxA
      data minA,maxA/60,190/   ! minA & maxA for any fission fragment

      integer mult0,mult1,mult2,mult
      double precision p0,p1,p2
      dimension mult0(0:Nk),mult1(0:Nk),mult2(0:Nk),mult(0:Nk) & ! Multiplicity
      ,p0(3,  mMax) &                  ! Momenta of n0 prefission neutrons
      ,p1(3,  mMax) &                  ! Momenta of n1 ejectiles from 1st fission fragment
      ,p2(3,  mMax)                    ! Momenta of n2 ejectiles from 2nd fission fragment

      double precision, dimension(3) :: xhat,yhat,zhat

      integer KRASH
      data KRASH/0/

      integer iZ0,iA0,m0,Nth,iE,NthFission,iN1&
      ,iN2,i,j,l,m1,m2,n0,n,k,n1,n2&
      ,iAp1,iAp2,iA,if,iZ,iN0,iA00,iN00,iZ00
      double precision Z0,A0,A,E0,eta,rng,Amin,Amax,dA,Zave,dZmax &
      ,xnormal,Z,dZ,Vc,Q12,Q1,Q2,epsf,Ekin,aA1ave,aA2ave,aAave,eps1 &
      ,aA1,eps2,aA2,e1,e2,alevel,T1,U1,U2,T2,deps1,deps2,D1,D2,gsM,W1 &
      ,W2,phi,costh,sinth,ex,ey,ez,cos,sin,Ekin1,Ekin2 &
      ,prob,aA0,angle,c,sina,cosa,cosphi,sinphi,d &
      ,E1rot,E2rot,P12,R, Rot12,Rotm,Rotp,RotRel,RotTot,S0 &
      ,Sf1,Sf2,S12,S1sq,S1x,S1y,S1z,S2x,S2y,S2z,S2sq,smx,smy &
      ,spx,spy,T,TKE,TS,Tsc,U12 &
      ,var1,var2,w,E0Rot,P00,Q0,SSx,SSy,SSz,Sxyz,cTS,xeps &
      ,alevel0,E00,En0,eps00,pn0,Q00,Q0max,rot00,Smin,Smax,W00 &
      ,Wi,Wn,msFREYA_SEPn
      logical :: SPINS = .FALSE. ! Include spins?

      integer iTKEfidx           ! maps index iK to index of file (Ein,dTKE)
      integer idxerg             ! index of interval where excitation
                                 ! energy of nucleus is

      double precision RRA,e,eeps,hSsq,S,BARRIER

!************************************************************************
      wn = wA + Dn               ! neutron mass
!
      kEv=kEv+1                  ! Event number (just for monitoring)
      if ((Einc.gt.EnMax).and.(react(ik).eq.1)) then
        write(error, 20) trim(label(iK)), Emax(iK), Einc
20      format("FREYA: Maximum incident neutron energy for ",a, &
              " is ", f7.3, " MeV. Einc=",f7.3," MeV.")
        call seterrorcase(1,error)
        RETURN
      endif
      LOOK=.FALSE.; if (kEv.eq.KRASH) LOOK=.TRUE.
!      if (0.0.lt.eps0.AND.eps0.lt.Bf(0,iK)) then
!#ifdef WRITEL6
!        write (L6, 21) eps0,trim(label(iK)), Bf(0,iK)
!#endif
!        iZ1=0; iZ2=0; n=0
!        write(error, 21) eps0,trim(label(iK)),Bf(0,iK)
!21      format("FREYA: The statistical excitation energy (", f8.4,&
!               "MeV) of the compound nucleus for ",&
!               a, " is below the fission barrier ",f8.4,"MeV, ", &
!               "so fission cannot occur.")
!        call seterrorcase(2,error)
!        RETURN
!      endif

!     reset the secondary particle IDs to -1
      id0=-1
      id1=-1
      id2=-1

      cTS=cTSk(iK)
      if (cTS.lt.0.0) then
        SPINS=.FALSE.
#ifdef WRITEL6
        write (L6,*) 'Angular momentum is NOT included =>'
        write (L6,*) 'initial angular momentum will be set to zero'
#endif
        S0=0.0; SSx=0.0; SSy=0.0; SSz=0.0
      else
#ifdef WRITEL6
        write (L6,"('Rescaling TS by the factor',f7.3)") cTS
        if (cTS.eq.0.0) write (L6,*) 'TS=0: Omit spin fluctuations ', &
                                     'at scission'
#endif
        SPINS=.TRUE.
      endif

! Start from the original compound nucleus at rest:
! This ignores the small push induced by the incoming neutron
      iZ0=iZk(iK); Z0=dble(iZ0) ! Charge number
      iA0=iAk(iK); A0=dble(iA0) ! mass number
      iN0=iA0-iZ0
      E0=PP0(0)                  ! total REST energy ("mass")

      IF (Einc.gt.0.0) THEN ! n-INDUCED: =========================
        En0=Einc                              ! Neutron kinetic energy
#ifdef WRITEL6
        write (L6,"('Neutron-induced fission, En =',f8.3,' MeV =>')") &
              En0
#endif
        Wi=iA0*wA+gsM(iZ0,iA0)                ! Target gs mass
        iZ00=iZ0                              !>   Charge & mass of the
        iA00=iA0                              ! >  initial prefission
        iN00=iA00-iZ00                        !>   compound nucleus
! En0:  The kinetic energy of the neutron
!       impinging on the stationary nucleus from a random direction.
        pn0=sqrt(2*wn*En0)                    ! Incident momentum NR
#ifdef WRITEL6
        write (L6,*) pn0
#endif
        pn0=sqrt(En0*(2*wn+En0))              ! Incident momentum
#ifdef WRITEL6
        write (L6,*) pn0
#endif
        E00=sqrt((Wi+wn+En0+pn0)*(Wi+wn+En0-pn0))     ! Compound tot rst mass
        W00=iA00*wA+gsM(iZ00,iA00)                    ! Compound gs mass
        eps00=E00-W00
#ifdef WRITEL6
        write (L6,*) 'E*:',eps00                      ! Compound excitation
#endif
        eps00=En0+msFREYA_SEPn(iK,iZ00,iA00)
#ifdef WRITEL6
        write (L6,*) 'E*:',eps00 ! more accurate?
#endif
        P00=pn0                                       ! Compound momentum
        E=0.5*P00**2/E00                              ! Compound kin erg (NR)
#ifdef WRITEL6
        write (L6,*) '              S      E*     Erot     Q      Ekin'
#endif
! Head-on collision (zero impact parameter):
        S00=0.0; rot00=0.5*S00**2/ROT(iA00)           ! min compound spin
        Q00=eps00-rot00                               ! max compound heat
#ifdef WRITEL6
        write (L6,"(' head-on: ',5f8.3)") S00,eps00,rot00,Q00,E
#endif
! Grazing collision (maximum impact parameter):
        if (SPINS) S00=RRA(iA00)*pn0/hbarc            ! max compound spin
        rot00=0.5*S00**2/ROT(iA00)                    ! compound rot erg
        Q00=eps00-rot00                               ! min compound heat
#ifdef WRITEL6
        write (L6,"(' grazing: ',5f8.3)") S00,eps00,rot00,Q00,E
#endif
      ELSE ! SPECIFIED E* & S: =========================
        iZ00=iZ0                              !>   Initial
        iA00=iA0                              ! >  compound 
        iN00=iA00-iZ00                        !>   nucleus.
        eps00=abs(Einc)                       !    Excitation (incl rot)
#ifdef WRITEL6
        write (L6,"(' Fission of a nucleus having E* =',f7.3,' MeV')") &
              eps00
#endif
        En0=eps00+gsM(iZ00,iA00)-Dn-gsM(iZ0,iA00-1)! Equiv n energy
#ifdef WRITEL6
        write (L6,"(' Equivalent head-on neutron En:',f9.4,' MeV')") En0
#endif
        Smax=0.0; if (SPINS) Smax=sqrt(2.0*eps00*ROT(iA00))
#ifdef WRITEL6
        write (L6,"('Maximum allowed angular momentum S:',f7.3)") Smax
#endif
        Q0max=EnMax+msFREYA_SEPn(iK,iZ00,iA00)
        Smin=0.0
        if (eps00.gt.Q0Max) Smin=sqrt(2*ROT(iA00)*(eps00-Q0max))
#ifdef WRITEL6
        write (L6,"('Minimum allowed angular momentum S:',f7.3)") Smin
#endif
        if (Einc.eq.0.0) then
#ifdef WRITEL6
          write (L6,*) 'Spontaneous fission'
#endif
!         There can be no angular momentum for E*=0:
          S00=0.0; SSx=0.0; SSy=0.0; SSz=0.0
        else
          if (SPINS) then
    3       continue
#ifdef WRITEL6
            write (L6,*) 'Specify the angular momentum S,sx,sy,sz:'
            write (L6,*) '   S: The magnitude of the angular ',&
                         'momentum \S;'
            write (L6,*) '   (sx,sy,sz): the unnormalized direction ', &
                         'of \S;'
            write (L6,*) '   (sx,sy,sz) = (0,0,0) => random direction.'
            write (L6,*) 'Enter S,sx,sy,sz:'
#endif
#ifdef WRITEL5
            read  (5,*) S00,SSx,SSy,SSz
#endif
            if (S00.gt.Smax) then
#ifdef WRITEL6
              write (L6,*) 'S cannot exceed',Smax,'!'
#endif
              goto 3
            endif
          endif
        endif
        Sxyz=sqrt(SSx**2+SSy**2+SSz**2)
        if (Sxyz.ne.0.0) then
              SSx=SSx/Sxyz; SSy=SSy/Sxyz; SSz=SSz/Sxyz!  normalize
        endif
        rot00=0.5*S00**2/ROT(iA00)                    ! Compound rot erg
        Q00=eps00-rot00                               ! Compound heat
        W00=iA00*wA+gsM(iZ00,iA00)                    ! Compound gs mass
        E00=W00+rot00+Q00                             ! Compound rest mass
        P00=0.0                                       ! Compound is at rest
        if (Q00.gt.Q0max) then
#ifdef WRITEL6
          write (L6, 25) Q00-Q0max
25        format("The heat exceeds the max allowed by:", f8.3, " MeV!")
#endif
          write(error, 24) trim(label(iK)), Q00-Q0max
24        format("FREYA: The heat exceeds the max allowed for ",a,&
                " by:", f8.3, " MeV.")
          call seterrorcase(8,error)
          RETURN
        endif
      ENDIF   !  =========================

! Calculate the angular momentum of the initial compound nucleus:
      S=S00                      ! (maximum) spin magnitude
      IF (Einc.gt.0.0) THEN      ! neutron-induced:
        if (ndir(1).eq.0.and.ndir(2).eq.0.and.ndir(3).eq.0) then
          ! Total momentum PP0 (random direction):
          costh=2*rng(iseed)-1.0; sinth=sqrt(1.0-costh**2)
          phi=twopi*rng(iseed); cosphi=cos(phi); sinphi=sin(phi)
          PP0(1)=P00*sinth*cosphi
          PP0(2)=P00*sinth*sinphi
          PP0(3)=P00*costh
        else
          ! Total momentum PP0 (user-specified direction):
          PP0(1)=P00*ndir(1)
          PP0(2)=P00*ndir(2)
          PP0(3)=P00*ndir(3)
          costh=PP0(3)/P00       ! P00>0 because Einc>0
          sinth=sqrt(1.-costh*costh)
          if (sinth.eq.0) then ! either pole
            sinphi=1.
            cosphi=0.
          else
            sinphi=PP0(2)/P00/sinth
            cosphi=PP0(1)/P00/sinth
          endif
        endif
! Random impact parameter n=(x,y,z) (normalized to unity, perp to
! PP0):
        angle=twopi*rng(iseed)
!       x= cos(angle)*costh*cosphi-sin(angle)*sinphi
!       y= cos(angle)*costh*sinphi+sin(angle)*cosphi
!       z=-cos(angle)*sinth
! Resulting angular momentum SS0 = S * n x PP0:
        S=S*sqrt(rng(iseed))    ! spin magnitude/direction:
        SS0(1)= S*(sin(angle)*costh*cosphi+cos(angle)*sinphi) ! Sx
        SS0(2)= S*(sin(angle)*costh*sinphi-cos(angle)*cosphi) ! Sy
        SS0(3)=-S* sin(angle)*sinth                           ! Sz
      ELSE                       ! specified:
        PP0(1)=0.0; PP0(2)=0.0; PP0(3)=0.0      
        if (Sxyz.eq.0.0) then    ! orient the spin randomly:
          costh=2*rng(iseed)-1.0;
          sinth=sqrt(1.0-costh**2)
          phi=twopi*rng(iseed)
          SS0(1)=S*sinth*cos(phi)
          SS0(2)=S*sinth*sin(phi)
          SS0(3)=S*costh
        else                     ! use the specified spin direction (SSx,SSy,SSz):
          SS0(1)=S*SSx; SS0(2)=S*SSy; SS0(3)=S*SSz
        endif
      ENDIF
      E0Rot=0.5*S**2/ROT(iA0)    ! Initial rotational energy
      Q0=eps0-E0Rot              ! Statistical energy Q (heat)
!     if (Q0.lt.0.0) pause 'Q<0!'
      PP0(0)=E0                  ! Initial rest energy

      m0=0                       ! There are no pre-fission ejectiles yet
      m1=0                       ! There are no ejectiles from fragment 1 yet
      m2=0                       ! There are no ejectiles from fragment 2 yet

!========================================================================
! Consider pre-fission neutron emission:
!
!       Nth (for N'th) denotes the number of emitted pre-fission neutrons,
!       so that Nth=0 indicates that there was no pre-fission emission
!       and Nth>0 indicates the number of pre-fission neutrons emitted.
!       [The notation "Nth" is supposed to mean the "N'th" neutron.]
! NOTE:         It is possible for the compound system to emit a number of 
!       neutrons and then become too de-excited to subsequently fission;
!       such events (which will further de-excite by photon emission) are
!       also characterized by a positive value of Nth, the number of emitted
!       neutrons, even though these events don't lead to fission (but "fuse").
!               In the present implementation of FREYA, such no-fission 
!       events are skipped (though counted by the integer value NthFission).
!       Thus the total number of fission events generated will equal the
!       specified sample size totEvents, while the combined number of fusion events
!       will be given by the variable NULL.  For future applications of FREYA
!       it may be of interest to include also these no-fission (or "fusion")
!       events into the analysis since they typically emit some neutrons.
!
!------------------------------------------------------------------------
! Consider pre-fission emission:
      do j=0,Nk                  ! max allowed ejectiles of specie k:
        mult0(j)=0               ! Source 0 may emit no photons
      enddo
      mult0(1)=Mthk(iK)          ! Source 0 may emit up to Mth neutrons
      if (eps00.lt.BARRIER(iK,iZ00,iA0)) mult0(1)=0
!     .........................................................
      call msFREYA_reseterrorflag_c()
      Nth = NthFission (iK,iZ0,iA0,eps0,SS0,PP0, &
                        mult0,m0,id0,p0)
      if (errorflagset().and.exitonerror()) RETURN
!       .........................................................
!!!	if (LOOK) write (L6,905) kEv,iA0,Nth,mult0(1),m0,eps0,iseed,rng(iseed)
!!!  905  format ('PreF:',i6,i4,i2,'=',i2,'=',i2,':',f10.5,i10,f11.7)
#ifdef WRITEL6
      if (LOOK) write (L6,"(i10,2i5,f10.5,i5,i10)") &
              kEv,iZ0,iA0,eps0,m0,iseed
#endif
      iE=1+int(eps0/dEf)       ! E* bin in the fissioning nucleus (iZ0,iA0)
!!!if (kEv.eq.KRASH) write (L6,"('S:',2i10)") kEv,iseed!!!
!========================================================================
! Choose the fission channel (first A, then Z, then TKE, and finally E*):
   11 IF (sampleAf(iK)) THEN   ! sample P(Af) directly: ********
        eta=rng(iseed); iA1=minAfk(iK)-1
        do while (eta.gt.YYAfk(iA1,iK)) ! OK because eta>0 always.
          iA1=iA1+1
        enddo
        A=dble(iA1)           ! A is needed below for Zave=iZ0*A/iA0
      ELSE                     ! sample 5-gaussian fit: ********
! Choose the particular fission mode n=0,1,2:
        eta=rng(iseed)
        n=0
        do while (eta.gt.FnNE(n,Nth,iE,iK))
          n=n+1
        enddo
! Check:   if (n.gt.Mxg) pause!!
! Choose A from the selected gaussian component:
   12   A=AnNE(n,Nth,iE,iK)+WnNE(n,Nth,iE,iK)*xnormal(iseed)
        Amin=dble(minA); Amax=A0-Amin
        if (A.lt.Amin.OR.A.gt.Amax) goto 12
        iA1=int(A+0.5); dA=A-iA1
        if (rng(iseed).lt.dA) iA1=iA1+1
      ENDIF                                    ! sample P(Af).	********
!!!if (kEv.eq.KRASH) write (L6,"('A:',2i10)") kEv,iseed!!!
!     write (L6,"(i5,' ->',i5,' +',i5)") iA0, iA1, iA0-iA1 
! Choose the charge numbers: (OBS: This must match what is in SCISSION)
      Zave=iZ0*A/iA0                           ! Same Z/A as mother nucleus
      dZmax=mdZ*ZDis(iK)                       !! Max |Z-Zave| included
   13 Z=Zave+ZDis(iK)*xnormal(iseed)           !! Should be constrained:
      if (abs(Z-Zave).gt.dZmax) goto 13        !! like this, for example.
      iZ1=int(Z); dZ=Z-iZ1
      if (rng(iseed).lt.dZ) iZ1=iZ1+1
! Call the lighter fragment #1 (and the heavier #2): Alight = A1 & A2 = Aheavy
      if (iA1.gt.iA0-iA1) then ! switch 1 <-> 2:
        iA1=iA0-iA1
        iZ1=iZ0-iZ1
      endif
      iN1=iA1-iZ1
! Complementary fragment #2 (heavy):
      iA2=iA0-iA1
      iZ2=iZ0-iZ1
      iN2=iA2-iZ2
#ifdef WRITEL6
      if (LOOK) write(L6,"(i8,i2,':',2i4,' ->',2i4,'  & ',2i4,3f10.4)")&
                kEv,iK,iZ0,iA0,iZ1,iA1,iZ2,iA2 &
               ,dble(iZ1)/iA1,dble(iZ2)/iA2,eps0
#endif
! Check that this partition is in the scission table:
      if (dTKEe(iK).eqv..FALSE.) then
        Vc = Esciss (iN1,iZ1,Nth,iK)           ! Energy of scission config
      else
        iTKEfidx = iKtodTKEfk(iK)
        idxerg=0                               ! default initialization
! find the right incident neutron energy
        if (Einc.lt.dTKE_Ek(iTKEfidx,1)) then
! energy lower than all array elements
          idxerg = 1
        else if (Einc.ge.dTKE_Ek(iTKEfidx,nEin(iTKEfidx))) then
! energy greater than all array elements
          idxerg = nEin(iTKEfidx)
        else
! energy between minimum and maximum
          do l=2,nEin(iTKEfidx)
            if (Einc.lt.dTKE_Ek(iTKEfidx,l)) then
              idxerg = l-1
              exit
            endif
          enddo
          prob = (Einc-dTKE_Ek(iTKEfidx,idxerg))/ &
                 (dTKE_Ek(iTKEfidx,idxerg+1)-dTKE_Ek(iTKEfidx,idxerg))
          if (rng(iseed).lt.prob) idxerg = idxerg+1
        endif
        Vc = Escissk (iN1,iZ1,Nth,iTKEfidx,idxerg)
      endif
!x    if (LOOK) write (L6,"(4i4,'=Nth, Vc=',f10.5)") iZ1,iN1,iA1,Nth,Vc
      if (Vc.eq.0.0) goto 11                   ! It isn't, so choose another!

      Q12 = Qsf(iN1,iZ1,Nth,iK) + eps0         ! Q-value for 0* -> 1+2
      Q1 = Q1sciss(iN1,iZ1,Nth,iK)             ! Distortion energy in #1
      Q2 = Q2sciss(iN1,iZ1,Nth,iK)             ! Distortion energy in #2
      epsf=Q12-Vc-Q1-Q2                        ! Total E* at scission
! Check whether this fission channel is available:
      if (epsf.le.0.0) then                    ! it is blocked:
#ifdef WRITEL6
        write (L6,922) kEv,iZ0,iA0,iZ1,iA2,iZ2,iA2,epsf
  922   format(i8,':',2i4,' ->',2i4,' +',2i4, &
               ' is blocked; trying again! E* =',f8.3) 
#endif
        goto 11                                ! Unavailable - try again!
      endif

!==========================================================================
! Change (JMV 2016.05.19), alevel0 was previously not updated to the
! right isotope iK. It did not matter, all alevel0 parameters were
! equal.
      alevel0=alevelk(iK)                
      aA0=dble(iA0)/alevel0          ! Macro level-density param for 0

! CHOOSE THE FRAGMENT 2D SPINS BEFORE ENERGY REDISTRIBUTION:      ---------
      IF (cTS.ge.0.0) THEN            ! select fragment spins Sf1 & Sf2:
!      ---------------------------------------------------------
        Tsc=sqrt(epsf/aA0)            ! Tsc: global temperature at scission
        TS=cTS*Tsc                    ! Spin temperature

! This is the new spin selection which conserves the total angular momentum.
! Convention for the dinuclear coordinate system:
!  z: direction from heavy CM to light CM (dinuclear axis);
!  y: direction along pre-scission compound angular momentum S0;
!  x: direction of relative orbital motion (=> z X x = y: right-handed).
! First define the scission coordinate system {xhat,yhat,zhat}:
        S12=SS0(1)**2+SS0(2)**2; S0=sqrt(S12+SS0(3)**2); S12=sqrt(S12)
        IF (S0.eq.0.0) THEN ! choose a randomly oriented xyz system:
          costh=2.0*rng(iseed)-1.0; sinth=sqrt(1.0-costh**2)
          phi=twopi*rng(iseed); sinphi=sin(phi); cosphi=cos(phi)
          xhat(1)=costh*cosphi; xhat(2)=costh*sinphi; xhat(3)=-sinth
          yhat(1)=-sinphi;      yhat(2)=cosphi;       yhat(3)= 0.0
          zhat(1)=sinth*cosphi; zhat(2)=sinth*sinphi; zhat(3)= costh
          angle=twopi*rng(iseed)
          c=xhat(1); xhat(1)=c*cos(angle)-yhat(1)*sin(angle)
                     yhat(1)=c*sin(angle)+yhat(1)*cos(angle)
          c=xhat(2); xhat(2)=c*cos(angle)-yhat(2)*sin(angle)
                     yhat(2)=c*sin(angle)+yhat(2)*cos(angle)
! or an arbitrary xyz system:
!         xhat(1)=1.0; xhat(2)=0.0; xhat(3)=0.0
!         yhat(1)=0.0; yhat(2)=1.0; yhat(3)=0.0
!         zhat(1)=0.0; zhat(2)=0.0; zhat(3)=1.0
        ELSE    ! choose a coordinate system having y-axis along S0:
          yhat(1)=SS0(1)/S0; yhat(2)=SS0(2)/S0; yhat(3)=SS0(3)/S0 ! ~SS0
          if (S12.eq.0.0) then
            xhat(1)=0.0; xhat(2)=1.0
          else
            xhat(1)=-SS0(2)/S12; xhat(2)=+SS0(1)/S12
          endif
          xhat(3)=0.0
          zhat(1)=xhat(2)*yhat(3)-xhat(3)*yhat(2)
          zhat(2)=xhat(3)*yhat(1)-xhat(1)*yhat(3)
          zhat(3)=xhat(1)*yhat(2)-xhat(2)*yhat(1)
! Then rotate a random angle in the xhat-zhat plane:
          angle=twopi*rng(iseed);  sina=sin(angle); cosa=cos(angle)
          c=xhat(1)*sina+zhat(1)*cosa
          xhat(1) =xhat(1)*cosa-zhat(1)*sina; zhat(1)=c
          c=xhat(2)*sina+zhat(2)*cosa
          xhat(2) =xhat(2)*cosa-zhat(2)*sina; zhat(2)=c
          c=xhat(3)*sina+zhat(3)*cosa
          xhat(3) =xhat(3)*cosa-zhat(3)*sina; zhat(3)=c
! Now zhat is the direction of the dinuclear axis
! and xhat is the direction of the orbital motion.
        ENDIF ! end of choosing (xhat,yhat,zhat) coordinate system.

! Consider only wriggling (+) and bending (-) but not tilting or twisting (~z):
! Do x and y directions separately: 
        d=4.0                                   ! Surface separation
        R=RRA(iA1)+RRA(iA2)+d                   ! Center separation
        w=wA*iA1*iA2/dble(iA0)                 ! Reduced mass (iA0=iA1+iA2)
        RotRel=w*R**2/hbarc**2                  ! Relative moment of inertia
        Rot12=ROT(iA1)+ROT(iA2)                 ! Rot1+Rot2
        RotTot=Rot12+RotRel                     ! Total moment of inertia
        Rotm=ROT(iA1)*ROT(iA2)/Rot12            ! Mom inertia for negative x,y modes
        Rotp=Rot12*RotTot/RotRel                ! Mom inertia for positive x,y modes
! Choose Sm & Sp for x and assign values to S1x & S2x:
   17   smx=sqrt(Rotm*TS)*xnormal(iseed)        ! P(sm) ~ exp(-sm**2/2*Rotm*TS)
        spx=sqrt(Rotp*TS)*xnormal(iseed)        ! P(sp) ~ exp(-sp**2/2*Rotp*TS)
        S1x = spx*ROT(iA1)/Rotp + smx
        S2x = spx*ROT(iA2)/Rotp - smx
! Choose Sm & Sp for y and assign values to S1y & S2y (add rigid rotation):
        smy=sqrt(Rotm*TS)*xnormal(iseed)        ! P(sm) ~ exp(-sm**2/2*Rotm*TS)
        spy=sqrt(Rotp*TS)*xnormal(iseed)        ! P(sp) ~ exp(-sp**2/2*Rotp*TS)
        S1y = spy*ROT(iA1)/Rotp + smy + S0*ROT(iA1)/RotTot
        S2y = spy*ROT(iA2)/Rotp - smy + S0*ROT(iA2)/RotTot
! Choose Sm & Sp for z and assign values to S1z & S2z:
!       smz=sqrt(Rotm *TS)*xnormal(iseed)       ! P(sm) ~ exp(-sm**2/2*Rotm *TS)
!       spz=sqrt(Rot12*TS)*xnormal(iseed)       ! P(sp) ~ exp(-sp**2/2*Rot12*TS)
!       S1z=spz*ROT(iA1)/Rot12+smz
!       S2z=spz*ROT(iA2)/Rot12-smz
! Caution: If the tilting mode is agitated, then the system should be tilted!!
! Coaxial modes are ignored:
        S1z=0.0; S2z=0.0
!       Sx is directed along xhat;
!       Sy is directed along yhat ~ SS0;
!       Sz is directed along zhat ~ R
! Combine spin components:
        do i=1,3
          SS1(i)=S1x*xhat(i)+S1y*yhat(i)+S1z*zhat(i)
          SS2(i)=S2x*xhat(i)+S2y*yhat(i)+S2z*zhat(i)
        enddo
! Calculate the final orbital angular momentum (from now on SS0 is L):
        do i=1,3
          SS0(i)=(SS0(i)-spx*xhat(i)-spy*yhat(i))*RotRel/RotTot
        enddo
        S0=sqrt(SS0(1)**2+SS0(2)**2+SS0(3)**2)
! Calculate fragment spin magnitudes:
        S1sq=S1x**2+S1y**2; Sf1=sqrt(S1sq)      ! Total spin of fragment 1
        S2sq=S2x**2+S2y**2; Sf2=sqrt(S2sq)      ! Total spin of fragment 2
! Class  E1rot=0.5*S1sq/Rot(iA1)                ! Rotational energy of fragm #1
! Class  E2rot=0.5*S2sq/Rot(iA2)                ! Rotational energy of fragm #2
        E1rot=hSsq(Sf1)/ROT(iA1)                ! Rotational energy of fragm #1
        E2rot=hSsq(Sf2)/ROT(iA2)                ! Rotational energy of fragm #2
        if (E1rot+E2rot.ge.epsf) goto 17        ! Erot exceeds epsf, try again!
        epsf=epsf-E1rot-E2rot                   ! Total statistical erg (heat)

#ifdef WRITEL6
        if (epsf.le.0.0) then
          write (L6,*) 'epsf<0:',epsf,iA1,iZ1
!         pause 'E*<0!'
        endif
#endif

! Collect spin information:
! Calculate transverse (xy) components of L & V:
! cL     xL=-S1x-S2x; xV=+yL*R/RotRel           ! Components along x
! cL     yL=-S1y-S2y; yV=-xL*R/RotRel           ! Components along y
!               ---------------------------------------------------------
      ELSE                                      ! put fragment spins Sf1 & Sf2 to zero:
!               ---------------------------------------------------------
! CAUTION:  I am not so sure that this option will work any longer!!
        Sf1=0.0; E1rot=0.0                      ! Fragment #1
        Sf2=0.0; E2rot=0.0                      ! Fragment #2
! cL    xL=0.0; xV=0.0
! cL    yL=0.0; yV=0.0
!               ---------------------------------------------------------
      ENDIF   ! select fragment spins Sf1 & Sf2.
!               ---------------------------------------------------------
! CHOOSE THE FRAGMENT 2D SPINS BEFORE ENERGY REDISTRIBUTION.      ---------

      T=sqrt(epsf/aA0)                ! Tentative temperature at scission

! CHOOSE THE EXCITATION ENERGIES OF THE INDIVIDUAL FISSION FRAGMENTS:
        Ekin=Vc                                ! Total ff kinetic energy 

! Consider first macroscopic sharing (in order to calculate alevel):
      aA1ave=dble(iA1)/alevel0                ! Macro level dens param for 1
      aA2ave=dble(iA2)/alevel0                ! Macro level dens param for 2
      aAave=aA1ave+aA2ave
      eps1=epsf*aA1ave/aAave; aA1=aA1ave       ! Approx excitation in 1. Fair-share heat in 1.
      eps2=epsf*aA2ave/aAave; aA2=aA2ave       ! Approx excitation in 2. Fair-share heat in 2.
! Change the split between light & heavy fragments:
! Correct for microscopic structure:
! If doing this, then the deformation dependence should be included:
      e1=eeps(eps1sc(iN1,iZ1,Nth,iK))          ! def at scission := def of gs
      e2=eeps(eps2sc(iN2,iZ2,Nth,iK))          ! def at scission := def of gs
      aA1=alevel(iK,iZ1,iA1,eps1,U1)           ! Micro lev dens param (e1 is ignored)
      aA2=alevel(iK,iZ2,iA2,eps2,U2)           ! Micro lev dens param (e2 is ignored)
      eps1=epsf*aA1/(aA1+aA2)                  !> Share the excitation in
      eps2=epsf*aA2/(aA1+aA2)                  !> proportion to heat capacity
!
      xeps=xepsk(iK)
! Change energy sharing away from equipartition:
! Change the sharing of the average heat (using the parameter xeps):
      eps1=xeps*eps1                           ! Average heat in 1
      eps2=epsf-eps1                           ! Average heat in 2
! Calculate level-density parameters:
      aA1=alevel(iK,iZ1,iA1,eps1,U1)
      aA2=alevel(iK,iZ2,iA2,eps2,U2)
! Calculate fragment temperatures:
      T1=sqrt(U1/aA1)                          !> Fragment temperatures:
      T2=sqrt(U2/aA2)                          !> will generally differ.

! Note on the model assumption:
! ----------------------------------------*
! We assume that the spins are sampled from thermal distributions based
! on the redistributed average fragment excitation energies; the subsequent
! thermal fluctuations in the excitation energies do not affect the spins,
! but are absorbed entirely by the relative fragment motion, i.e. TKE.
! Also note that the compensating relative orbital fragment motion is not
! kept track of separately but is included into the TKE.
! ----------------------------------------------------------------------*

! Choose the thermal fluctuations of the excitation energies:
! When the excitation energy is small, then the truncation of the normal
! distribution leads to a significant reduction in the variance of the
! resulting distribution of eps; this can be remedied by enhancing the
! width of the gaussian sampled, but the enhancement depends on E*
! and it needs to be determined in a more elaborate manner!
!     c=1.18  ! Enhance the thermal fluctuations (to match variance 2UT):
!     c=1.0
      c=ck(iK)
   16 continue

! Choose the thermal fluctuation of eps1:
      var1=2*T1*U1                            ! Canonical variance of eps1
   14 deps1=c*sqrt(var1)*xnormal(iseed)       ! Fluctuation of eps1
      if (abs(deps1).ge.eps1) goto 14         ! We require dis(E*) < ave(E*)
      eps1=eps1+deps1                         ! Heat content of fragment #1
                                              ! Excitation of fragm #2

! Choose the thermal fluctuation of eps2:
! cx    deps2=sqrt(6*U2*T2)*(2*rng(iseed)-1.0)
      var2=2*T2*U2                            ! Canonical variance of eps2
  15  deps2=c*sqrt(var2)*xnormal(iseed)       ! Fluctuation of eps2
      if (abs(deps2).ge.eps2) goto 15         ! We require dis(E*) < ave(E*)

      eps2=eps2+deps2                         ! Heat content of fragment #2
                                              ! Excitation of fragm #2
#ifdef WRITEL6
      if (eps2.le.0.0) then
        write (L6,"(i5,':',f10.4)") k,T
        write (L6,"('1:',6f10.4,i5)")eps1-deps1,deps1,eps1,aA1,U1,T1,iA1
        write (L6,"('2:',6f10.4,i5)")eps2-deps2,deps2,eps2,aA2,U2,T2,iA2
!       pause 'eps2<0!'
      endif
#endif

! Correct & complete the fragment kinetic energy:
      Ekin=Ekin-deps1-deps2                    ! Adjust Ekin and check:
      if (Ekin.le.0.0) goto 16                 ! Ekin < 0: Try again!

      eps1=eps1+Q1                             ! Add scission heat
      eps2=eps2+Q2                             ! Add scission heat
      D1=gsM(iZ1,iA1)                          ! Mass excess of #1
      D2=gsM(iZ2,iA2)                          ! Mass excess of #2
      W1=iA1*wA+D1                             ! g-s mass of ff #1 (light)
      W2=iA2*wA+D2                             ! g-s mass of ff #2 (heavy)
      PP1(0)=W1+E1rot+eps1                     ! Total rest erg of ff #1
      PP2(0)=W2+E2rot+eps2                     ! Total rest erg of ff #2
! Calculate the fission fragment momenta:
! The following kinematics is done non-relativistically for simplicity
! as well as numerical accuracy (since Ekin << Etot):
      P12=sqrt(2*Ekin/(1.0/PP1(0)+1.0/PP2(0)))! Fragment momentum magnitude
! Without orbital rotation, the fragments emerge along R ~ zhat:
      do i=1,3
        PP1(i)=P12*zhat(i); PP2(i)=-PP1(i)
      enddo

      IF (S0.ne.0.0) THEN     ! Coulomb rotation:
! Construct the exit xyz coordinate system:
! The z axis remains unchanged: along the dinuclear axis.
! The y axis is along the exit orbital angular momentum:
        do i=1,3
          yhat(i)=SS0(i)/S0               ! exit y axis
        enddo
        xhat(1)=yhat(2)*zhat(3)-yhat(3)*zhat(2) !>      The orbital velocity
        xhat(2)=yhat(3)*zhat(1)-yhat(1)*zhat(3) ! >     is along this axis
        xhat(3)=yhat(1)*zhat(2)-yhat(2)*zhat(1) !>      at the scission time.

        U12=P12/w                               ! Relative fragment speed
        angle=twopi/4-atan(esq*iZ1*iZ2/(S0*hbarc*U12))
        sina=sin(angle); cosa=cos(angle)
        do i=1,3
          PP1(i)=(xhat(i)*sina+zhat(i)*cosa)*P12; PP2(i)=-PP1(i)
        enddo
      ENDIF                   ! Coulomb rotation.

      ex=PP1(1)/P12; ey=PP1(2)/P12; ez=PP1(3)/P12     ! Direction of fragm 1
      Ekin1=iA2*Ekin/iA0; PP1(4)=PP1(0)+Ekin1 ! Total energy of fragment 1
      Ekin2=iA1*Ekin/iA0; PP2(4)=PP2(0)+Ekin2 ! Total energy of fragment 2
      TKE=Ekin1+Ekin2                         ! Total fragment kinetic energy

!     C12d=s12sciss(iN1,iZ1,Nth)

!=======================================================================
#ifdef WRITEL6
      if (LOOK) write (L6,"('C:',2i10,2f10.5)") k,iseed,eps1,eps2
#endif

      do if=1,2           ! Begin loop over both fission fragments:
!       =====================================
        if (if.eq.1) then ! 1st decay chain:
          iZ=iZ1
          iA=iA1
          mult1(0)=-1     ! No restriction on photon emission (=0 => no photons)
          mult1(1)=-1     ! No restriction on neutron emission
!          mult1(0)=0     ! no photon emission!!
!          mult1(1)=0     ! no neutron evaporation!!
          id1(mMax)=1
          call msFREYA_reseterrorflag_c()
          call DecayS(iK,1,iZ,iA,eps1,SS1,PP1,mult1,m1,id1,p1)
          if (errorflagset().and.exitonerror()) RETURN
          iAp1=iA         ! Mass number of product nucleus #1
        else              ! 2nd decay chain:
          iZ=iZ2
          iA=iA2
          mult2(0)=-1     ! No restriction on photon emission
          mult2(1)=-1     ! No restriction on neutron emission
!          mult2(0)=0     ! no photon emission!!
!          mult2(1)=0     ! no neutron evaporation!!
          id2(mMax)=2
          call msFREYA_reseterrorflag_c()
          call DecayS(iK,1,iZ,iA,eps2,SS2,PP2,mult2,m2,id2,p2)
          if (errorflagset().and.exitonerror()) RETURN
          iAp2=iA         ! Mass number of product nucleus #1
        endif
!       =====================================
      enddo               ! End of loop over both fission fragments.



! OLDER CODE TO BE REMOVED
!  !     ====================================     Decay chain #1:
!        mult1(0)=-1      ! No restriction on photon emission (=0 => no photons)
!        mult1(1)=-1      ! No restriction on neutron emission
!        m1=0             ! There are no ejectiles from fragment 1 yet
!  ! cancel  id1(mMax)=1     what is it for?
!  !!!if (kEv.eq.KRASH) write (L6,"('D:',2i5,f10.5,i10)") iZ1,iA1,eps1,iseed
!        call DDECAY (iK,kEv,iZ1,iA1,eps1,PP1,mult1,m1,id1,p1)
!  !!!if (kEv.eq.KRASH) write (L6,"('D:',2i5,f10.5,i10)") iZ1,iA1,eps1,iseed
!  !     =====================================    Decay chain #2:
!        mult2(0)=-1      ! No restriction on photon emission (=0 => no photons)
!        mult2(1)=-1      ! No restriction on neutron emission
!  ! cancel  id2(mMax)=2     what is it for?
!        m2=0             ! There are no ejectiles from fragment 2 yet
!  !!!if (kEv.eq.KRASH) write (L6,"('D:',2i5,f10.5,i10)") iZ2,iA2,eps2,iseed
!        call DDECAY (iK,kEv,iZ2,iA2,eps2,PP2,mult2,m2,id2,p2)
!  !      =====================================

! Combine all three ejectile classes:   ---------------------------------
      n0=0                      ! Number of ejectiles from source 0 that are neutrons:
      do m=1,m0                 ! Pre-fission ejectiles:
        if (id0(m).eq.1) n0=n0+1! Ejectile is a neutron
        id(m)=id0(m)            ! Ejectily identity
        s=0.0
        do j=1,3
          p(j,m)=p0(j,m); s=s+p0(j,m)**2
        enddo
      enddo
      id0(m0+1)=-1    ! reset particle ID
#ifdef WRITEL6
      if (n0.ne.mult0(1)) &
        write (L6,*) k,': 0 mult mismatch:',n0,mult0(1)
#endif
! Check neutron multiplicity from fragment #1:
      n1=0                      ! Number of ejectiles from source 1 that are neutrons:
      do m=1,m1                 ! Ejectiles from the LIGHT fragment:
        if (id1(m).eq.1) n1=n1+1! This ejectile is a neutron
        id(m0+m)=id1(m); s=0.0
        do j=1,3
          p(j,m0+m)=p1(j,m); s=s+p1(j,m)**2
        enddo
      enddo
      id1(m1+1)=-1    ! reset particle ID
#ifdef WRITEL6
      if (n1.ne.mult1(1)) &
        write (L6,*) k,': 1 mult mismatch:',n1,mult1(1)
#endif
! Check neutron multiplicity from fragment #2:
      n2=0                      ! Number of ejectiles from source 2 that are neutrons:
      do m=1,m2                 ! Ejectiles from the HEAVY fragment:
        if (id2(m).eq.1) n2=n2+1! This ejectile is a neutron
        id(m0+m1+m)=id2(m); s=0.0
        do j=1,3
          p(j,m0+m1+m)=p2(j,m); s=s+p2(j,m)**2
        enddo
      enddo
      id2(m2+1)=-1    ! reset particle ID
#ifdef WRITEL6
      if (n2.ne.mult2(1)) &
        write (L6,*) k,': 2 mult mismatch:',n2,mult2(1)
#endif
      n=n0+n1+n2                ! Total neutron multiplicity = mult(1):
      m=m0+m1+m2                ! Total ejectile multiplicity
      do j=0,1                  ! Nk
        mult(j)=mult0(j)+mult1(j)+mult2(j)
      enddo

! Neutron multiplicity & energy versus TKE:
#ifdef WRITEL6
      if (iA2.lt.iA1) write (L6,*) 'A2<A1:',k  ! Check that 1 is lightest
#endif

! Calculate all neutron and gamma kin. energies for use in fission library:
      do i=1,m
        if (id(i).eq.1) then
          p(4,i) = 5.344305e-4 *(p(1,i)**2+p(2,i)**2+p(3,i)**2)  !  non-rel
        else
          p(4,i) = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
        endif
      enddo


! #ifdef NEUTRONSONLY
! ! Collect only the neutrons:
!       n=0
!       do i=1,m
!         if (id(i).eq.1) then
!           n=n+1
!           do j=1,4
!             p(j,n)=p(j,i)
!             id(n)=1
!           enddo
!         endif
!       enddo
! #else
! ! Collect all ejectiles (neutrons and photons):
!       n=m
! #endif

#ifdef WRITEL6
      if (LOOK) write (L6,"('X:',i10)") iseed!!!
#endif
      if (LOOK) then
        ! pause!!!
      endif
      RETURN
      END

!************************************************************************
      FUNCTION NthFission (iK,iZ0,iA0,eps0,SS,PP,mm,m,id,p)
! considers pre-fission neutron emission ("N'th-chance fission"),
! first one possible pre-equilibrium neutron, then sequential evaporation.
! It generates one event that is ready to fission
! after having evaporated at most mm(1)=Mth neutrons.
!  INPUT:       -------------------------
! iK:           Index identifying the isotope considered (iK=1,...,mxK);
! iA0,iZ0,eps0  A & Z & E* of the initial (pre-emission) nucleus;
! SS(1:3)       Angular momentum of the initial nucleus;
! mm(0:Nk)      Max multiplicity of each ejectile type; 
!   mm(k).GE.0: At most that many type-k emissions will be admitted;
!   mm(k).LT.0: There are no restrictions on type-k emissions.
! m:            Controls PreEq emission: m<0 => is turned OFF, else is ON.
! Note:         To carry out at most one neutron emission,
!               specify mm(1)=1 and mm(k.NE.1)=0 at entry.??
!  OUTPUT:      -------------------------
!NthFission:    The number of neutrons emitted before fission (=Nth=m)
! iA0,iZ0,eps0  A & Z & E* of the post-emission nucleus, ready to fission;
! SS(1:3)       Angular momentum of the fissioning nucleus;
! PP(0:4)       Mass (0), momentum (1:3), energy (4) of fissioning nucleus
! mm(0:Nk)      Actual multiplicity of each ejectile type k=0:Nk
! m             Total number of ejectiles this event: m = mm(0) +...+ mm(Nk);
! id(1:m)       their type (id = 0,...,Nk) [only id = 0(g) & 1(n) here!];
! p(1:3,1:m)    their cartesian momentum components.
! Called from FREYA main code as 
!       Nth = NthFission (k,iZ0,iA0,eps0,PP0,mult0,m0,id0,p0)
!       where Nth becomes the number of pre-fission neutrons emitted.
!
!========================================================================
! Nth:  Number of pre-fission neutrons emitted (Nth = 0,...,Mth):
!       Nth (for nth) denotes the number of emitted pre-fission neutrons,
!       so that Nth=0 indicates that there was no pre-fission emission
!       and Nth>0 indicates the number of pre-fission neutrons emitted.
!       [The notation "Nth" is supposed to mean the "N'th" neutron.]
! OBS:       The code is written so work properly alsofor Mth=0, where it serves
!       to initialize PP(0:4) and put mm(1):=0.
! NOTE:         It is possible for the compound system to emit a number of 
!       neutrons and then become too de-excited to subsequently fission;
!       such events (which will further de-excite by photon emission) are
!       also characterized by a positive value of Nth, the number of emitted
!       neutrons, even though these events dont lead to fission (but "fuse").
!               In the present implementation of FREYA, such no-fission 
!       events are skipped (though counted by the integer value NthFission).
!       Thus the total number of fission events generated will equal the
!       specified sample size totEvents, while the combined number of fusion events
!       will be given by the variable NULL.  For future applications of FREYA
!       it may be of interest to include also these no-fission (or "fusion")
!       events into the analysis since they typically emit some neutrons.
!
!       NthFission has been augmented by a call to PreEq which may lead to 
!       the emission of one pre-equilibrium neutron (Apr2010).
!
!========================================================================
      use freyaout
      use freyaEv
      use freyaB
      use freyaG
      use freyacmass
      use freyaconsts
      use freyaseed
      use freyacLOOK
      use freyaQmin
      use freyacROT
      use freyaparams
      use freyainterfaces, only: msfreya_reseterrorflag_c

      implicit none

      integer NthFission
      integer iK,iZ0,iA0,mm,m,id
      double precision eps0,PP,SS,p
!      double precision eps0,PP,SS,ROT,p
      dimension id(mMax)           ! Identities of up to mMax/1 ejectiles
      dimension SS(3)              ! Angular Momentum
!      dimension ROT(300)           ! Moments of inertia

!     Channel counters
      integer,save::null
      integer,save::more
      integer,dimension(:,:),allocatable,save::kf ! # fissions after emission of Nth neutrons
      integer,dimension(:,:),allocatable,save::kn ! # neutrons before fission or fusion = Nth
      integer,dimension(:,:),allocatable,save::kg ! # fusions after emission of Nth neutrons
      integer,dimension(:),allocatable,save::kx   ! # fusions for specie iK (these are replaced)

      dimension &
       mm(0:Nk) &               ! Multiplicity of each ejectile type
      ,p(3,mMax) &              ! Momenta of up to mMax ejectiles
      ,PP(0:4)                  ! Momentum & energy of the nucleus (0:Mass)

      integer Mth,Nth,iA,n0,n,i1
      double precision D0,gsM,W0,E0,eps,Pf,rng,d1,d2,Pn,eta
      double precision S2,Erot,Q

      logical exitonerror, errorflagset

! Caution:  There is yet no check to ensure that m.LE.mMax=30!
!       This is probably OK during the pre-fission stage.

       NthFission=0             ! # neutrons emitted before fission occurs (default initialization)

       if (.not.allocated(kf)) then
         allocate(kf(0:Mxth+1,mxK))
         kf=0
       endif
       if (.not.allocated(kn)) then
         allocate(kn(0:Mxth+1,mxK))
         kn=0
       endif
       if (.not.allocated(kg)) then
         allocate(kg(0:Mxth+1,mxK))
         kg=0
       endif
       if (.not.allocated(kx)) then
         allocate(kx(mxK))
         kx=0
       endif

! Start from the original compound nucleus at rest:
#ifdef WRITEL6
       if (LOOK) write (L6,"(i9,': Entering NthF:'2i5,f10.5,2i10)") &
                       kEv,iZ0,iA0,eps0,mm(1),iseed
#endif
       D0=gsM(iZ0,iA0)           ! Mother mass excess
       W0=iA0*wA+D0              ! Mother g.s. mass
       E0=W0+eps0                ! Mother total mass
       Mth=mm(1)                 ! Max # pre-fiss neutrons
       mm(1)=0                   ! Actual # pre-fiss neutrons

       PP(0)=E0                                        ! Rest energy (=mass*)
       PP(4)=sqrt(PP(0)**2+PP(1)**2+PP(2)**2+PP(3)**2) ! Total energy

!      PP0=E0;   PP1=PP(1); PP2=PP(2); PP3=PP(3); PP4=PP(4)
!  10  PP(0)=E0; PP(1)=PP1; PP(2)=PP2; PP(3)=PP3; PP(4)=PP4
       Nth=0   !    -> Number of neutrons emitted before fission/fusion
       iA=iA0; eps=eps0; n0=0
       if (Mth.eq.0) goto 50

! Consider pre-equilibrium emission of up to ONE neutron (gives Nth=0 or 1):
       call msFREYA_reseterrorflag_c()
       call PreEq (iK,iZ0,iA,eps,PP,Nth,id,p) ! Pre-Equilibrium
       if (errorflagset().and.exitonerror()) RETURN
       n0=n0+Nth                 ! # pre-eq this event (0 or 1)
! Cancel: goto 50 ! to include Pre-Eq ONLY, without any subsequent Pre-Fiss evap!

!------------------------------------------------------------------------
! The pre-fission evaporation procedure has been modified
! The new (iterative) procedure is as follows:
!    ->   Is evaporation possible (from this nucleus at this energy)?
!            NO:  -> this nucleus is fission ready
!         Would fission be possible from the evaporation daughter?
!            NO:  -> this nucleus is fission ready
!         Consider competition between evaporation and fission
!         Calculate Pn=Pnf(Nf): probability for evaporation
!            if eta > Pn -> this nucleus is fission ready
!         Carry out one evaporation using evapS1n (Nth->Nth+1) 
!         Repeat the procedure
! Thus fission will always occur
!------------------------------------------------------------------------

       DO WHILE (Nth.lt.Mth)  ! Emit up to Mth neutrons:        ========
         S2=SS(1)**2+SS(2)**2+SS(3)**2
         Erot=0.5*S2/ROT(iA0) ! Rotational energy
         Q=eps-Erot           ! Statistical excitation
!-----------------------------------------------------------------------
         if (LOOK) then
#ifdef WRITEL6
           write (L6,*) kEv,': Checking counts in NthFission:',Nth
#endif
           kg(Mth+1,iK)=0; kn(Mth+1,iK)=0; kf(Mth+1,iK)=0
           do n=0,Mth
             kg(Mth+1,iK)=kg(Mth+1,iK)+kg(n,iK)
             kf(Mth+1,iK)=kf(Mth+1,iK)+kf(n,iK)
             kn(Mth+1,iK)=kn(Mth+1,iK)+kn(n,iK)
           enddo
#ifdef WRITEL6
           write (L6,921) 'k,  Nth:',kEv,(n,n=0,Mth)
           write (L6,921) 'Fissions:',kf(Mth+1,iK),(kf(n,iK),n=0,Mth)
           write (L6,921) 'Neutrons:',kn(Mth+1,iK),(kn(n,iK),n=0,Mth)
           write (L6,921) 'Compound:',kg(Mth+1,iK),(kg(n,iK),n=0,Mth), &
             null
#endif
!           pause!!
  921      format(1x,a9,i8,':',6i8)
#ifdef WRITEL6
           write (L6,"(i3,':',i5,4f10.5)") &
                 Nth,iA,eps,Bf(Nth,iK),Sn(Nth,iK),eps-Sn(Nth,iK)-Qmin
#endif
         endif
!-----------------------------------------------------------------------
! Is evaporation possible?
!
         IF (Q.lt.Sn(Nth,iK)+Qmin) THEN
! Evaporation is NOT possible:  -----------------------------------------
#ifdef WRITEL6
           if (LOOK) write(L6,*) 'eps-Erot < Sn(Nth)+Qmin: no evap'
#endif
           goto 50                      ! This nucleus is fission ready
! Evaporation is not possible.  -----------------------------------------
         ENDIF

!------------------------------------------------------------------------
! Would fission be possible from the evaporation daughter?
!
         IF (Q-Sn(Nth,iK)-Emin.lt.Bf(Nth+1,iK)) THEN
! Evaporation daughter cannot fission:  ---------------------------------
#ifdef WRITEL6
           if (LOOK) write (L6,"(i3,5f9.5)") &
                           Nth,eps-Sn(Nth,iK),Bf(Nth+1,iK)+Emin
#endif
           goto 50             ! This nucleus is fission ready
! Evaporation daughter cannot fission.  ---------------------------------
         ENDIF
!------------------------------------------------------------------------
! Consider competition between evaporation and fission
! Calculate the evaporation probability Pn=Pnf(Q;Nth):
         i1=int(Q/dEnf); d1=Q-i1*dEnf; d2=dEnf-d1        ! Lookup eps in table
         Pn=(d2*Pnf(i1,Nth,iK)+d1*Pnf(i1+1,Nth,iK))/dEnf ! Calculate Pn(Q;Nth)
#ifdef WRITEL6
         if (LOOK) then
           write (L6,"(2i5,6f10.5)") Nth,i1,Q,d1, &
                 Pnf(i1,Nth,iK),Pn,Pnf(i1+1,Nth,iK),d2
           write (L6,*) i1,d1+d2-dEnf,' =0?',Pn,iseed
         endif
#endif
! Does fission occur?
         eta=rng(iseed)
#ifdef WRITEL6
         if (LOOK) write (L6,*) 'eta=',eta,Pn,iseed
#endif
         if (eta.gt.Pn) goto 50     ! This event is fission ready
!------------------------------------------------------------------------
! Carry out evaporation of one add'l neutron:
#ifdef WRITEL6
         if (LOOK) write (L6,*) iA,Q,Sn(Nth,iK),Nth    ! before evaporation
#endif
!        -----------------------------------------
         call EvapS1n (kEv,iK,iZ0,iA,eps,SS,PP,Nth,p)  ! => Nth -> Nth+1
!        -----------------------------------------
#ifdef WRITEL6
         if (LOOK) write (L6,*) iA,eps,Sn(Nth,iK),Nth  ! after evaporation
#endif
         id(Nth)=1                                     ! This ejectile is a neutron
       ENDDO                                           ! Emit up to Mth neutrons. =========
!      =================================================================

! Check if additional pre-fission evap + fiss would be possible:
       if (Nth.eq.Mth) then                            ! This event might emit more than Mth
         S2=SS(1)**2+SS(2)**2+SS(3)**2                 ! neutrons but fission is enforced now!
         Q=eps-0.5*S2/ROT(iA)                          ! Statistical energy
         if (Q.lt.Sn(Mth,iK)+Qmin) then                ! evaporation is possible:
           if (Q-Sn(Mth,iK)-Emin.lt.Bf(Mth+1,iK)) then ! subseq fiss is possible:
             more=more+1                               ! count such cases
           endif                                       ! subseq fiss is possible.
         endif                                         ! evaporation is possible.
       endif

!=======================================================================
! The nucleus is now done with pre-fission emission and is ready to fission:

   50  iA0=iA                             ! Mass number of the fission-ready nucleus
       eps0=eps                           ! Excitation in the fission-ready nucleus
       mm(1)=mm(1)+Nth                    ! Pre-fission neutron multiplicity. Only n =>
       m=mm(1)                            ! m: Total pre-fission multiplicity / same here
       kx=kx+n0                           ! Count pre-equilibrium neutrons
       do n=1,Nth
         kn(n-1,iK)=kn(n-1,iK)+1          ! Count # prefiss evaporations
       enddo
! OBS:  The counting of ejectiles could be done differently to better prepare
!        for emission of other ejectile species (gammas, protons, ..).
       NthFission=Nth                     ! # neutrons emitted before fission occurs
       kf(Nth,iK)=kf(Nth,iK)+1            ! Count the number of fission events (=totEvents?)
!-----------------------------------------------------------------------
       if (LOOK.OR.kEv.eq.totEvents) then
#ifdef WRITEL6
         write (L6,*) kEv,': Checking final counts in NthFission:',Nth
#endif
         kg(Mth+1,iK)=0; kn(Mth+1,iK)=0; kf(Mth+1,iK)=0
         do n=0,Mth
           kg(Mth+1,iK)=kg(Mth+1,iK)+kg(n,iK)
           kf(Mth+1,iK)=kf(Mth+1,iK)+kf(n,iK)
           kn(Mth+1,iK)=kn(Mth+1,iK)+kn(n,iK)
         enddo
#ifdef WRITEL6
         write (L6,921) 'k,  Nth:',kEv,(n,n=0,Mth)
         write (L6,921) 'Fissions:',kf(Mth+1,iK),(kf(n,iK),n=0,Mth),more
         write (L6,921) 'Neutrons:',kn(Mth+1,iK),(kn(n,iK),n=0,Mth),m
         write (L6,921) 'Compound:',kg(Mth+1,iK),(kg(n,iK),n=0,Mth),null
         write (L6,111)
#endif
       endif
!-----------------------------------------------------------------------
#ifdef WRITEL6
       if (LOOK) write (L6,"(i9,': Leaving NthF:'2i5,f10.5,2i10)") &
                       kEv,iZ0,iA0,eps0,mm(1),iseed
#endif
       if (LOOK) then
         ! pause
       endif
!!!	if (kEv.gt.217399) write (L4,"(i10,':',i10)") kEv,iseed!!!
  111  format (80('-'))

       RETURN
       END        ! NthFission. 

!************************************************************************
      SUBROUTINE msfreya_getids_c (id0_c,id1_c,id2_c) &
      bind (C, name="msfreya_getids_c_")
! called after msfreya_event_c() to retrieve ids of pre-fission particles
! secondaries emitted from first and second fission fragments

      use freyaparameters
      use freyaID
      use iso_c_binding, only: C_INT

      implicit none

      integer (kind=c_int), dimension(mMax) :: id0_c
      integer (kind=c_int), dimension(mMax) :: id1_c
      integer (kind=c_int), dimension(mMax) :: id2_c

      id0_c=int(id0,kind=c_int)
      id1_c=int(id1,kind=c_int)
      id2_c=int(id2,kind=c_int)

      RETURN
      END

!************************************************************************
       SUBROUTINE PreEq (iK,iZ,iA,eps0,PP,m,id,p)
! called from NthFission; the subroutine PreEq
! contains the same arguments as NthFission, but m MUST be set here. PreEq
! considers possible emission of ONE pre-equilibrium neutron (so far), based on
! calculations by Erich Ormand (Aug11 for 240Pu, 239U, 236U).
! Caution:  the array mm(0:Nk) is NOT affected here but is updated in NthFission
!       after the pre-fission emission process is completed; this is done
!       because the process may be aborted and restarted occasionally (namely
!       when more than the specified max number of neutrons would be emitted).
! Input parameter values of PreEq:
! iK:           Index identifying the isotope considered (iK=1,...,mxK);
! iZ, iA        Charge & mass number of combined system (target + neutron)
! eps0          Excitation energy of combined system
! PP(0:4)       Momentum & energy of the combined system (0:Mass)
! m     <0:     Calling PreEq with m<0 initially turns PreEq emission OFF;
!   initial     In the initial call (iK=0) m should be the number of cases mxK;
!   general     m is a variable (incremented by 1 if there is PreEq emiss).
! Output parameter values of PreEq:
! iZ, iA        Charge & mass number of resulting compound nucleus
! eps0          Excitation energy of resulting compound nucleus
! PP(0:4)       Momentum & energy of the resulting compound nucleus (0:Mass)
! mm(0:Nk)      Multiplicity of each ejectile type (updated in NthFission)
! m  initial:   m is not touched, only read;
!    general:   m is incremented by 1 if a pre-equilibrium neutron is emitted;
!               right now PreEq is called with m=0 resultimng in m=0 or 1;
!               extension to multiple pre-equilibrium emission is not done yet.
! id(mMax)      Identities of up to mMax ejectiles (up to m=0 or 1 for now)
! p(3,mMax)     Momenta of up to mMax ejectiles (up to m=0 or 1 for now)
!
! NOTE: The energy argument in PreEq is the nuclear EXCITATION energy E*,
!       whereas the pre-equilibrium probability and the neutron spectra are
!       tabulated in terms of the kinetic energy of the incident neutron En0;
!       the incident energy is En0=eps0-Sn (Sn: neutron separation energy).

!========================================================================
       use freyaout
       use freyacmass
       use freyaconsts
       use freyaseed
       use freyaparameters
       use freyaK
       use freyaMth
       use freyaPath
       use freyaerrors

       use freyainterfaces

       implicit none

       integer iK,iZ,iA,m,id
       double precision eps0,PP,p
       dimension id(1)               ! Identities of up to mMax/1 ejectiles

       integer mxEin
       parameter (mxEin=200)         ! Max number of incident neutron energies

       logical, save :: dataread=.FALSE.
                                     ! static variable storing whether file
                                     ! 'ZA.xs' and 'ZA.PreEq' were read
       character(len=maxDATAPATH+100) dataf1, & ! name of data file 'ZA.xs'
                                      dataf2    ! name of data file 'ZA.PreEq'
       logical exists1, &            ! true if file 'ZA.xs' exists
               exists2               ! true if file 'ZA.PreEq' exists
       integer, save :: nPreEq=0     ! Number of files with pre-equilibrium
                                     ! neutron emission probabilities
                                     ! and spectra
       integer iiK                   ! index identifying the isotope
                                     ! considered (iiK=1,...,mxK);

       integer,save,dimension(:),allocatable::iKtoname
                                     ! mapping btw iK and index of file
                                     ! set 'ZA.xs' and 'ZA.PreEq'
       integer,save,dimension(:),allocatable::mEin
       integer,save,dimension(:,:),allocatable::maxE
       double precision,save,dimension(:,:),allocatable::En
       double precision,save,dimension(:,:),allocatable::xs
       double precision,save,dimension(:,:),allocatable::E
       double precision,save,dimension(:,:,:),allocatable::Y
       double precision,save,dimension(:,:,:),allocatable::F
       double precision,save,dimension(:),allocatable::dE
       double precision,save,dimension(:),allocatable::dEn

       dimension p(3,mMax) &    ! Momenta of up to mMax ejectiles
       ,PP(0:4)                 ! Momentum & energy of the nucleus (0:Mass)

       integer name             ! index of ZA in the list 'react.dat' of 
                                ! isotopes having pre-emission neutron data
       integer n,i
       double precision Ein,En0,dm,dp,P0,rng,D0,gsM,W0,E0,Df,Wf,eta &
       ,F1,F2,E1,Y1,E2,Y2,eps12,pk,d,Ef,Ek,epsk,phi,costh,sinth &
       ,px,py,pz,g0,VVX,VVy,VVz,V2,pdotV,g1,c,Eftot,msFREYA_SEPn
       ! double precision dErg

       character(len=70) headline ! Headline of data file
       double precision wN
       data wN/0.0/
       logical OFF
       OFF=.FALSE.

       if (OFF) RETURN          ! Pre-Eq is turned OFF

!========================================================================
       if (dataread.eqv..FALSE.) THEN 
                                ! INITIALIZE:   *********
         dataread=.TRUE.
         if (m.lt.0) then       !----------------------------------------
#ifdef WRITEL6
           write (L6,*) '** Turning Pre-Equilibrium emission OFF!'
#endif
           OFF=.TRUE.
           RETURN
         endif                  !----------------------------------------
#ifdef WRITEL6
         write (L6, & 
           "('PreEq: Preparing pre-equilibrium neutron emission:')")
#endif
         wN=wA+Dn               ! Neutron mass

!
! allocate memory for array mapping iK to name
!
         allocate(iKtoname(mxK))
!
         do iiK=1,mxK
           if (Mthk(iiK).gt.0) then
             dataf1=trim(freyadir)//"/"//trim(preprob(iiK))
             dataf2=trim(freyadir)//"/"//trim(prespec(iiK))
             inquire(FILE=dataf1, EXIST=exists1)
             inquire(FILE=dataf2, EXIST=exists2)
             if (exists1.eqv..TRUE..and.exists2.eqv..TRUE.) then
               nPreEq=nPreEq+1
             endif
           endif
         enddo
! now that we know how many pairs of 'ZA.xs' and 'ZA.PreEq' files 
! we need to read, i.e. nPreEq, allocate memory
         allocate(En(0:mxEin,nPreEq))
         allocate(xs(0:mxEin,nPreEq))
         allocate(maxE(0:mxEin,nPreEq))
         allocate(E(0:mxE,nPreEq))
         allocate(Y(0:mxE,0:mxEin,nPreEq))
         allocate(F(0:mxE,0:mxEin,nPreEq))
         allocate(mEin(nPreEq))
         allocate(dE(nPreEq))
         allocate(dEn(nPreEq))

         name=0
         do iiK=1,mxK
           exists1=.FALSE.
           exists2=.FALSE.
           if (Mthk(iiK).gt.0) then
             dataf1=trim(freyadir)//"/"//trim(preprob(iiK))
             dataf2=trim(freyadir)//"/"//trim(prespec(iiK))
             inquire(FILE=dataf1, EXIST=exists1)
             inquire(FILE=dataf2, EXIST=exists2)
           endif
           if ((exists1.eqv..FALSE.).and.(trim(preprob(iiK)).ne.'-')) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf1),' does not exist'
#endif
             write(error, 13) trim(dataf1)
13           format("FREYA: data file ", a, " could not be found. ", &
                    "Check react.dat file in FREYA data directory.")
             call seterrorcase(5,error)
             RETURN
           endif
           if ((exists2.eqv..FALSE.).and.(trim(prespec(iiK)).ne.'-')) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf2),' does not exist'
#endif
             write(error, 13) trim(dataf2)
             call seterrorcase(6,error)
             RETURN
           endif

           if ((exists1.eqv..FALSE.).or.(exists2.eqv..FALSE.)) cycle
           name=name+1
           iKtoname(iiK)=name
! prepare PreEq for this specie:
#ifdef WRITEL6
           write (L6,"('Preparing PreEq for nucleus',2i5)") & 
                  iZk(iiK),iAk(iiK)
#endif

!     READ PRE-EQUILIBRIUM PROBABILITIES:       *************************
! The PreEq probability table *.xs should have the following format:
!line 1 headline        up to 70 characters
!line 2 mEin, dEn       mEin incident-energy bins of width dEn
!1+mEin n, En, Pn       n=0(1)mEin, En=n*dEn, Pn: prob for PreEq emiss
!------------------------------------------------------------------------
#ifdef WRITEL6
           write (L6,*)'Reading PreEq "probabilities" data from file ', & 
              trim(preprob(iiK))
#endif
           OPEN (L7,file=trim(freyadir)//'/'//trim(preprob(iiK)), & 
                 status='old')
           read  (L7,'(a70)') headline
#ifdef WRITEL6
           write (L6,'("HEADLINE: ",a70)') headline
#endif
           read (L7,*) mEin(name),dEn(name)             ! Number of values Ein>0
           do n=0,mEin(name)
             read  (L7,*) i,En(n,name),xs(n,name)       ! #, Ein, Pn(Ein)
!            write (L6,"(2i5,f10.4,f10.6)") &
!            n,i,En(n,name),xs(n,name)
           enddo
           CLOSE (L7)
!       READ PRE-EQUILIBRIUM NEUTRON SPECTRA:   *************************
! The PreEq cross section table *.PreEq should have the following format:
!line 1 headline        up to 70 characters
!line 2 mEin, dEn       mEin incident-energy bins of width dEn (same as for xs)
!1+mEin sections:
!       each headed by  n, En, maxE(n), dE      (n=0(1)mEin and En=n*dEn);
!       followed by maxE(n) lines listing {E(i), Y(i,n)}
!       where E(i=i*dE, i=1(1)maxE(n), is the pre-eq neutron energy;
! Note:	All spectra have the SAME bin width dE common to all incident En!
!------------------------------------------------------------------------
#ifdef WRITEL6
           write (L6,*) 'Reading spectra PreEq(En) from file ', & 
                        trim(prespec(iiK))
#endif
           OPEN (L7,file=trim(freyadir)//'/'//trim(prespec(iiK)), & 
                 status='old')
           read  (L7,'(a70)') headline
#ifdef WRITEL6
           write (L6,'("HEADLINE: ",a70)') headline
#endif
           read  (L7,'(a70)') headline
#ifdef WRITEL6
           write (L6,'("HEADLINE: ",a70)') headline
#endif
!          It is required that mEin & dEn be the same as for .xs above!
#ifdef WRITEL6
           write (L6,"(i5,' incident energies En in steps of',f8.3)") & 
                 mEin(name),dEn(name)
#endif
           do n=0,mEin(name)
             do i=1,mxE
               Y(i,n,name)=0.0; F(i,n,name)=1.0
             enddo
             read  (L7,*) i,Ein,maxE(n,name),dE(name)   ! Must have i=n
!            write (L6,"(i5,f8.3,i5,f10.6)") i,Ein,maxE(n,name),dE(name)
             i=0; E(i,name)=0.0; Y(i,n,name)=0.0; F(i,n,name)=0.0
             do i=1,maxE(n,name)
               read (L7,*) E(i,name),Y(i,n,name)
               F(i,n,name)=F(i-1,n,name) & 
                  +(Y(i-1,n,name)+Y(i,n,name))*dE(name)/2  ! Integrate Y
             enddo
             maxE(n,name)=maxE(n,name)+1; i=maxE(n,name)
             if (i.eq.1) Y(0,n,name)=2.0/dE(name)
             F(i,n,name)=F(i-1,n,name) & 
                +(Y(i-1,n,name)+Y(i,n,name))*dE(name)/2
! OBS:  Erichs values of Y are normalized to unity, INT{Y(E)dE}=1,
!       EXCEPT for n=5 where the value is only half of what it should be!
!       the following renormalization takes care of this automatically:
             do i=0,maxE(n,name)
               Y(i,n,name)=Y(I,n,name)/F(maxE(n,name),n,name)
               F(i,n,name)=F(i,n,name)/F(maxE(n,name),n,name)
!               write(L6,*) i,E(i,name),Y(i,n,name),F(i,n,name)
             enddo
! For each incident energy En, the spectrum Y(i,n) is now tabulated for
! the energies E=0(dE)maxE(n)*dE with Y(0) and Y(maxE(n)) being zero;
! the associated distribution function (the incomplete integral of Y)
! F(i,n) then increases steadily from F(0,n)=0 to F(maxE(n),n)=1.0).
           enddo
           CLOSE (L7)

#ifdef WRITEL6
           write(L6,"('Done preparing PreEq for',3i5)") & 
                 iZk(iiK),iAk(iiK),mEin(name)
#endif
         enddo
       endif             ! INITIALIZE.  *********
!========================================================================
       ! The statement below will not work in general
       ! NAME=(iZ-90)/2   ! WARNING: fast way to assign NAME; it is NOT general!
       name=0
       name=iKtoname(iK)
       if (name.eq.0) write(*,*) 'Reminder to handle this, ', & 
                                 'which should never happen'

! Check if equivalent incoming neutron energy is above zero:
       En0=eps0-msFREYA_SEPn(iK,iZ,iA)          ! Incident neutron energy
       if (En0.lt.0.0) RETURN                   ! E* is too low for Pre-Eq
! Calculate the probability for pre-eq emission:
       ! This code does no work when dEn(name) is not constant, in which
       ! case we need to compute dE for each interval
       n=int((En0-En(0,name))/dEn(name))          ! => En(n) < En0 < En(n+1)
       ! More general scenario where dEn(name) is variable
       ! n=mEin(name)                             ! default value
       ! do i=0,mEin(name)
       !   if (En0.lt.En(i,name)) then
       !     n=i-1                                ! => En(n) <= En0 < En(n+1)
       !     exit
       !   endif
       ! enddo
       n=min(n,mEin(name)-1)                    ! => n=mEin if En0.GE.En(mEin)
       ! This code does no work when dEn(name) is not constant
       dm=(En0-En(n,name))/dEn(name)            ! Coefficients for the
       dp=1.0-dm
       ! More general scenario where dEn(name) is variable
       ! dErg=En(n+1,name)-En(n,name)
       ! dm=(En0-En(n,name))/dErg                 ! Coefficients for the
       ! dp=(En(n+1,name)-En0)/dErg               ! interpolation wrt En;
       P0=0.0
       if (xs(n,name).gt.0.0) then
         P0=xs(n,name)*dp+xs(n+1,name)*dm       ! interpolated probability
       endif
       if (rng(iseed).gt.P0) RETURN             ! No pre-equil emission
!       -----------------------------------------------------------------
! Carry out the emission of one pre-equilibrium neutron:
       D0=gsM(iZ,iA)                            ! Mother mass excess
       W0=iA*wA+D0                              ! Mother g.s. mass
       E0=W0+En0                                ! Mother total mass
       Df=gsM(iZ,iA-1)                          ! Daughter mass excess
       Wf=(iA-1)*wA+Df                          ! Daughter g.s. mass
       m=m+1                                    ! Emit one neutron
       id(m)=1                                  ! This ejectile is a neutron
       iA=iA-1                                  ! Mother loses one neutron
! Choose the kinetic energy from the normalized spectrum:  This is not
! corrected for the (small) recoil effect:  The energy value selected here is
! considered to be the combined kinetic energy, whereas the given differential
! cross sections are for the energy spectrum of the neutron only: small error.
       eta=rng(iseed); i=1; F2=0.0
       do while (eta.gt.F2)
         i=i+1
         F2=F(i,n,name)*dp+F(i,n+1,name)*dm     ! Interpol val of F(E2)
       enddo
       F1=F(i-1,n,name)*dp+F(i-1,n+1,name)*dm   ! Interpol val of F(E1)
!       We now have: F1=F(i-1) .lt. eta .le. F2=F(i) and E1.lt.eps12.le.E2:
       E1=E(i-1,name); Y1=Y(i-1,n,name)*dp+Y(i-1,n+1,name)*dm
       E2=E(i  ,name); Y2=Y(i  ,n,name)*dp+Y(i,  n+1,name)*dm
! We assume that the spectrum Y(eps12) is locally linear:
!       Y(eps12) = (Y1*(E2-eps12)+Y2*(eps12-E1))/dE for E1 < eps12 < E2;
! then F(eps12) is locally quadratic, so F(eps12):=eta can be solved to give
       if (Y1.eq.Y2) then
         eps12=E1+(eta-F1)/Y1/dE(name)
       else
         d=dE(name)
         eps12=E1+d*(sqrt(Y1**2+2*(eta-F1)*(Y2-Y1)/d)-Y1)/(Y2-Y1)
       endif
       eps0=En0-eps12                           ! Excitation in daughter
       Ef=Wf+eps0                               ! Total mass of daughter
! Calculate the magnitude of the ejectile momentum pk in CM:
!      pk=0.5*sqrt((E0-Ef-wN)*(E0+Ef+wN)*(E0-Ef+wN)*(E0+Ef-wN))/E0 -> 
       pk=0.5*sqrt(eps12*(eps12+2*wN)*(E0+Ef+wN)*(E0+Ef-wN))/E0
       Ek=sqrt(wN**2+pk**2)                     ! Ejectile total energy
       epsk=Ek-wN                               ! Kinetic energy of ejectile
! Choose the direction of the ejectile motion:
       phi=twopi*rng(iseed)                     ! Azimuthal angle phi
       costh=2.0*rng(iseed)-1.0                 ! cos(polar angle theta)
       sinth=sqrt(1.0-costh**2)                 ! sin(polar angle theta)
       px=pk*sinth*cos(phi)                     !>
       py=pk*sinth*sin(phi)                     ! > Ejectile momentum
       pz=pk*costh                              !>
! Boost the momenta to the overall reference frame:
       g0 =PP(4)/PP(0)                          ! Lorentz boost coefficient
       VVx=PP(1)/PP(4)                          !>
       VVy=PP(2)/PP(4)                          ! > Velocity of mother nucleus
       VVz=PP(3)/PP(4)                          !>
       V2=VVx**2+VVy**2+VVz**2                  ! V * V
       pdotV = px*VVx + py*VVy + pz*VVz         ! p * V
       g1=0.0; if (V2.gt.0.0) g1=(g0-1.0)*pdotV/V2
! Boost the ejectile:
       c = g0*Ek + g1                           ! Boost mass
       p(1,m) = c*VVx + px                      !>
       p(2,m) = c*VVy + py                      ! > Momentum of ejectile
       p(3,m) = c*VVz + pz                      !>
! Boost the daughter:
       PP(0)=Ef                                 ! Total rest energy of daughter
       Eftot=sqrt(Ef**2+pk**2)                  ! Daughter total energy in CM
       Eftot = g0*Eftot - g1                    ! Boost its total energy
       PP(1) = Eftot*VVx - px                   !>
       PP(2) = Eftot*VVy - py                   ! > Daughter momentum
       PP(3) = Eftot*VVz - pz                   !>
       PP(4)=sqrt(PP(0)**2+PP(1)**2+PP(2)**2+PP(3)**2) ! 4: Energy of daughter
       RETURN
!========================================================================
       END

!************************************************************************
       SUBROUTINE EVAPS1n (kEv,iK,iZ,iA,eps,SS,PP,m,p)
!  simulates ONE single neutron evaporation from a compound nucleus whose
!  excitation energy eps is above the neutron separation energy Sn so that
!  evaporation is energetically allowed.
!  Although the treatment is similar to what is being done for the fragments,
!  it is easier to have a separate routine for the pre-fission evaporation,
!  since certain quantities (e.g. level densities) may differ.
!  The energetics is done differently: since Sn has been specified (and it
!  is not identical to the calculated value), we calculate the daughter M
!  by subtracting mN+Q* from the mother M*; the resulting M* of the fissioning
!  nucleus is then carried forward consistently (in PP(0)).
!
! INPUT:
! kEv           Event number
! iZ,iA,eps     Charge number, mass number, excitation of the mother nucleus
! SS(1:3)       Fragment spin components of the mother nucleus (in hbar)
! PP(0:4)       Four-momentum of the mother nucleus (0:mass,1-3:P,4:energy)
! m.GE.0        The number of previously emitted neutrons in this decay chain
! OUTPUT:
! iZ,iA,eps     Charge number, mass number, excitation of the daughter nucleus
! PP(0:4)       Four-momentum of daughter nucleus (0:mass,1-3:mom,4:energy) 
! SS(1:3)       Fragment spin components of the daughter nucleus (in hbar)
! m > 0         m(out)=m(in)+1: the accumulated number of emitted neutrons
! p(1:3,1:m)    The cartesian momentum components of this evaporated neutron
!
       use freyaout
       use freyacmass
       use freyaconsts
       use freyaseed
       use freyakmass
       use freyaparameters
       use freyacLOOK
       use freyaQmin
       use freyacROT
       use freyaparams

       implicit none

       integer kEv,iK,iZ,iA,m,iAf
       double precision eps,SS,PP,p

       dimension & 
        p(3,mMax) &                ! Momenta of up to mMax ejectiles
       ,PP(0:4) &                  ! Four-momentum of the nucleus (0:M,4:E)
       ,SS(1:3)                    ! Fragment Spin Components

       double precision wn,Ei,Wi,S2,epsi,Qi,Yi,Sn,Qn,Wf,EiRot,EfRot,Yf, &
            aA,U,T,Ekin,eta1,eta2,Ef,pn,RA,phi, &  
            costh,sinth,x,y,z,ax,ay,az,bx,by,bz,cx,cy,cz, &
            cosphi,sinphi, &
            eta,a,b,c,px,py,pz,p2,ox,oy,oz,wx,wy,wz, &
            En,epsn,SSx,SSy,SSz,Sf2,Efkin,Qf,Di,Df, &
            g0,VVx,VVy,VVz,V2,pdotV,g1,cc,Eftot, &
            gsm 
       double precision alevel,rng,log,rra
       double precision egs
       data egs/0.0/            ! Deformation of the daughter nucleus

!  100  format ('EVAP1n m:',i6,i4,f8.3,f12.3,3f8.3,2f8.3)!!
!  101  format ('EVAP1n d:',i6,i4,f8.3,f12.3,3f8.3,2f8.3)!!
!  102  format ('Check:',3f10.6,f15.6)

!       do i=1,4
!         q(i)=PP(i)
!       enddo

       wn=wA+Dn                                ! Mass of ejectile #1: neutron
       Ei=PP(0)                                ! Total mother rest mass
       Wi=iA*wA+gsM(iZ,iA)                     ! Mother gs mass
       S2=SS(1)**2+SS(2)**2+SS(3)**2           ! Spin squared of mother
       EiRot=0.5*S2/ROT(iA)                    ! Rotational energy of mother
       epsi=eps        ! not needed            ! Mother excitation: rot+heat
       Qi=epsi-EiRot                           ! Mother heat: Q
       Yi=Wi+EiRot        ! not needed!        ! Mother yrast mass: gs+rot
       iAf=iA-1                                ! Daughter mass number
       Sn=gsM(iZ,iAf)+Dn-gsM(iZ,iA)            ! Neutron separation energy
       Qn=eps-Sn        ! or Qi-Sn?            ! Q* value for n emission (S=0)
       Wf=iAf*wA+gsM(iZ,iAf)                   ! Daughther gs mass
       EfRot=0.5*S2/ROT(iAf)                   ! Rot erg of daughter if SS=c
       Yf=Wf+EfRot        ! not needed!        ! Daughter yrast mass: gs+rot
       Qn=Qn+EiRot-EfRot                       ! Adjusted Q value (for S>0)
!        Neutron emission occurs if the Q value exceeds the specified Qmin.
       if (Qn.lt.Qmin) goto 19                 ! This evap is too soft -> quit

       aA=alevel(iK,iZ,iAf,Qn,U)               !!! new argument list, U is output
       T=sqrt(U/aA)                            ! Max temperature in daughter
   12  Ekin=Qn                                 ! Select energy release:
!       write (6,"(i8,'#',i3,f10.5)") kEv,m,Qn!; pause 'EvapS1n'
       do while (Ekin.lt.Emin.OR.Ekin.ge.Qn) 
         eta1=-log(rng(iseed))
         eta2=-log(rng(iseed))
         Ekin=(eta1+eta2)*T                    ! Tot kin energy (2-body CM)
       enddo
       eps=Qn-Ekin                             ! Excitation in daughter
       Ef=Yf+eps                               ! Daughter rst erg: gs+rot+heat
! Calculate the magnitude of the ejectile momentum pn:
!       pn=0.5*sqrt((Ei-Ef-w)*(Ei+Ef+w)*(Ei-Ef+w)*(Ei+Ef-w))/Ei -> 
       pn=0.5*sqrt(Ekin*(Ekin+2*wn)*(Ei+Ef+wn)*(Ei+Ef-wn))/Ei
! Choose the emission point r=(x,y,z) randomly on the spherical nuclear surface:
       RA=RRA(iA)                              ! Nuclear radius
       phi=twopi*rng(iseed)                    ! Azimuthal angle phi
       cosphi=cos(phi)                         ! cos(azimuthal angle phi)
       sinphi=sin(phi)                         ! sin(azimuthal angle phi)
       costh=2.0*rng(iseed)-1.0                ! cos(polar angle theta)
       sinth=sqrt(1.-costh**2)                 ! sin(polar angle theta)
       x=RA*sinth*cosphi                       !>
       y=RA*sinth*sinphi                       ! > r = (x,y,z)
       z=RA*costh                              !>
! Construct the local reference system abc:
       ax= cosphi*costh; ay=sinphi*costh; az=-sinth  ! a=(ax,ay,az)
       bx=-sinphi;       by=cosphi;       bz= 0.0    ! b=(bx,by,bz)
       cx= cosphi*sinth; cy=sinphi*sinth; cz= costh  ! c is normal to surf
!        write (6,"('x:',4f10.6)")  x, y, z,sqrt( x**2+ y**2+ z**2)
!        write (6,"('o:',4f10.6)") ox,oy,oz,sqrt(ox**2+oy**2+oz**2)
!        write (6,"('w:',4f10.6)") wx,wy,wz,sqrt(wx**2+wy**2+wz**2)
! Choose the direction u of the ejectile motion (biased with u*c):
       phi=twopi*rng(iseed)                          ! Azimuthal angle phi (0,2pi]
       eta=rng(iseed)                                ! cos(polar angle theta)**2
       costh=sqrt(eta)                               ! cos(polar angle theta) (0,1]
       sinth=sqrt(1.0-eta)                           ! sin(polar angle theta) [0,1)
       a=sinth*cos(phi); b=sinth*sin(phi); c=costh   ! dir of p in abc syst
! Components of emitted neutron in the local comoving system:
       px=pn*(a*ax+b*bx+c*cx)                        !>  Ejectile momentum
       py=pn*(a*ay+b*by+c*cy)                        ! > p=(px,py,pz) [pc: energy]
       pz=pn*(a*az+b*bz+c*cz)                        !>  in comoving frame
       p2=px**2+py**2+pz**2
! Calculate the rotation vector omega/c = (ox,oy,oz) [omega ~ S for a sphere]:
       ox=SS(1)/ROT(iA)/hbarc                        !>  This is omega/c
       oy=SS(2)/ROT(iA)/hbarc                        ! > given in inverse fermi:
       oz=SS(3)/ROT(iA)/hbarc                        !>  [omega/c]=(1/T)/(L/T)=1/L
! Collective velocity w/c=(wx,wy,wz) at the emission point r=(x,y,z): w = o x r
       wx=oy*z-oz*y; wy=oz*x-ox*z; wz=ox*y-oy*x!   [w/c]=1
! Components of the emitted neutron after rotational boost:
       px = px + wn*wx                               !>  Ejectile momentum
       py = py + wn*wy                               ! > p=(px,py,pz) [pc: energy]
       pz = pz + wn*wz                               !>  after rotational boost
       p2 = px**2+py**2+pz**2                        ! Neutron momentum squared
       En=sqrt(wn**2+p2)                             ! Total neutron energy in CM
       epsn=En-wn                                    ! Neutron kinetic energy in CM
! Calculate tentative spin of daughter nucleus:
       SSx=SS(1)-(y*pz-z*py)/hbarc                   !>
       SSy=SS(2)-(z*px-x*pz)/hbarc                   ! >  S' = S - r x p
       SSz=SS(3)-(x*py-y*px)/hbarc                   !>
       Sf2=SSx**2+SSy**2+SSz**2                      ! (daughter spin)**2

! Check energy bound (non-relativistically):
!       write (6,"('S:',4f10.4)") SSx,SSy,SSz,sqrt(Sf2)
       EfRot=0.5*Sf2/ROT(iAf)                        ! Tentative daughter rot energy
!        write (6,*) EfRot,sqrt(S2),ROT(iAf),iAf
       Yf=Wf+EfRot                                ! Tent. yrast mass of daughter
!        write (6,*) Yf,iAf*wA+gsM(iZ,iAf),EfRot
       Efkin=(p2/wn+p2/Wf)/2                        ! Tentative final kin energy
!        Qf=EiRot+epsi-Sn-EfRot-Efkin                ! same as:        
       Qf=Qn+(EiRot-EfRot)-Efkin                ! Tentative heat in daughter
#ifdef WRITEL6
       if (Qf.lt.0.0) then ! this emission violates energy conservation:
                write (L6,"(i8,' EvapS1n: Q<0 - trying again ..')") kEv
                Di=gsM (iZ,iA)                        ! Mother mass excess
                Df=gsM (iZ,iAf)                        ! Daughter mass excess
                write (L6,"(2i8,4f10.4)") kEv,iAf,Qf,epsn,EiRot,EfRot
                write (L6,"(2f12.4)") PP(0)
                write (L6,"(2f12.4)") Ei
                write (L6,"(2f12.4)") iA*wA+Di,        iAf*wA+Df
                write (L6,"(2f12.4)") Ei-epsi,        Ef-eps
                write (L6,"(2f12.4)") epsi,        Qf
                write (L6,"(2f12.4)") EiRot,        EfRot
                write (L6,"(2f12.4)") Yi,        Yf
                write (L6,"(2f12.4)") Ei,        Ef
       endif
#endif
       if (Qf.le.0.0) goto 12                        ! try again!
       iA=iAf                                        ! Mass number of daughter
       eps=Qf                                        ! Actual heat in daughter
       Ef=Yf+eps                                ! Actual yrast mass of daughter
       SS(1)=SSx; SS(2)=SSy; SS(3)=SSz; S2=Sf2        ! Actual daughter spin

! Boost the momenta to the overall reference frame:
       g0 =PP(4)/PP(0)                 ! Lorentz boost coefficient
       VVx=PP(1)/PP(4)                 !>
       VVy=PP(2)/PP(4)                 ! > Velocity of mother nucleus
       VVz=PP(3)/PP(4)                 !>
       V2=VVx**2+VVy**2+VVz**2         ! V * V
       pdotV = px*VVx + py*VVy + pz*VVz! p * V
       g1=0.0; if (V2.gt.0.0) g1=(g0-1.0)*pdotV/V2
! Boost the ejectile:
       m=m+1                           ! One more neutron emitted
       cc = g0*En + g1                  ! Boost mass
       p(1,m) = cc*VVx + px             !>
       p(2,m) = cc*VVy + py             ! > Momentum of ejectile #m
       p(3,m) = cc*VVz + pz             !>
! Boost the daughter:
       PP(0)=Ef                        ! Total daughter rest mass
       Eftot=sqrt(Ef**2+pn**2)         ! Daughter total energy in CM
       cc = g0*Eftot - g1              ! convienient coefficient
       PP(1) = cc*VVx - px             !>
       PP(2) = cc*VVy - py             ! > Momentum of daughter
       PP(3) = cc*VVz - pz             !>
       PP(4) = g0*(Eftot - pdotV)      ! Total daughter energy
!        PP(4)=sqrt(PP(0)**2+PP(1)**2+PP(2)**2+PP(3)**2)        ! Total daughter energy
!        write (6,"(i4,4x,i4,6f10.4)") kEv,iA,eps,(SS(i),i=1,3)!!
   19   RETURN
       END

!************************************************************************
       SUBROUTINE DDECAY (iK,kEv,iZ,iA,eps,PP,mmk,m,id,p)
! called as       DDECAY (k,iZ1,iA2,eps1,PP1,mult1,m1,id1,p1) for fragm #1;
! called as       DDECAY (k,iZ2,iA2,eps2,PP2,mult2,m2,id2,p2) for fragm #2.
! simulates sequential decay of an excited fission fragment:
! INPUT:
! iK:           index identifying the isotope considered (iK=1,...,mxK);
! kEv  :        Event number (for monitoring only).
! iZ,iA,eps:    Charge number, mass number, excitation of mother nucleus.
! PP(0:4):      Rest mass (0), mom (1-3) and tot erg (4) of mother nucleus;
! mmk(0:Nk)      Max allowed multiplicity of each ejectile type k=0:Nk:
!   mmk(k).GE.0:   at most that many type-k emissions will be admitted;
!   mmk(k).LT.0:   no restrictions on the number of type-k emissions.
! m.GE.0        Total number of ejectiles so far (i.e. prior to this call)
! Note:         To carry out at most one neutron emission,
!               specify mm(1)=1 and mm(k.NE.1)=0 at entry.
! OUTPUT:
! iZ,iA,eps:    Charge number, mass number, excitation of daughter nucleus,
! OBS:          eps is now the post-neutron excitation (BEFORE photon emiss)
!               which is more informative since E*=0 after the photon cascade.
! PP(0:4):      Rest mass (0), mom (1-3) and tot erg (4) of daughter nucleus;
! maxk(0:Nk)    The actual multiplicity of each ejectile type;
! m.GE.0        total number of ejectiles now: m = mm(0) +...+ mm(Nk);
! id(1:m):      their type (id=0:Nk) [only id=0(g) & id=1(n) here!];
! p(1:4,1:m):   their cartesian momentum components (1-3) & tot erg (4).
! -----------------------------------------------------------------------
!OBS:  The four-vector convention follows standard FREYA:
!        0              rest mass
!       1-3             momentum
!        4              total energy
! -----------------------------------------------------------------------
!! This version considers neutron emission FOLLOWED by photon emission;	!!
!! generally the different decay modes should compete at each stage,	!!
!! so the present sequential scheme is only preliminary and schematic.	!!

       use freyaout
       use freyaconsts
       use freyaseed
       use freyakmass
       use freyaparameters
       use freyacLOOK
       use freyacmass
       use freyaQmin
       use freyaparams

       implicit none

       integer iK,kEv,iZ,iA,mmk,m,id
       double precision eps,PP,p

       integer kZ,kA
       double precision Qm,epsm
       dimension id(mMax) &     ! Identity of up to mMax ejectiles
       ,mmk(0:Nk) &             ! Multiplicity of each ejectile type
       ,p(4,mMax) &             ! Momenta of up to mMax ejectiles
       ,PP(0:4) &               ! Four-momentum of the nucleus (0:M,4:E)
       ,kZ(0:Nk) &              ! Charge number of k-type ejectiles
       ,kA(0:Nk)                ! Mass number of k-type ejectiles
!     c	,Gamma(0:Nk)		! Partial width for each mode k
       dimension Qm(mMax),epsm(mMax) ! to print if multiplicity exceeds mMax

! Data for Nk=1:
       data kZ /0,0/            ! Charge number of various ejectile types
       data kA /0,1/            ! Mass number of various ejectile types

       integer iA0,max,n,k,iZf,iAf,Nsoft
       double precision epsi0,Di,gsM,Wi,Ei,w,Zf,epsi,Af,Df,Wf,Qk,aA,T,U,sqrt & 
       ,eps12,eta1,rng,log,eta2,epsf,Ef,pk,En,epsk,phi,costh & 
       ,sinth,px,cos,py,sin,pz,g0,VVx,VVy,VVz,alevel,V2,pdotV & 
       ,g1,cc,Eftot,eta3,Ek,Esoft
       double precision epsmin ! Maximum heat left after statistial photons.

#ifdef WRITEL6
       if (LOOK) write (L6,*) kEv,': DECAY iseed =',iseed,eps
#endif
       iA0=iA; epsi0=eps                        ! Initial A & E*
       epsi=eps                                 ! Mother excitation
       Di=gsM (iZ,iA)                           ! Mother mass excess
       Wi=iA*wA+Di                              ! Mother ground-state mass
! Choose among the competing decay channels k = gamma, n, p, d, t, tau, alpha:
!       eta=rng(iseed)                          ! to select ejectile type k = 0,...,Nk
!      Set epsmin
       epsmin=epsmink(iK)
! Consider ONLY neutrons (followed by photons) in this version:

! NEUTRON emission:     -------------------------------------------------
       k=1                                      ! k=1: Ejectile is a neutron
       max=mmk(k)                               ! max # k-emissions in total
       n=0                                      ! -> actual k multiplicity
   10  if (n.eq.max) goto 19                    ! Max # k-emissions is reached
       Ei=Wi+epsi                               ! Rest energy of mother nucleus
       w=wk(k)                                  ! Ejectile mass
       iZf=iZ-kZ(k); Zf=dble(iZf)              ! Daughter charge number
       iAf=iA-kA(k); Af=dble(iAf)              ! Daughter mass                                
       Df=gsM (iZf,iAf)                         ! Daughter mass excess
       Wf=iAf*wA+Df                             ! Daughter ground-state mass
       Qk=Ei-Wf-w                               ! Q value for this mode
!      Qk=Di+epsi-Df-Dk(k)                      ! Q value (~200x more accurate)
       if (Di+Df.ge.500.) then
#ifdef WRITEL6
         write (L6,*) iZ,iA,' -> ',iZf,iAf,Di,Df,epsi,m
#endif
         ! pause!!
       endif
! Check the Q-value (for neutron emission):
       if (Qk.le.Qmin) goto 19                  ! Such a decay is impossible

! Choose the tentative kinetic energy of relative motion (ie tot kin erg in CM):
       n=n+1                                    ! Increment k multiplicity
       m=m+1                                    ! Increment total multiplicity
       if (m.gt.mMax) then
#ifdef WRITEL6
         write (L6,*) 'Emitting neutrons in DECAY:'
         write (L6,*) 'Ejectile multiplicity exceeds specified max =', & 
           mMax,'!'
         write (L6,*) 'Suggestion: Increase the parameter mMax!'
         write (L6,*) 'Event',kEv,':',iZ,iA-iZ,iA,iA0,epsi0
         do m=1,mMax
           write (L6,*) 'm =',m,':',wk(id(m)),Qm(m),epsm(m)
         enddo
#endif
         STOP!!
       endif
       id(m)=k                                  ! Store identity of ejectile
! Calculate the level-density parameter:
       aA=alevel(iK,iZf,iAf,Qk,U)               ! alevel(iZ,iA,eps,U)

       T=sqrt(U/aA)                             ! Max temp in daughter
       eps12=Qk
       do while (eps12.lt.Qmin.OR.eps12.ge.Qk) 
         eta1=-log(rng(iseed))
         eta2=-log(rng(iseed))
         eps12=(eta1+eta2)*T                    ! Tot kin energy (in 2-body CM)
       enddo
       Qm(m)=Qk; epsm(m)=eps12                  ! to print if mult > mMax
       epsf=Qk-eps12                            ! Excitation in daughter
       Ef=Wf+epsf                               ! Total rest energy of daughter
! Calculate the magnitude of the ejectile momentum pk:
!       pk=0.5*sqrt((Ei-Ef-w)*(Ei+Ef+w)*(Ei-Ef+w)*(Ei+Ef-w))/Ei -> 
       pk=0.5*sqrt(eps12*(eps12+2*w)*(Ei+Ef+w)*(Ei+Ef-w))/Ei
       En=sqrt(w**2+pk**2)                      ! Neutron total energy in CM
       epsk=En-w                                ! Kinetic energy of ejectile
! Choose the direction of the ejectile motion:
       phi=twopi*rng(iseed)                     ! Azimuthal angle phi
       costh=2.*rng(iseed)-1.                   ! cos(polar angle theta)
       sinth=sqrt(1.-costh**2)                  ! sin(polar angle theta)
       px=pk*sinth*cos(phi)                     !>
       py=pk*sinth*sin(phi)                     ! > Ejectile momentum
       pz=pk*costh                              !>
! Boost the momenta to the overall reference frame:
       g0 =PP(4)/PP(0)                          ! Lorentz boost coefficient
       VVx=PP(1)/PP(4)                          !
       VVy=PP(2)/PP(4)                          ! > Velocity of mother nucleus
       VVz=PP(3)/PP(4)                          !
       V2=VVx**2+VVy**2+VVz**2                  ! V * V
       pdotV = px*VVx + py*VVy + pz*VVz         ! p * V
       g1=0.0; if (V2.gt.0.0) g1=(g0-1.0)*pdotV/V2
! Boost the ejectile:
       p(4,m)=g0*(En+pdotV)                     ! Total energy of ejectile
       cc=g0*En+g1                              ! cc: convenient coefficient
       p(1,m) = cc*VVx + px                     !>
       p(2,m) = cc*VVy + py                     ! > Momentum of ejectile
       p(3,m) = cc*VVz + pz                     !>
! Boost the daughter:
       PP(0)=Ef                                 ! 0: Total daughter rest energy
       Eftot=sqrt(Ef**2+pk**2)                  ! Total daughter energy in CM
       PP(4) = g0*(Eftot-pdotV)                 ! 4: Total daughter energy
       cc=g0*Eftot-g1                           ! cc: convenient coefficient
       PP(1) = cc*VVx - px                      !>
       PP(2) = cc*VVy - py                      ! > 1-3: Daughter momentum
       PP(3) = cc*VVz - pz                      !>
! Daughter becomes next mother:
       iZ=iZf
       iA=iAf
       Di=Df
       Wi=Wf                                    ! g-s mass of new mother
       epsi=epsf                                ! Excitation of new mother
#ifdef WRITEL6
       if (LOOK) write (L6,*) n,': DECAY iseed =',iseed,epsi
#endif
       goto 10                                  ! Repeat decay procedure
! Complete neutron decay chain:
   19  mmk(k)=n                                 ! k-type multiplicity
       eps=epsi      ! -> OUTPUT                ! Excitation after n evap
       Ei=Wi+epsi                               ! Total residue rest energy

! PHOTON emission:      -------------------------------------------------
! Method:
! At each step, calculate the nuclear temperature from E* (= eps_max)
! and then select the photon energy eps and reduce E* -> E* - eps.
! If eps < eps_min then do not register this photon but iterate the procedure;
! Stop when eps_max < eps_min.
       Nsoft=0; Esoft=0.0                       ! Photons below eps_min
       k=0                                      ! k=0: Ejectile is a photon
       max=mmk(k)                               ! Specified max # k-emissions
       n=0                                      ! -> actual k multiplicity
   90  if (n.eq.max) goto 99                    ! Max # k-emissions is reached
       Qk=epsi                                  ! Q-value for photon emission
       if (Qk.lt.epsmin) then 
         goto 99                                ! Only soft photons left: stop
       else                                     ! select the photon energy:
         aA=alevel(iK,iZ,iA,Qk,U)
         T=sqrt(U/aA)! T before
         epsk=Qk                                ! Max kin erg energy release
         do while (epsk.ge.Qk)                  ! We must demand epsk < Qk 
           eta1=-log(rng(iseed))
           eta2=-log(rng(iseed))
           eta3=-log(rng(iseed))
           epsk=(eta1+eta2+eta3)*T              ! Kinetic energy (in 2-body CM)
         enddo
       endif
! Count:
       n=n+1                                    ! Increment k multiplicity
       m=m+1                                    ! Increment total multiplicity
       if (m.gt.mMax) then
#ifdef WRITEL6
         write (L6,*) 'Emitting photons in DECAY:'
         write (L6,*) 'Ejectile multiplicity exceeds specified max =', & 
           mMax,'!'
         write (L6,*) 'Suggestion: Increase the parameter mMax!'
         write (L6,*) 'Event',kEv,':',iZ,iA-iZ,iA,iA0,epsi0
         do m=1,mMax
           write (L6,*) 'm =',m,':',wk(id(m)),Qm(m),epsm(m)
         enddo
#endif
         STOP!!
       endif
       id(m)=k                                  ! Store identity of ejectile
       Qm(m)=Qk; epsm(m)=epsk                   ! to print if mult > mMax
       epsf=Qk-epsk                             ! Excitation in daughter
       Ef=Wi+epsf                               ! Total mass of daughter
! Calculate the magnitude of the ejectile momentum pk:
       pk=0.5*epsk*(Ei+Ef)/Ei                   ! Recall: mass(photon) = 0
       Ek=pk                                    ! Ejectile total energy
       epsk=Ek                                  ! Kinetic energy of ejectile
! Choose the direction of the ejectile motion:
       phi=twopi*rng(iseed)                     ! Azimuthal angle phi
       costh=2.*rng(iseed)-1.                   ! cos(polar angle theta)
       sinth=sqrt(1.-costh**2)                  ! sin(polar angle theta)
       px=pk*sinth*cos(phi)                     !>
       py=pk*sinth*sin(phi)                     ! > Ejectile momentum
       pz=pk*costh                              !>
! Boost the momenta to the overall reference frame:
       g0 =PP(4)/PP(0)                          ! Lorentz boost coefficient
       VVx=PP(1)/PP(4)                          !>
       VVy=PP(2)/PP(4)                          ! > Velocity of mother nucleus
       VVz=PP(3)/PP(4)                          !>
       V2=VVx**2+VVy**2+VVz**2                  ! V * V
       pdotV = px*VVx + py*VVy + pz*VVz         ! p * V
       g1=0.0; if (V2.gt.0.0) g1=(g0-1.0)*pdotV/V2
! Boost the ejectile:
       p(4,m)=g0*(Ek+pdotV)
       cc=g0*Ek+g1                              ! cc: convenient coefficient
       p(1,m) = cc*VVx + px                     !>
       p(2,m) = cc*VVy + py                     ! > Momentum of ejectile
       p(3,m) = cc*VVz + pz                     !>
! Boost the daughter:
       PP(0)=Ef                                 ! 0: Total daughter rest energy
       Eftot=sqrt(Ef**2+pk**2)                  ! Total daughter energy in CM
       PP(4) = g0*(Eftot-pdotV)                 ! 4: Total daughter energy
       cc=g0*Eftot-g1                           ! cc: convenient coefficient
       PP(1) = cc*VVx - px                      !>
       PP(2) = cc*VVy - py                      ! > 1-3: Daughter momentum
       PP(3) = cc*VVz - pz                      !>
! Check photon energy: Do not register it if it is too soft: 
       if (epsk.lt.epsmin) then                 ! Too soft to register:
         Nsoft=Nsoft+1                          ! Number of soft photons
         Esoft=Esoft+epsk                       ! Energy of soft photons
         n=n-1; m=m-1                           ! Reset ejectile multiplicity
       endif
! Daughter becomes next mother:
       epsi=epsf; Ei=Ef
#ifdef WRITEL6
       if (LOOK) write (L6,*) m-n,': GAMMA iseed =',iseed,epsi
#endif
       goto 90                                  ! Repeat decay procedure
! Complete photon decay chain:
   99  mmk(k)=n                                 ! k-type multiplicity
!       write (L6,"('DECAY:',i4,' +',i4,' =',i4)") mmk(1),mmk(0),m!!
       RETURN
       END

!************************************************************************
