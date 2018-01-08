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
        SUBROUTINE DecayS (iK,phase,iZ,iA,eps,SS,PP,mult,m,id,p)
! called as 	   DecayS (k,iZ,iA,eps1,SS1,PP1,mult1,m1,id1,p1) for fragm #1
! simulates sequential decay of an excited fission fragment:
! ===============================================
!phase=0: Initialize residual temperatures and specify gmin via eps:
! iZ,iA:	Charge & mass numbers of the combined system
! eps:		specifies gmin, the minimum energy for photon recording.
! m:	m>0:	specifies mp, the power of eta in the photon spectrum [std:2]
! 	m<0:	the usual spectrum is modified by the factor eps**|m|
!		and sampled by biasing with (eps/epsmax)**|m|.
! ===============================================
!phase>0: Carry out one decay chain:
! INPUT:      	---------------------------------
! iK:           index identifying the isotope considered (iK=1,...,mxK);
! phase:	0: initialization
!               >0: carry out one decay chain
!               <0: analysis
! iZ,iA,eps:	Charge & mass numbers and HEAT (stat. exc.) of mother nucleus
! SS(1:3):	Initial fragment spin components (in units of hbar)
! PP(0:4):	Four-momentum of the mother nucleus (0:M, 1-3:P, 4:E)
! mult(0:Nk)	Max allowed multiplicity of each ejectile type k=0:Nk:
!  mult(k).GE.0: at most that many type-k emissions will be admitted;
!  mult(k).LT.0: no restrictions on the number of type-k emissions.
! m.GE.0	Total number of ejectiles so far (i.e. prior to this call)
! Note:		To carry out at most one neutron emission,
!		specify mult(1)=1 and mult(k.NE.1)=0 at entry.
! OUTPUT:	---------------------------------
! iZ,iA,eps:	Charge & mass numbers and post-evap excitation of daughter.
! OBS:		eps is now the post-neutron excitation (BEFORE photon emiss)
!		which is more informative since E*=0 after the photon cascade;
!		eps contains the rotational erg, just as the input eps does.
! SS:		Final fragment spin components (should be zero)
! PP(0:4):	Four-momentum of the daughter nucleus (0:M, 1-3:P, 4:E).
! mult(0:Nk)	The actual multiplicity of each ejectile type;
! m.GE.0	Total number of ejectiles now: m = m + mult(0) +...+ mult(Nk);
! id(1:m):	their type (id=0:Nk) [only id=0(g) & id=1(n) here!];
! p(1:3,1:m):	their cartesian 3-momentum components.
! ===============================================
!phase<0: Finalize residual temperatures.

!! This version considers neutron emission FOLLOWED by photon emission;	!!
!! generally the different decay modes should compete at each stage,	!!
!! so the present sequential scheme is only preliminary and schematic.	!!
!
!  After neutron emission has ceased, statistical photons are first emitted
!  until the yrast line is reached (to within epsmin) and then stretched E2
!  protons are emitted, each one causing the spin to decrease by two units;
!  the final excitation is then radiated in the form of one last photon.
       use freyaparameters
       use freyaout
       use freyaconsts
       use freyaseed
       use freyacmass
       use freyakmass
       use freyaQmin
       use freyacROT
       use freyaSource
       use freyaRIPL
       use freyacount
       use freyaTdist
       use freyalocal
       use freyaaveE
       use freyaDecayS
       use freyaPath
       use freyaerrors
       use freyaparams
       use freyaEv
       use freyaGD

       implicit none

       integer iK,phase,iZ,iA,mult,m,id
       double precision eps,SS,PP,p

       dimension id(mMax) &    ! Identity of up to mMax ejectiles
       ,mult(0:Nk) &           ! Multiplicity of each ejectile type
       ,p(3,mMax) &            ! Momenta of up to mMax ejectiles
       ,PP(0:4) &              ! Four-momentum of the nucleus (0:M,4:E)
       ,kZ(0:Nk) &             ! Charge number of k-type ejectiles
       ,kA(0:Nk) &             ! Mass number of k-type ejectiles
       ,SS(3)                  ! Fragment angular momentum
       dimension Qm(mMax),epsm(mMax)   ! to print if multiplicity exceeds mMax
       double precision ox,oy,oz,p2,pdotv,phi,pk,px,py,pz,qf,qk,qm,&
            ra,riple,ripln,s0,s2,sf2,sinit,sinphi,sinth,&
            sn,softe,softn,spost,ssx,ssy,ssz,state,&
            statn,sx,sy,sz,t,tmax,u,v2,vvx,vvy,vvz,w,wf,wi,wx,wy,&
            wz,y,yf,yi,z,zf,gsm
       integer i,iaf,iCOLL,iCOUNT,iLAST,index,iripl,iSOFT,iSTAT,isum,&
               iZf,j,k,k0,ka,kinc,kk,koll,kZ,l,l0,linc,ll,&
            max,missing,mnZ,mxz,n,n0,ninc,nn,ntry
       double precision cc,d,df,di,ds,ef,efkin,&
            EfRot,Eftot,ei,EiRot,ek
       double precision En,eps0,eps12,epsf,epsi,epsi0,epsk,epsm,&
            epstot,eta,eta1,eta2,ev,&
            g0,g1,GDRmax,gk,gk1,gk2,gk3
       data Sinit,Spost/0.0,0.0/       ! Spin magnitudes before & after evap
       data statN,statE/0.0,0.0/       ! Number & energy of stat photons
!x	data tranN,tranE/0.0,0.0/	       ! Number & energy of tran photons
       data missing/0/                 ! Nuclei not in RIPL table
       data dS/1.0/                    ! Reduction of S in statistical decays
       data iSOFT,iSTAT,iRIPL,iCOLL,iLAST/0,1,2,3,4/
!x	data collN,collE/0.0,0.0/	       ! Number & energy of coll photons
       double precision lastN,lastE                ! Number & energy of last photons
!x	data lastN,lastE/0.0,0.0/          ! Number & energy of last photons
!x	data softN,softE/0.0,0.0/          ! Number & energy of soft photons
       LOGICAL preRIPL
       logical :: noRIPL = .FALSE.     ! Do not use RIPL tables!

!Calculate P(T) (to compare with Madland & Nix):
       dimension h(0:61); character(len=1) h
       character(len=1) :: blank = ' '

! NEW: -------------------------
       data tmax/10.0E-9/      ! Slowest half-life allowed (sec) [Talou]
       double precision E,A,S,x,aA,Af,ax,ay,az,bx,by,bz,b
       character(len=70) char70     ! headlines
!!	dimension El(1),lf(1),F(1)	!! temporary!!
! NEW. -------------------------

! We assume Gamma(n)>>Gamma(g) for Q>Qmin and Gamma(n)<<Gamma(g) for Q<Qmin:

       data EiRot,EfRot,U/0.0,0.0,0.0/

! Data for Nk=6:
!      data kZ /0,0,1,1,1,2,2/	! Charge number of various ejectile types
!      data kA /0,1,1,2,3,3,4/	! Mass number of various ejectile types
! Data for Nk=1:
       data kZ /0,0/           ! Charge number of various ejectile types
       data kA /0,1/           ! Mass number of various ejectile types

!cancel	GDRf(E)=1.0
       double precision bias0,biasN
       double precision c,cx,cy,cz
       double precision collN,collE
       double precision cosphi,costh
       double precision GDRw,GDRm,GDRf,hSsq,Q1,Q2,Q3,Q4,Q5,Q6,RRA
       double precision alevel,rng
       double precision epsmin ! Maximum heat left after statistical photons.

       data bias0,biasN/0.0,0.0/

       character(len=maxDATAPATH+100) dataf
               ! name of data file
       logical exists
       logical, save :: dataread=.FALSE.
               ! static variable storing whether data files were read in

!Copy argument values into local variables to avoid changing them:
!      phase,iZ,iA,eps,SS,PP,mm,m,id,p:
!      iZ=iZx; iA=iAx; eps=epsx
!      do i=1,3; SS(i)=SSx(i); enddo
!      do i=0,4; PP(i)=PPx(i); enddo
!      m=mx; do i=0,Nk; mult(i)=mmx(i); enddo
!      p(*)=q(*); id(*)=idx(*) not necessary since p & id initialized here

!************************************************************************
       IF (phase.eq.0) THEN ! INITIALIZE:	=========================
!Change maximum isomer half-life considered:
!      	 tmax=10.0E-9	! Standard
!	 tmax=10.0E-9	! Verbinski
         tmax=1.5E-9  ! Billnert
!	 tmax=5.0E-9	! Wang
         mp=m    ! Power of eta in photon spectrum [std:2]
#ifdef WRITEL6
         write (L6,"('DecayS:  Power of eps, mp =',i2)") mp
         write (L6,25) GDRm(0.3*iA),GDRw()
25       format('DecayS: GDRm & GRDw for A=0.3*A0:',2f8.3,' MeV')
         write (L6,26) GDRm(0.5*iA),GDRw()
26       format('DecayS: GDRm & GRDw for A=0.5*A0:',2f8.3,' MeV')
         write (L6,27) GDRm(0.7*iA),GDRw()
27       format('DecayS: GDRm & GRDw for A=0.7*A0:',2f8.3,' MeV')
         write (L6,*) 'Include GDR modulation? [0:no/1:yes]'
#endif
!        read  (L5,*) c
         c=1
         GDwidth=c*GDwidth
         if (GDwidth.eq.0.0) GDRmax=0.0
!        ------------------------------------------------------
#ifdef WRITEL6
         write (L6,*) 'Omit RIPL cascades [1/0]?'
#endif
!        read  (L5,*) n
         n=0
         if (n.ne.0) noRIPL=.TRUE.
!	 ---------------------------------------------------------
#ifdef WRITEL6
         write (L6,"('Reduction of S in statistical decays:',f6.3)") dS
#endif
!        read  (L5,*) dS
         dS=1.0
!        ---------------------------------------------------------
         dT=0.05; nT=0    ! bin width, total counts
         gmin=eps         ! Set gmin for photons
         do i=1,mT
           Ti(i)=0.0; T1(i)=0.0; T2(i)=0.0; T3(i)=0.0
         enddo
!Calculate masses of various ejectile types [here only 0 & 1]:
         do k=0,Nk
           wk(k)=kA(k)*wA+Dk(k)     ! Mass of ejectile type k
         enddo
!Counters per event (not needed & already initialized above):
         countN(iSOFT)=0.0; countE(iSOFT)=0.0  ! 0: Soft photons
         countN(iSTAT)=0.0; countE(iSTAT)=0.0  ! 1: Statistical photons
         countN(iRIPL)=0.0; countE(iRIPL)=0.0  ! 2: RIPL photons
         countN(iCOLL)=0.0; countE(iCOLL)=0.0  ! 3: Collective photons
         countN(iLAST)=0.0; countE(iLAST)=0.0  ! 4: Last photons
         Sinit=0.0    ! Pre-evaporation fragment spin magnitude
         Spost=0.0    ! Post-evaporation fragment spin magnitude
         ave0=0.0; ave1=0.0; ave2=0.0; ave3=0.0
         do i=1,3; aves(i)=0.0; enddo; ave=0.0
         do i=1,3; do j=1,3; vars(i,j)=0.0; enddo; enddo
         do i=0,Ncos; bin(i)=0.0; ben(i)=0.0; enddo
         Q0s=0.0; Q1s=0.0; Q2s=0.0; Q3s=0.0; Q4s=0.0; Q5s=0.0; Q6s=0.0
         NQf=0; aveQf=0.0
         Sst1=0.0; Sst2=0.0; Syr1=0.0; Syr2=0.0
         Est1=0.0; Est2=0.0; Eyr1=0.0; Eyr2=0.0
         minZ=100; maxZ=1
         do i=1,100; minA(i)=300; maxA(i)=1; enddo
         d0=0.0; d1=0.0; d2=0.0        ! RIPL energy shift
         E0=0.0; E1=0.0; E2=0.0        ! Highest RIPL energy
         do n=0,Nmx; isomer(n)=0; enddo
!        =================================================================
         if (dataread.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,"(40('-'),10x,'Starting reading RIPL tables:')")
!Copy over the RIPL data tables prepared by 'Gamma/levels.f':

! nA:	nA(Z): Number of isotopes for each Z:
           write (L6,*) 'Copying file Decays/nA.dat into nA(Z) ..'
#endif
           dataf=trim(freyadir)//'/Decays/nA.dat'
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 14) trim(dataf)
14           format("FREYA: data file ", a, " could not be found.")
             call seterror(14,error)
             RETURN
           endif
           OPEN  (L2,file=dataf,status='old')
           read  (L2,"(a70)") char70
#ifdef WRITEL6
           write (L6,"(a70)") char70
#endif
           read  (L2,"(2i5)") mnZ,mxZ
#ifdef WRITEL6
           write (L6,"(2i5,': minZ & maxZ')") mnZ,mxZ
#endif
           do iZ=mnZ,mxZ; read (L2,"(2i5)") i,nA(i); enddo
           CLOSE (L2)

! iAZ:	iAZ(Z,m): Mass numbers for each Z:
#ifdef WRITEL6
           write (L6,*) 'Copying file Decays/iAZ.dat into iAZ(Z) ..'
#endif
           dataf=trim(freyadir)//'/Decays/iAZ.dat'
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 15) trim(dataf)
15           format("FREYA: data file ", a, " could not be found.")
             call seterror(15,error)
             RETURN
           endif
           OPEN  (L2,file=dataf,status='old')
           read  (L2,"(a70)") char70
#ifdef WRITEL6
           write (L6,"(a70)") char70
#endif
           read  (L2,"(2i5)") mnZ,mxZ
#ifdef WRITEL6
           write (L6,"(2i5,': minZ & maxZ')") mnZ,mxZ
#endif
           do iZ=mnZ,mxZ
             do n=1,nA(iZ); read (L2,"(3i5)") i,j,iAZ(i,j); enddo
           enddo
           CLOSE (L2)

! nAZ:	Index for nucleus (Z,A):
#ifdef WRITEL6
           write (L6,*) 'Copying file Decays/nAZ.dat into nAZ(Z) ..'
#endif
           dataf=trim(freyadir)//'/Decays/nAZ.dat'
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 16) trim(dataf)
16           format("FREYA: data file ", a, " could not be found.")
             call seterror(16,error)
             RETURN
           endif
           OPEN (L2,file=dataf,status='old')
           read  (L2,"(a70)") char70
#ifdef WRITEL6
           write (L6,"(a70)") char70
#endif
           read  (L2,"(2i5)") mnZ,mxZ
#ifdef WRITEL6
           write (L6,"(2i5,': minZ & maxZ')") mnZ,mxZ
#endif
           do iZ=mnZ,mxZ
             m=iAZ(iZ,nA(iZ))-iAZ(iZ,1)+1
             do n=1,m; read (L2,"(3i5)") i,j,nAZ(i,j); enddo
           enddo
           CLOSE (L2)

! nl:	nl(N): Number of energy levels in nucleus #N:
#ifdef WRITEL6
           write (L6,*) 'Copying file Decays/nl.dat into nl(N) ..'
#endif
           dataf=trim(freyadir)//'/Decays/nl.dat'
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 17) trim(dataf)
17           format("FREYA: data file ", a, " could not be found.")
             call seterror(17,error)
             RETURN
           endif
           OPEN  (L2,file=dataf,status='old')
           read  (L2,"(a70)") char70
#ifdef WRITEL6
           write (L6,"(a70)") char70
#endif
           read  (L2,"(i5)") Ninc
#ifdef WRITEL6
           write (L6,"(i5,': Number of Nuclei')") Ninc
#endif
           do i=1,Ninc; read (L2,"(4i5)") nl(N),N,iZN(N),iAN(N); enddo
           CLOSE (L2)

! lN:	lN(N) points to the last level of nucleus #N:
#ifdef WRITEL6
           write (L6,*) 'Copying file Decays/lN.dat into lN(N) ..'
#endif
           dataf=trim(freyadir)//'/Decays/lN.dat'
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 18) trim(dataf)
18           format("FREYA: data file ", a, " could not be found.")
             call seterror(18,error)
             RETURN
           endif
           OPEN  (L2,file=dataf,status='old')
           read  (L2,"(a70)") char70
#ifdef WRITEL6
           write (L6,"(a70)") char70
#endif
           read  (L2,"(i5)") Ninc
#ifdef WRITEL6
           write (L6,"(i5,': Number of Nuclei')") Ninc
#endif
           do i=0,Ninc; read (L2,"(i5,i10)") N,lN(N); enddo
           CLOSE (L2)

! Ek:	E(l,N) & nk(l,N) & t(l,N): 
#ifdef WRITEL6
           write (L6,*) 'Copying Decays/Ek.dat into E(N,l) & nk(N,l)', &
                        ' & t(N,l)..'
#endif
           dataf=trim(freyadir)//'/Decays/Ek.dat'
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 22) trim(dataf)
22           format("FREYA: data file ", a, " could not be found.")
             call seterror(19,error)
             RETURN
           endif
           OPEN  (L2,file=dataf,status='old')
           read  (L2,"(a70)") char70
#ifdef WRITEL6
           write (L6,"(a70)") char70
#endif
           read  (L2,"(i10)") linc
#ifdef WRITEL6
           write (L6,"(i10,': Number of levels')") linc
#endif
           do i=1,linc
             read (L2,"(i10,f10.6,2i10,e11.2)") l,ElN(l),nkN(l), &
                  nkNl(l),tlN(l)
           enddo
           CLOSE (L2)

! Fij:	nf & Fij:
#ifdef WRITEL6
           write (L6,*) 'Copying Decays/Fij.dat to nf(N,l,k) &', &
                       ' Fij(N,l,k) ..'
#endif
           dataf=trim(freyadir)//'/Decays/Fij.dat'
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 23) trim(dataf)
23           format("FREYA: data file ", a, " could not be found.")
             call seterror(20,error)
             RETURN
           endif
           OPEN  (L2,file=dataf,status='old')
           read  (L2,"(a70)") char70
#ifdef WRITEL6
           write (L6,"(a70)") char70
#endif
           read  (L2,"(i10)") kinc
#ifdef WRITEL6
           write (L6,"(i10,': Number of decays')") kinc
#endif
           do NN=1,Ninc
             l0=lN(NN-1)                          ! Last level of previous nucleus
             read (L2,*) N,n0,nl(N)!,':',iZN(N),iAN(N),' ==============='
             do ll=1,nl(NN)
               read (L2,*) N,n0,l,nkN(n0+l)!,':',iZN(N),iAN(N),' ---'
               do kk=1,nkN(l0+ll)                 ! nk(N,l)=nkN(lN(N-1)+l)
                  k0=nkNl(l0+ll-1)+kk             ! Loc of this decay: k(N,l,k)
                  read (L2,"(i9,2i5,i9,i5,f12.8,2i5)") &
                        N,l,k,k0,nf(k0),Fij(k0),iZ,iA
               enddo
             enddo
           enddo
           CLOSE (L2)
           nRIPL=0                                ! Number of RIPL decays
#ifdef WRITEL6
           write (L6,"(40('-'),10x,'Finished reading RIPL tables.')")
#endif
           dataread = .TRUE.
         endif ! dataread == false

#ifdef WRITEL6
!x        write (L6,*) 'Now checking the RIPL tables:'
!x  799   write (L6,*) 'Enter Z & A:'
!x        read  (L5,*) iZ,iA
!x        i=iA-iAZ(iZ,1)+1; index=nAZ(iZ,i)	! nucleus index
!x        write (L6,"(2i4,':',3i8)") iZ,iA,index,iAZ(iZ,1),i
!x        if (iZ.gt.0) goto 799
!x        pause 'Done checking the RIPL tables.'
#endif
!        ===============================================================
         RETURN
       ENDIF              ! INITIALIZE.	=========================

!***********************************************************************
       IF (phase.lt.0) THEN ! FINALIZE:		=========================
#ifdef WRITEL6
         write (L6,"(80('='))")
#endif
         Ev=dble(-kEv)                                ! Number of events
#ifdef WRITEL6
         if (noRIPL) write (L6,*) 'No RIPL cascades were included'
         write (L6,"('Reduction of S in statistical decays:',f6.3)") dS
         write (L6,"('DecayS:  photon threshold =',f6.3,i10)") gmin,-kEv
         write (L6,"('DecayS:  maximum halflife =',1pe9.2,' s')") tmax
         write (L6,"('DecayS:',i9,' nuclei were not in RIPL')") missing
#endif
         softN=countN(iSOFT); softE=countE(iSOFT)
         statN=countN(iSTAT); statE=countE(iSTAT)
         riplN=countN(iRIPL); riplE=countE(iRIPL)
         collN=countN(iCOLL); collE=countE(iCOLL)
         lastN=countN(iLAST); lastE=countE(iLAST)
#ifdef WRITEL6
         write (L6,48)
48       format('Photons:   stat +  coll +  last +  RIPL =  total')
         write (L6,"(8x,4f8.4,f9.4)") statN/Ev,collN/Ev,lastN/Ev, &
               riplN/Ev,(statN+collN+lastN+riplN)/Ev
         if (GDwidth.ne.0.0) then
           write (L6,"('GDR sampling efficiency:',f8.5)") bias0/biasN
         else
           write (L6,"('No GDR modulation in this run')")
         endif
#endif
         k=-2*kEv
#ifdef WRITEL6
         write (L6,"('<pre-evap spin>: ',f9.4)") Sinit/k
         write (L6,"('<post-evap spin>:',f9.4)") Spost/k
         write (L6,49)
49       format('Photons:',3x,'Ngamma',9x,'Egamma',5x,'(per event)') 
         write (L6,30)softN/Ev,softE/Ev,softE/softN,int(softN+0.5),iSOFT
30       format(8x,'Nsoft=',f7.4,'; Esoft=',f7.4,'; E/N=',f7.4,2i9)
         write (L6,31)statN/Ev,statE/Ev,statE/statN,int(statN+0.5),iSTAT
31       format(8x,'Nstat=',f7.4,'; Estat=',f7.4,'; E/N=', f7.4,2i9)
         if (riplN.gt.0.0) then
           write (L6,32) riplN/Ev,riplE/Ev,riplE/riplN, &
                         int(riplN+0.5),iRIPL
32         format(8x,'Nripl=',f7.4,'; Eripl=',f7.4,'; E/N=',f7.4,2i9)
         else
           write (L6,33) riplN/Ev,riplE/Ev,int(riplN+0.5),iRIPL
33         format(8x,'Nripl=',f7.4,'; Eripl=',f7.4,13x,2i9)
         endif
         if (collN.gt.0.0) then
           write (L6,34) collN/Ev,collE/Ev,collE/collN, &
                         int(collN+0.5),iCOLL
34         format(8x,'Ncoll=',f7.4,'; Ecoll=',f7.4,'; E/N=',f7.4,2i9)
         else
           write (L6,35) collN/Ev,collE/Ev,int(collN+0.5),iCOLL
35         format(8x,'Ncoll=',f7.4,'; Ecoll=',f7.4,13x,2i9)
         endif
         if (lastN.gt.0.0) then
           write (L6,36) lastN/Ev,lastE/Ev,lastE/lastN,int(lastN+0.5), &
                 iLAST
36         format(8x,'Nlast=',f7.4,'; Elast=',f7.4,'; E/N=',f7.4,2i9)
         else
           write (L6,37) lastN/Ev,lastE/Ev,int(lastN+0.5),iLAST
37         format(8x,'Nlast=',f7.4,'; Elast=',f7.4,13x,2i9)
         endif
#endif
         if (collN+lastN.gt.0.0) &
#ifdef WRITEL6
           write (L6,38) (collN+lastN)/Ev,(collE+lastE)/Ev, &
                  (collE+lastE)/(collN+lastN),int(collN+lastN+0.5), &
                  iCOLL,iLAST
38         format('  Coll + Last:',f7.4,8x,f7.4,6x,f7.4,i9,i8,'+', i1) 
#endif

         Sst1=Sst1/k; Sst2=sqrt(Sst2/k-Sst1**2)
         Est1=Est1/k; Est2=sqrt(Est2/k-Est1**2)
         Syr1=Syr1/k; Syr2=sqrt(Syr2/k-Syr1**2)
         Eyr1=Eyr1/k; Eyr2=sqrt(Eyr2/k-Eyr1**2)
#ifdef WRITEL6
         write (L6,"('S at start of stat: ',f16.4,' +/-',f9.4)") Sst1, &
                                                                Sst2
         write (L6,"('E at start of stat: ',f16.4,' +/-',f9.4)") Est1, &
                                                                Est2
         write (L6,"('S at start of yrast:',f16.4,' +/-',f9.4)") Syr1, &
                                                                Syr2
         write (L6,"('E at start of yrast:',f16.4,' +/-',f9.4)") Eyr1, &
                                                                Eyr2
#endif
#ifdef WRITEL8
         write (L8,'(80("-"))')
         write (L8,*) 'DecayS: min & max Ap encountered for each Zp:'
         do i=minZ,maxZ
           write (L8,"(i3,':',2i5)") i,minA(i),maxA(i)
         enddo
#endif
!Calculate P(T):	---------------------------------------------------------
#ifdef WRITEL8
         write (L8,'(80("-"))')
         write (L8,*) 'Residual temperature distribution P_n(T):'
         write (L8,*) ' A: After any neutron emission'
         write (L8,*) " n: After the n'th emission in any ", &
                      "evaporation chain"
         write (L8,*) 'Normalization: Sum_{n>0} P_n(T) = P_{all n}(T)'
         write (L8,*) '   n   T(MeV)  P_n(T)'
#endif
         do i=1,mT
           do l=0,50; h(l)=blank; enddo
           T3(i)=T3(i)/nT; l=int(500*T3(i)); h(l)='3' ! 1h3
           T2(i)=T2(i)/nT; l=int(500*T2(i)); h(l)='2' ! 1h2
           T1(i)=T1(i)/nT; l=int(500*T1(i)); h(l)='1' ! 1h1
           Ti(i)=Ti(i)/nT; l=int(500*Ti(i)); h(l)='A' ! 1hA
#ifdef WRITEL8
           write (L8,800) i,(i-0.5)*dT,Ti(i),(h(l),l=1,50)
#endif
         enddo
  800    format (i5,f8.3,f10.5,1x,'|',50a1)
!Calculate P(T).	---------------------------------------------------------
         IF (bin(0).gt.0.0) THEN       !	=================
#ifdef WRITEL8
           write (L8,"(72('-'))")
           write (L8,*) 'DecayS:'
           write (L8,40)
40         format(17x,'Angular distribution of post-fission neutrons')
           write (L8,41)
41         format(17x,'relative to the direction of the emitter spin:')
           write (L8,"(' P(cos(theta)):  ',55('-'))")
!          write (L8,*) ' #  cos   yield 0.5                      1.0'	!  50*
!    c	,          '                      1.5'			!  50*
           write (L8,*) ' #  cos   yield 0.8                      1.0'  ! 125*
!    c	,          '                      1.2'			! 125*
           write (L8,*) ' #  cos   yield 0.9                      1.0', &   ! 250*
                       '                      1.1'             ! 250*
#endif
           c=0.0
           do i=1,Ncos
             costh=2.0*(i-0.5)/Ncos-1.0
             bin(i)=Ncos*bin(i)/bin(0); c=c+bin(i)
             ben(i)=Ncos*ben(i)/ben(0)
             do n=1,61; h(n)=' '; enddo; h(25)='|'; h(50)='|'
!		n=26+ 50*(ben(i)-1.0)	! 0.5 - 1.5
               n=int(26+250*(bin(i)-1.0))   ! 0.9 - 1.1
               if (n.lt. 1) then; h( 1)='-'
               else if (n.gt.60) then; h(61)='+'
               else; h(n)='x'; endif
!		n= 1+ 25*bin(i)		! 0.0 - 2.0 same as below
!		n=26+ 25*(bin(i)-1.0)	! 0.0 - 2.0 same as above
!		n=26+ 50*(bin(i)-1.0)	! 0.5 - 1.5
!		n=26+125*(bin(i)-1.0)	! 0.8 - 1.2
                n=int(26+250*(bin(i)-1.0))   ! 0.9 - 1.1
               if (n.lt. 1) then; h( 1)='<'
               else if (n.gt.60) then; h(61)='>'
               else; h(n)='*'; endif
#ifdef WRITEL8
               write (L8,"(i3,f6.2,f7.3,'  |',61a1)") i,costh,bin(i), &
                     (h(n),n=1,61)
#endif
           enddo
           i=Ncos/2; d=(bin(i)+bin(i+1))/(bin(1)+bin(Ncos))
#ifdef WRITEL8
           write (L8,42) c/Ncos,d
42         format(' Average:',f7.3,4x,'In-plane/out-of-plane:', &
                  f8.3,2x,20('-'))
#endif
           n=int(bin(0)+0.5); i=int(ben(0)+0.5); c=bin(0)/ben(0)
#ifdef WRITEL8
           write (L8,"(' x (tentative)',i10)") i
           write (L8,"(' *  (actual)  ',i10)") n
           write (L8,"('Evaporation success:',i9,' of',i9,':',f9.4)") &
                 n,i,c
#endif
! Angular momentum of tentative neutrons:	-------------------------
           do i=1,3; aves(i)=aves(i)/ave; enddo
           do i=1,3; do j=1,3; vars(i,j)=vars(i,j)/ave; enddo; enddo
           c=sqrt(vars(1,1)+vars(2,2)+vars(3,3)-aves(1)**2-aves(2)**2 &
             -aves(3)**2)
#ifdef WRITEL8
           write (L8,43)
43         format('Angular momentum of tentative neutrons:', f12.4)
           write (L8,"('avei:',f9.4,2x,3f9.4,'; sigmai:',f9.4)") &
                 (aves(i),(vars(i,j),j=1,3), &
                 sqrt(vars(i,i)-aves(i)**2), i=1,3)
#endif
           if (NQf.gt.0) aveQf=aveQf/NQf
#ifdef WRITEL8
           write (L8,"(i8,' rejects, <Qf> =',f8.4)") NQf,aveQf
           write (L8,"(72('-'))")
#endif
         ENDIF                         !      =================
!Check (not needed but interesting):
         if (ave0.gt.0.0) then
#ifdef WRITEL8
           write (L8,*) '<eps>:      emit     rotate    final    ', &
                       'emit+rot'
#endif
           w=wk(1)
           ave1=ave1/(2*w*ave0); ave2=ave2/(2*w*ave0); ave3=ave3/ &
                (2*w*ave0)
#ifdef WRITEL8
           write (L8,"('DecayS: ',4f10.4,'  (MeV)')") ave1,ave2,ave3, &
                 ave1+ave2
           write (L8,"(' speed: ',3f10.4,10x,'   (c) ')") &
                 sqrt(2*ave1/w),sqrt(2*ave2/w),sqrt(2*ave3/w)
#endif
         endif
!Check  (not needed but interesting).
         if (Q0s.gt.0.0) then
#ifdef WRITEL8
           write (L8,*) '            <P1>      <P2>      <P3>      ', &
                       '<P4>      <P5>      <P6>  '
#endif
           Q1s=Q1s/Q0s; Q2s=Q2s/Q0s; Q3s=Q3s/Q0s
           Q4s=Q4s/Q0s; Q5s=Q5s/Q0s; Q6s=Q6s/Q0s
#ifdef WRITEL8
           write (L8,"('<P(1-6)>:',6f10.5)") Q1s,Q2s,Q3s,Q4s,Q5s,Q6s
           write (L8,"('<cos(theta)**2> =',f10.5,' (actual)')") &
                 (2*Q2s+1)/3
           write (L8,"('<cos(theta)**2> =',f10.5,' (isotropic)')") 1.0/3
#endif
         endif
!Cut-off in photon energy:
#ifdef WRITEL6
         write (L6,"('Minimum photon energy recorded:',f8.4)") gmin
#endif
!Calculate the distribution of the number of used RIPL energy levels:
#ifdef WRITEL6
         write (L6,"(27x,'ave      disp       var')")
#endif
         if (h0.gt.0.0) then
           h1=h1/h0; h2=h2/h0-h1**2
#ifdef WRITEL6
           write (L6,"('Number of RIPL levels:',3f10.4)") h1,sqrt(h2),h2
#endif
         endif
!Calculate the distribution of the maximum energy of RIPL levels used:
         if (E0.gt.0.0) then
           E1=E1/E0; E2=E2/E0-E1**2
#ifdef WRITEL6
           write (L6,"('Maximum RIPL energy:  ',3f10.4)") E1,sqrt(E2),E2
#endif
         endif
!Calculate the distribution of the photon energy shifts needed by RIPL:
         if (d0.gt.0.0) then
           d1=d1/d0; d2=d2/d0-d1**2
#ifdef WRITEL6
           write (L6,"('RIPL energy shift:    ',3f10.4)") d1,sqrt(d2),d2
#endif
         endif
!Count number of RIPL decays:
#ifdef WRITEL6
         write (L6,"(i9,' RIPL decays =',f8.4,' per event')") &
               nRIPL,nRIPL/Ev
#endif
!Compilation of isomers encountered:
         if (isomer(0).gt.0) then
#ifdef WRITEL8
           write (L8,'(80("-"))')
           write (L8,44) isomer(0),tmax
44         format(i9,' RIPL levels lived longer than',1pe9.2,' s:')
           write (L8,*) '  Zp   Ap     index     N '
#endif
           isum=0; isomer(0)=0           ! to list only RIPL nuclei
           do n=1,Nmx
             iZ=iZN(n); iA=iAN(n)
             i=iA-iAZ(iZ,1)+1; index=nAZ(iZ,i) ! nucleus index
             if (isomer(index).gt.0) isum=isum+isomer(index)
#ifdef WRITEL8
             if (isomer(index).gt.0) write (L8,"(2i5,':',3i8)") &
               iZ,iA,index,isomer(index),isum
#endif
           enddo
         endif
         RETURN
       ENDIF         ! FINALIZE.		=========================
!************************************************************************

!================================     DecayS:	=========================

       iA0=iA                                  ! Ur-mother mass number
       epsi=eps                                ! Ur-mother excitation (heat)
       epsi0=epsi                              ! Ur-mother excitation (heat)
!	Di=DROP(iZ,iA)                         ! Ur-mother mass excess (LD)
       Di=gsM (iZ,iA)                          ! Ur-mother mass excess
       Wi=iA*wA+Di                             ! Ur-mother ground-state mass
       S2=SS(1)**2+SS(2)**2+SS(3)**2           ! Spin squared of urmother
       Sinit=Sinit+sqrt(S2)                    ! -> <pre-evap fragment spin>
!Class	EiRot=0.5*S2/ROT(iA)                   ! Rotational energy of urmother
       EiRot=hSsq(sqrt(S2))/ROT(iA)            ! Rotational energy of urmother
       Yi=Wi+EiRot                             ! Urmother yrast mass: gs+rot
!Choose among the competing decay channels k = gamma, n, p, d, t, tau, alpha:
!	eta=rng(iseed)	! to select ejectile type k = 0,...,Nk
!Consider ONLY neutrons (followed by photons) in this version:
! NEUTRON emission:	-------------------------------------------------
       k=1                                     ! k=1: Ejectile is a neutron
       max=mult(k)                             ! max # k-emissions in total
       n=0                                     ! -> actual k multiplicity
   10  if (n.eq.max) goto 19                   ! Max # k-emissions is reached
!      By specifying max one can control how many neutrons may be evaporated;
!      for example: none at all, or at most one, which can be instructive.
       Ei=Yi+epsi                              ! Mother rest erg: gs+rot+heat
       w=wk(k)                                 ! Ejectile mass
       iZf=iZ-kZ(k); Zf=dble(iZf)             ! Daughter charge number
       iAf=iA-kA(k); Af=dble(iAf)             ! Daughter mass number
!      Df=DROP(iZf,iAf)			! Daughter mass excess (LD)
       Df=gsM (iZf,iAf)                        ! Daughter mass excess
       Wf=iAf*wA+Df                            ! Daughter ground-state mass
       Sn=Df+Dn-Di                             ! Neutron separation energy
! The Q value corresponds to emission of an ultrasoft point particle
! along the spin direction which leaves S unaffected (so Sf=Si)
!   => EfRot = Sf2/ROTf > EiRot = Si2/ROTi because ROTf < ROTi.
       S2=SS(1)**2+SS(2)**2+SS(3)**2           !> If SS is conserved then this
!Class EfRot=0.5*S2/ROT(iAf)                   !> is the rot erg of daughter,
       EfRot=hSsq(sqrt(S2))/ROT(iAf)
       Yf=Wf+EfRot                             ! Daughter yrast mass: gs+rot
!Change:	Qk=Ei-Yf-w				! Q value for evaporation (S=0)
!      Qk=epsi-Sn & Sn=Df+Dn-Di => Qk=Di+epsi-Df-Dk(k)
       Qk=epsi-Sn                              ! Q value (~200x more accurate)
#ifdef WRITEL6
!      write (L6,"(i4,3f12.3,8x,4f8.3)") 
!      c      -iAf,Wf,EfRot,Yf,epsi,Di+epsi-Df-Dk(k),Qk!!,sqrt(2*w*Qk)!!
!      write (L6,"(i4,3f12.3,24x,2f8.3)") -iAf,Wf,EfRot,Yf,Qk!!,sqrt(2*w*Qk)!!
       if (Di+Df.ge.500.) then
         write (L6,*) iZ,iA,' -> ',iZf,iAf,Di,Df,epsi,m
!        pause !!
       endif
#endif
!Check Gamma(neutron)/Gamma(gamma) (using Woods-Saxon with 1 MeV width):
! 	step=1.0				! WS width
! 	if (rng(iseed).lt.exp(-(Qk-Qmin)**2/step**2)) goto 19
! 
!Cut-off approximation to P(n)=Gamma(n)/[Gamma(n)+Gamma(g)] is used here:
!************************************************************************
!      Neutron emission occurs if the Q value exceeds the specified Qmin.
       if (Qk.lt.Qmin) goto 19                 ! Evap is forbidden -> no more evap
!      NOTE: Qk is usually negative (~-Sn/2) when evaporation is abandoned,
!      but for numerical reasons evap is also abandoned for 0 < Qk < Qmin.
!************************************************************************

!Carry out emission of one (more) neutron (this is possible because Qk>0):
       n=n+1                                   ! Increment k multiplicity
       m=m+1                                   ! Increment total multiplicity
!Checking that there is array space for one more ejectile:
       if (m.gt.mMax) then
#ifdef WRITEL6
         write (L6,*) 'Emitting neutrons in DecayS:'
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
!Checking that there is array space for one more ejectile.
       id(m)=k                                 ! Store identity of ejectile
!      egs=ecR(gsC(iZf,iAf)/RRA(iAf))          ! Eccentricity of ground state
!      aA=alevel(iK,iZf,iAf,Qk,egs)
!Calculate the level-density parameter:
       aA=alevel(iK,iZf,iAf,Qk,U)          ! alevel(iK,iZ,iA,eps,ecc,U) => U
       T=sqrt(U/aA)                            ! Maximum temp in daughter

!Calculate P(T):	=========================================================
!      if (id(mMax).eq.1.AND.n.ge.1)
       nT=nT+1; i=int(1+T/dT)                  ! i: temperature bin
       if (i.le.mT) then
         Ti(i)=Ti(i)+1.0                       ! After any neutron evap;
         if (n.eq.1) T1(i)=T1(i)+1.0           ! after 1st neutron evap;
         if (n.eq.2) T2(i)=T2(i)+1.0           ! after 2nd neutron evap;
         if (n.eq.3) T3(i)=T3(i)+1.0           ! after 3rd neutron evap.
       endif
!Calculate P(T).	=========================================================

!************************************************************************
! NEUTRON emission:	-------------------------------------------------
!Choose the tentative kinetic energy of relative motion (ie tot kin erg in CM):
!      This is all tentative because the emission affects the rotational 
!      energy, so the particular selection of p may not be permitted by
!      energy conservation, in which case the selection must be repeated.
!			-------------------------------------------------
       Ntry=0
   12  eps12=Qk        ! Maximum possible (relative) neutron energy:
       do while (eps12.lt.Emin.OR.eps12.ge.Qk) ! Numerics: Ek > Emin required 
!      eps12 < Emin:	The kinetic energy is below the specified limit!
!      eps12 > Qk:	The kinetic exceeds the amount available!
!      Note: eps12 is the relative kinetic energy of ejectile and daughter
         eta1=-log(rng(iseed))                ! P(eta1)=exp(-eta1)
         eta2=-log(rng(iseed))                ! P(eta2)=exp(-eta2)
         eps12=(eta1+eta2)*T                   ! Tot kin energy (in 2-body CM)
       enddo
       Qm(m)=Qk; epsm(m)=eps12                 ! to print if mult > mMax
       Qf=Qk-eps12                             ! Tentative heat content of daughter > 0
       Ef=Yf+Qf                                ! Tentative daughter rest erg: gs+rot+heat
! The above tentative value of Qf ignores the change in rotational energy!
!Calculate the magnitude of the local ejectile momentum pk (relativistic):
!      pk=0.5*sqrt((Ei-Ef-w)*(Ei+Ef+w)*(Ei-Ef+w)*(Ei+Ef-w))/Ei -> 
       pk=0.5*sqrt(eps12*(eps12+2*w)*(Ei+Ef+w)*(Ei+Ef-w))/Ei ! is pk*c: erg
!Choose the emission point r=(x,y,z) randomly on the spherical nuclear surface:
       RA=RRA(iA)                              ! Nuclear radius
       phi=twopi*rng(iseed)                    ! Azimuthal angle phi
       sinphi=sin(phi)                         ! sin(azimuthal angle phi)
       cosphi=cos(phi)                         ! cos(azimuthal angle phi)
       costh=2.0*rng(iseed)-1.0                ! cos(polar angle theta)
       sinth=sqrt(1.-costh**2)                 ! sin(polar angle theta)
       x=RA*sinth*cosphi                       !>
       y=RA*sinth*sinphi                       ! > r = (x,y,z)
       z=RA*costh                              !>
!Construct the local reference system abc:
       ax= cosphi*costh; ay=sinphi*costh; az=-sinth
       bx=-sinphi; by=cosphi; bz=0.0
       cx= cosphi*sinth; cy=sinphi*sinth; cz= costh
! This is an alternative abc system that works just as well:
!      ax=SS(2)*cz-SS(3)*cy; ay=SS(3)*cx-SS(1)*cz; az=SS(1)*cy-SS(2)*cx
!      a=sqrt(ax**2+ay**2+az**2); ax=ax/a; ay=ay/a; az=az/a
!      bx=ay*cz-az*cy; by=az*cx-ax*cz; bz=ax*cy-ay*cx
!Check that the abc system is orthonormal:
!      ab=ax*bx+ay*by+az*bz; a=sqrt(ax**2+ay**2+az**2); a*b & |a|
!      ac=ax*cx+ay*cy+az*cz; b=sqrt(bx**2+by**2+bz**2); a*c & |b|
!      bc=bx*cx+by*cy+bz*cz; c=sqrt(cx**2+cy**2+cz**2); b*c & |c|
#ifdef WRITEL6
!      write (L6,*) ' a:', a,' b:',b,' c:',c,' =0?'
!      write (L6,*) 'ab:',ab,'ac:',ac,'bc:',bc,' =0?'; pause
#endif
!Choose the direction u of the ejectile motion (biased with u*c):
       phi=twopi*rng(iseed)                    ! Azimuthal angle phi (0,2pi]
       eta=rng(iseed)                          ! cos(polar angle theta)**2
       costh=sqrt(eta)                         ! cos(polar angle theta) (0,1]
       sinth=sqrt(1.-eta)                      ! sin(polar angle theta) [0,1)
       a=sinth*cos(phi); b=sinth*sin(phi); c=costh ! dir of p in abc syst
!Components of emitted neutron in local comoving system:
       px=pk*(a*ax+b*bx+c*cx)                  !>  Ejectile momentum
       py=pk*(a*ay+b*by+c*cy)                  ! > p=(px,py,pz) [pc: energy]
       pz=pk*(a*az+b*bz+c*cz)                  !>  in comoving frame
       p2=px**2+py**2+pz**2
       ave0=ave0+1.0                           ! Number of tentative neutrons
       ave1=ave1+p2                            ! -> pre-boost kinetic energy
!Calculate the rotation vector omega = (ox,oy,oz) [omega ~ S for a sphere]:
       ox=SS(1)/ROT(iA)/hbarc                  !>   This is omega/c
       oy=SS(2)/ROT(iA)/hbarc                  ! >  given in inverse fermi:
       oz=SS(3)/ROT(iA)/hbarc                  !>   [omega/c]=(1/T)/(L/T)=1/L
!Collective velocity w =(wx,wy,wz) at the emission point r=(x,y,z): w = o x r
       wx=oy*z-oz*y; wy=oz*x-ox*z; wz=ox*y-oy*x! w=(wx,wy,wz) is w/c: dimless
       ave2=ave2+(wx**2+wy**2+wz**2)*w**2      ! -> boost energy
!Components of emitted neutron after rotational boost:
       px = px + w*wx                          !>  Ejectile momentum
       py = py + w*wy                          ! > p=(px,py,pz) [pc: energy]
       pz = pz + w*wz                          !>  after rotational boost
       p2 = px**2+py**2+pz**2
       ave3=ave3+p2                            ! -> post-boost kinetic energy
       En=sqrt(w**2+p2)                        ! Total ejectile energy in CM
!      epsk=En-w                               ! Kinetic energy of ejectile
       epsk=0.5*p2/w*(1-0.25*p2/w**2)          ! Neutron kinetic energy in CM


!Calculate ang dist of the tentative neutrons relative to the mother spin SS:
       if (p2*S2.gt.0.0) then   ! -----------------------
         costh=(px*SS(1)+py*SS(2)+pz*SS(3))/sqrt(p2*S2)
         i=int(1+0.5*Ncos*(1.0+costh)); ben(0)=ben(0)+1.0
         ben(i)=ben(i)+1.0
       endif                    ! -----------------------
!Calculate ang dist of the tentative neutrons relative to the mother spin SS.

!Calculate the angular momentum of the tentative ejectile:
       sx=(y*pz-z*py)/hbarc                                       !>	Angular momentum
       sy=(z*px-x*pz)/hbarc                                       ! >  	s = r x p
       sz=(x*py-y*px)/hbarc                                       !>	of the ejectile
       aves(1)=aves(1)+sx; aves(2)=aves(2)+sy; aves(3)=aves(3)+sz !
       vars(1,1)=vars(1,1)+sx*sx; vars(1,2)=vars(1,2)+sx*sy       !>
       vars(1,3)=vars(1,3)+sx*sz; vars(2,1)=vars(2,1)+sy*sx       ! >
       vars(2,2)=vars(2,2)+sy*sy; vars(2,3)=vars(2,3)+sy*sz       !  > correls
       vars(3,1)=vars(3,1)+sz*sx; vars(3,2)=vars(1,2)+sz*sy       ! >
       vars(3,3)=vars(3,3)+sz*sz; ave=ave+1.0                     !>
!Calculate the tentative spin of the daughter nucleus:
       SSx=SS(1)-sx                                               !>	Tentative spin
       SSy=SS(2)-sy                                               ! > 	S' = S - s
       SSz=SS(3)-sz                                               !>	of the daughter
       Sf2=SSx**2+SSy**2+SSz**2                                   !	(daughter spin)**2
!Check energy bound (non-relativistically):
!Class	EfRot=0.5*Sf2/ROT(iAf)                                    ! Tentative daughter rot energy
       EfRot=hSsq(sqrt(Sf2))/ROT(iAf)                             ! Tentative daughter rot energy
       Yf=Wf+EfRot                                                ! Yrast energy of the daughter
       Efkin=(p2/w+p2/Wf)/2                                       ! Tentative final kin energy
!      Qf=EiRot+epsi-Sn-EfRot-Efkin                               ! same as:	
       Qf=Qk+(EiRot-EfRot)-Efkin                                  ! Tentative heat in daughter
       if (Qf.le.0.0) then ! this emission violates energy conservation:
         NQf=NQf+1; aveQf=aveQf+Qf
#ifdef WRITEL6
!C         write (L6,"(i8,':',5f10.4,i8)") kEv,Qf,epsi,EfRot,Efkin,T,NQf
#endif
! This problem occurs only relatively rarely but is difficult to avoid,
! because (I think) it arises from the inconsistency in the moments of inertia:
! The moment of inertia of the daughter nucleus plus that of the emitted
! nucleon is NOT equal to the moment of inertia of the mother nucleus,
! at least in part because the nuclear moments of inertia are reduced
! relative to the rigid-body values, whereas that of the nucleon is equal to
! the "rigid-body" value (i.e. the unreduced free nucleon mass is used).
! So the 'brute-force' fix is to give up after a number (100) tries;
! because this problem occurs only when the Q-value is very small,
! there should be no discernable effect on the calculated results.
         Ntry=Ntry+1                       ! Counts number of tries
         if (Ntry.le.100) goto 12          ! Qf<0: select another p!
         n=n-1                             ! Backshift k multiplicity
         m=m-1                             ! Backshift total multiplicity
         goto 19                           ! Abandon neutron emission
       endif                               ! this emission violates energy conservation.

!Calculate ang dist of the actual neutrons relative to the mother spin SS:
       if (p2*S2.gt.0.0) then              ! -----------------------
         costh=(px*SS(1)+py*SS(2)+pz*SS(3))/sqrt(p2*S2)
         i=int(1+0.5*Ncos*(1.0+costh))
#ifdef WRITEL6
!      if (i.lt.1.OR.i.gt.Ncos) write (L6,*) kEv,costh,i,Ncos
#endif
         bin(0)=bin(0)+1.0; bin(i)=bin(i)+1.0
         Q0s=Q0s+1.0
         Q1s=Q1s+Q1(costh)
         Q2s=Q2s+Q2(costh)
         Q3s=Q3s+Q3(costh)
         Q4s=Q4s+Q4(costh)
         Q5s=Q5s+Q5(costh)
         Q6s=Q6s+Q6(costh)
       endif                                   ! -----------------------
!Calculate ang dist of the actual neutrons relative to the mother spin SS.

!      Ef*=EfRot+Qf                            ! Total excitation in daughter
       epsf=Qf                                 ! Actual heat in daughter
!Comment out the next line to omit the spin recoil from evaporation:
       SS(1)=SSx; SS(2)=SSy; SS(3)=SSz; S2=Sf2 ! Actual daughter spin

! Boost the momenta to the overall reference frame:
       g0 =PP(4)/PP(0)                         ! Lorentz boost coefficient
       VVx=PP(1)/PP(4)                         !>
       VVy=PP(2)/PP(4)                         ! > Velocity of mother nucleus
       VVz=PP(3)/PP(4)                         !>
       V2=VVx**2+VVy**2+VVz**2                 ! V * V
       pdotV = px*VVx + py*VVy + pz*VVz        ! p * V
       g1=0.0; if (V2.gt.0.0) g1=(g0-1.0)*pdotV/V2
! Boost the ejectile:
       cc=g0*En+g1                             ! cc: convenient coefficient
       p(1,m) = cc*VVx + px                    !>
       p(2,m) = cc*VVy + py                    ! > Momentum of ejectile
       p(3,m) = cc*VVz + pz                    !>  (in the absolute frame)
! Boost the daughter:
       PP(0)=Ef                                ! 0: Total daughter rest energy
!      Eftot=sqrt(Ef**2+pk**2)                 ! Total daughter energy in CM
       Eftot=sqrt(Ef**2+p2)                    ! Total daughter energy in CM
       PP(4) = g0*(Eftot-pdotV)                ! 4: Total daughter energy
       cc=g0*Eftot-g1                          ! cc: convenient coefficient
       PP(1) = cc*VVx - px                     !>
       PP(2) = cc*VVy - py                     ! > 1-3: Daughter momentum
       PP(3) = cc*VVz - pz                     !>
!checkg	gk1=PP(1)/PP(4); gk2=PP(2)/PP(4); gk3=PP(3)/PP(4)   ! check gamma 
! Daughter becomes next mother, f -> i:
       iZ=iZf                                  ! Charge number (n => same Z)
       iA=iAf                                  ! Mass number
       Di=Df                                   ! Mass excess
       Wi=Wf                                   ! Ground-state mass
       EiRot=EfRot                             ! Rotational energy
       Yi=Yf                                   ! Yrast "mass": gs+rot
       epsi=epsf                               ! Statistical energy: heat
       goto 10                                 ! Repeat evaporation procedure
!Complete the neutron evaporation chain:
   19  mult(k)=n                               ! k-type multiplicity
       Wf=Wi                                   ! gs mass after evaps
       eps=epsi                                ! -> OUTPUT		! Heat remaining before photons
       EfRot=EiRot                             ! Residue rotational erg

!************************************************************************
! PHOTON emission:	-------------------------------------------------
! Method:
!	First radiate the statistical energy, then follow the yrast band down.
!	Note:  All photons are emitted isotropically (in the emitter frame).
! STATISTICAL photons:
! At each step, calculate the nuclear temperature from Estar (= epsmax),
! then select the photon energy epsk and reduce the excitation: E* -> E*-epsk;
! if epsk < gmin then don't record this photon but do iterate the procedure;
! stop the statistical cascade when epsmax < epsmin; S remains unchanged.
! YRAST photons:
! Truncate the spin S to an even integer J = 2*[S/2] < S and emit a photon
! with the corresponding freed-up energy: epsk = eps+S**2/2Irot-J**2/2Irot.
! As long as J>0, perform an yrast emission: epsk = [J**2-(J-2)**2]/2Irot
! and J:=J-2; iterate the yrast cascade until J=0.
!			-------------------------------------------------
       if (iZ.lt.minZ) minZ=iZ
       if (iZ.gt.maxZ) maxZ=iZ
       if (iA.lt.minA(iZ)) minA(iZ)=iA
       if (iA.gt.maxA(iZ)) maxA(iZ)=iA
       eps0=eps                                ! Keep value for DecayS exit
! eps0:	Statistical excitation (heat) after neutron evaporation.
! eps:	Statistical excitation (heat) at any stage of the deexcitation chain.
! NOTE:	eps is typically ~Sn/2 at this stage.
       S=sqrt(S2)                              ! |S| after evaporation
       Ei=Wf+EfRot                             ! Total mass before photons
       A=dble(iA)                             ! Mass number of emitter
       k=0                                     ! k=0: Ejectile is a photon
       n=0                                     ! -> actual k multiplicity
       Spost=Spost+S                           ! <post-evap spin magnitude>

!Caculate the RIPL index of this product nucleus:
       i=iA-iAZ(iZ,1)+1; index=nAZ(iZ,i)       ! nucleus index
       if (index.eq.0) missing=missing+1
#ifdef WRITEL6
!      if (index.eq.0) write (L6,"('Not in RIPL:',2i4)") iZ,iA
#endif
       if (noRIPL) index=0
!      index=0!!!	  Uncomment this line to restore pre-RIPL DecayS

!========================================================================
! The pre-RIPL treatment below was copied from the old DecayS routine,
! except that S -> S-dS after each statistical emission.
!Current method:
! Statistical photon emission until eps < gmin then yrast cascade;
! whenever the total excitation drops below the highest RIPL energy,
! the remaining decays are based on the RIPL table.
! Note:
!   epsmin: the minimum heat for which stasticial emission is done;
!   gmin:   the detection threshold below which photons are unrecorded.

       S0=S                        ! |S| before any radiation
       max=mult(k)                 ! Specified max # k-emissions
       J=-1                        ! Flag: no yrast emiss yet
       koll=0                      ! Number of collective photons
       epstot=eps+EfRot            ! Combined stat & rot excitation energy
       if (epstot.lt.gmin) goto 29 ! Any decays will be too soft to record

       if (index.eq.0.0) then
         E=0.0                     ! Energy of the highest considered RIPL level
         l=0
       else
         l0=lN(index-1)      ! Last level of previous nucleus
         l=nl(index)         ! Number of tabulated levels in this nucleus,
                             ! equal to level number of the highest level.
!Consider only levels l having E(l) < Emax:
!        Emax=0.01
!        do while (l.gt.2.AND.ElN(l0+l).gt.Emax); l=l-1; enddo
         E=ElN(l0+l)         ! Energy of the highest considered RIPL level
       endif
       h0=h0+1.0; h1=h1+l; h2=h2+l**2
       E0=E0+1.0; E1=E1+E; E2=E2+E**2
       preRIPL=.TRUE.        ! Always make at least one pre-RIPL emission

   20  if (n.eq.max) goto 29 ! Max # k-emissions is reached
!      By specifying max one can control how many photons may be evaporated
!      which can be instructive; for example: none at all, or at most one.
       epsmin=epsmink(iK)
       gmin=gmink(iK)
!===============================================================================
       IF (preRIPL) THEN ! pre-RIPL decay process:
!===============================================================================
         if (eps.gt.epsmin) then   ! Statistical emission:		     -----------
           iCOUNT=iSTAT    ! This is a statistical photon
           if (n.eq.0) then ! first statistical photon
!            Sst1=Sst1+S; Sst2=Sst2+S**2
             Est1=Est1+EfRot+eps; Est2=Est2+(EfRot+eps)**2
           endif
           aA=alevel(iK,iZ,iA,eps,U)     ! => effective excitation U
           T=sqrt(U/aA)                  ! T before photon emission
           epsk=eps                      ! Max kin erg energy release
           GDRmax=1.0; if (epsk.lt.GDRw()) GDRmax=GDRf(A,epsk)
           do while (epsk.ge.eps)        ! We must demand epsk < eps
   21        eta=0.0
             biasN=biasN+1.0             ! Number of GDR tries
             do i=0,mp                   ! mp=2: black-body emission
               eta=eta-log(rng(iseed))  ! P(eta)~eta**mp/exp(eta)
             enddo
!------------------------------------------------------------------------
!Consider the Giant-Dipole Resonance form factor:
! The energy of the proposed emission is E=eta*T; its maximum is Emax=epsk.
! The form factor is GDRf(E) = Gamma*2*E**2 /[(E**2-Gamma**2)**2+Gamma**2*E**2].
!                               <  GDRf(E<<Gamma) = E**2/Gamma**2
! which has the forms           <  GDRf(E= Gamma) = 1
!                               <  GDRf(E>>Gamma) = Gamma**2/E**2
! For E>0 it lies in (0,1], so it can be interpreted as a probability: the
! statistical decay is accepted if ran < GDRf(E), hence rejected if ran > GDRf(E).
! But this is (highly) inefficient when E << Gamma because then GDRf(E) << 1.
! Therefore we normalize by its maximum value, Paccept(E)=GDRf(E)/GDRmax,
! where GDRmax=GDRf(Emax) for Emax < Gamma and GDRmax=1.0 otherwise:
             if (rng(iseed)*GDRmax.gt.GDRf(A,eta*T)) goto 21
!------------------------------------------------------------------------
             bias0=bias0+1.0             ! Number of GDR successes
!0           eta=-log(rng(iseed))       ! P(eta)=exp(-eta)
!1           eta=eta-log(rng(iseed)) 	 ! P(eta)~eta/exp(eta)
!2           eta=eta-log(rng(iseed))    ! P(eta)~eta**2/exp(eta)
!3           eta=eta-log(rng(iseed))    ! P(eta)~eta**3/exp(eta)
             epsk=eta*T                  ! Kinetic energy (in 2-body CM)
           enddo
           eps=eps-epsk                  ! Remaining heat in product
!Convert some rotational energy to statistical excitation:
           EiRot=EfRot                   ! Erot before radiation
           epstot=eps+EiRot              ! is not relevant here
           S=S-dS                        ! E1,2 reduces S by 1,2 units
           if (S.lt.0.0) S=0.0           ! No more rotation
           EfRot=hSsq(S)/ROT(iA)         ! Erot after radiation
#ifdef WRITEL6
!          write (L6,*) S,EfRot,ROT(iA)
#endif
           eps=eps+EiRot-EfRot           ! Statistical energy after
         else if (S.ge.2.0) then         ! Collective (yrast) emission:	     -----------
           iCOUNT=iCOLL                  ! This is a  collective photon
!Classical epsk = [S**2-(S-2)**2]/2*Irot = 2*(S-1)/Irot (using Erot~S**2):
!Class	  epsk=2*(S-1.0)/ROT(iA)         ! yrast emission (classical)
           EiRot=EfRot                   ! rotational energy before
           if (koll.eq.0) then ! first collective photon:
             Syr1=Syr1+S; Syr2=Syr2+S**2
             Eyr1=Eyr1+EiRot; Eyr2=Eyr2+EiRot**2
           endif
           koll=koll+1
           S=S-2.0                       ! reduce the spin: J:=J-2
!Class     EfRot=0.5*S**2/ROT(iA)        ! rotational energy after
           EfRot=hSsq(S)/ROT(iA)         ! rotational energy after
           epsk = EiRot-EfRot            ! photon energy
         else if (S.gt.0.0) then         ! Emit a last photon:
           iCOUNT=iLAST                  ! This is the last photon
           epsk=EfRot+eps                ! photon energy
           S=0.0; eps=0.0; EfRot=0.0
         else
           goto 29
         endif                           ! Done with stat & coll emission.    -----------
         epstot=eps+EfRot                ! Total excitation energy
         if (epstot.lt.E) then           ! This is the last pre-RIPL decay - match it:
           eps=epstot                    ! Excitation energy E*
           do while (E.gt.eps)
             l=l-1; E=ElN(l0+l)          ! Find the level l just below E*
           enddo
! 	  l is now the level just below E*
           d=eps-E                       ! Add'l photon energy needed to match
           d0=d0+1.0; d1=d1+d; d2=d2+d**2
           epsk=epsk+d                   ! Increase the photon energy
           eps=E                         ! Put E* to E(l) = E* - d
           preRIPL=.FALSE.               ! No more pre-RIPL emissions
#ifdef WRITEL6
!       	  write (L6,"(3i5,2f9.5,i5,f9.5,' +',f9.5,' =',f9.5)")
!      c	    iZ,iA,index,eps,E,l,epsk-d,d,epsk
#endif
         endif                           ! This is the last pre-RIPL decay - match it.
!===============================================================================
       ELSE                              ! RIPL decay process:
!===============================================================================
!Check the half-life of the mother level:
         if (tlN(l0+l).gt.tmax) then             ! Isomeric state:
           isomer(index)=isomer(index)+1         ! Count for this nucleus
           isomer(0)=isomer(0)+1                 ! Count for all nuclei
#ifdef WRITEL6
!	   write (L6,"(2i4,':',i5,f9.5,1pe11.2,' ISOMER',3i7)")
!     c	  	  iZ,iA,l,ElN(l0+l),tlN(l0+l),index,isomer(index),isomer(0)
#endif
           goto 29                               ! Done
         endif                                   ! Isomeric state.
         iCOUNT=iRIPL                            ! This is a RIPL photon
         eta=rng(iseed)
         j=1                                     ! First decay branch
         do while (Fij(nkNl(l0+l-1)+j).lt.eta)
#ifdef WRITEL6
!	   write (L6,"(3x,2i5,2f10.6,i5)")
!     c	   l,j,Fij(nkNl(l0+l-1)+j),eta,nf(nkNl(l0+l-1)+j)
#endif
           j=j+1                                 ! Next decay branch
         enddo
#ifdef WRITEL6
!	 write (L6,"(3x,2i5,2f10.6,i5)")
!     c	 l,j,Fij(nkNl(l0+l-1)+j),eta,nf(nkNl(l0+l-1)+j)
#endif
         l=nf(nkNl(l0+l-1)+j)                    ! Daughter level number
         epsk=E-ElN(l0+l)                        ! Energy release: Ei-Ef
         E=ElN(l0+l); eps=E                      ! Daughter level energy
         epstot=eps
#ifdef WRITEL6
!	 write (L6,"('Selected decay branch:',6x,2i5,2f9.5)") j,l,E,epsk
#endif
         nRIPL=nRIPL+1                           ! One more RIPL decay
!===============================================================================
       ENDIF                                     ! End of decay process.
!===============================================================================

!Count this photon and do the emission kinematics:
       n=n+1                                   ! Increment k multiplicity
       m=m+1                                   ! Increment total multiplicity
#ifdef WRITEL6
!      write (L6,"(2i5,' emits:',i5,f9.5)") iZ,iA,n,epsk
#endif
       if (m.gt.mMax) then
#ifdef WRITEL6
         write (L6,*) 'Emitting photons in DecayS:'
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
       id(m)=k                                   ! Store identity of ejectile
       Qm(m)=eps; epsm(m)=epsk                   ! To print if mult > mMax
       Ef=Wf+EfRot+eps                           ! Total mass after: gs+rot+heat
!Calculate the magnitude of the ejectile momentum pk:
       pk=0.5*epsk*(Ei+Ef)/Ei                    ! Recall: mass(photon) = 0
       Ek=pk                                     ! Photon tot=kin energy in CM
!Choose the direction of the ejectile motion:
       phi=twopi*rng(iseed)                      ! Azimuthal angle phi
       costh=2.*rng(iseed)-1.                    ! cos(polar angle theta)
       sinth=sqrt(1.-costh**2)                   ! sin(polar angle theta)
       px=pk*sinth*cos(phi)                      !>
       py=pk*sinth*sin(phi)                      ! > Ejectile momentum
       pz=pk*costh                               !>
! Boost the momenta to the overall reference frame:
       g0 =PP(4)/PP(0)                           ! Lorentz boost coefficient
       VVx=PP(1)/PP(4)                           !>
       VVy=PP(2)/PP(4)                           ! > Velocity of mother nucleus
       VVz=PP(3)/PP(4)                           !>
       V2=VVx**2+VVy**2+VVz**2                   ! V * V
       pdotV = px*VVx + py*VVy + pz*VVz          ! p * V
       g1=0.0; if (V2.gt.0.0) g1=(g0-1.0)*pdotV/V2
! Boost the photon:
       cc=g0*Ek+g1                               ! cc: convenient coefficient
       epsk = g0*(Ek+pdotV)                      ! Photon tot=kin energy in LAB
       p(1,m) = cc*VVx + px                      !>
       p(2,m) = cc*VVy + py                      ! > Momentum of photon in LAB
       p(3,m) = cc*VVz + pz                      !>
! Boost the daughter:
       PP(0)=Ef                                  ! 0: Total daughter rest energy
       Eftot=sqrt(Ef**2+pk**2)                   ! Total daughter energy in CM
       PP(4) = g0*(Eftot-pdotV)                  ! 4: Total daughter energy
       cc=g0*Eftot-g1                            ! cc: convenient coefficient
       PP(1) = cc*VVx - px                       !>
       PP(2) = cc*VVy - py                       ! > 1-3: Daughter momentum
       PP(3) = cc*VVz - pz                       !>

!Check photon LAB energy: Do not record this photon if it is too soft: 
       if (epsk.lt.gmin) then                    ! Too soft to record:
         iCOUNT=iSOFT                            ! This is a soft photon
         n=n-1; m=m-1                            ! Unrecord this photon
       endif
!Count number & energy of this photon:
       countN(iCOUNT)=countN(iCOUNT)+1.0         ! Number
       countE(iCOUNT)=countE(iCOUNT)+epsk        ! Energy

! Daughter state becomes next mother state:
       Ei=Ef                                     ! Total mass (i.e. gs+rot+heat)
       if (eps.gt.gmin) goto 20                  ! Consider another decay

!Complete photon decay chain:	!---------------------------------------
   29  mult(k)=n                                 ! k-type multiplicity
       SS(1)=0.0; SS(2)=0.0; SS(3)=0.0           ! Product nucleus has S:=0
! Note: This is generally not correct, because of the nuclear spin;
! 	furthermore, an isomeric product nucleus usually has S>0.
       goto 99                                   !---------------------------------------
!Check conservation of photon energy:		 !===============
       gk1=PP(1)/PP(4)-gk1; gk2=PP(2)/PP(4)-gk2; gk3=PP(3)/PP(4)-gk3
       gk=0.5*PP(4)*(gk1**2+gk2**2+gk3**2)       ! NR kin erg 
       do n=1,m
         if (id(n).eq.0) then
           gk=gk+sqrt(p(1,n)**2+p(2,n)**2+p(3,n)**2)
         endif
       enddo
!      if (mod(kEv,10).eq.0) pause 'DecayS check'!===============
   99  eps=eps0                                  ! Post-evap heat
       RETURN
       END                                       ! DecayS
