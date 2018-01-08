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
       module freyainterfaces
!
!      interfaces between FORTRAN and C for bindings
!      ---------------------------------------------
!
       interface
         subroutine msfreya_event_c(iK,Einc,eps0,PP0,iZ1,iA1, &
          PP1,iZ2,iA2,PP2,m,p,id,ndir) bind (C, name="msfreya_event_c_")
           use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
           use freyaparameters
           implicit none
           integer (kind=c_int), value :: iK
           real (kind=c_double), value :: Einc,eps0
           ! Four-momentum of the fissioning nucleus
           real (kind=c_double), dimension(0:4) :: PP0
           ! Four-momentum of the 1st evaporating nucleus
           real (kind=c_double), dimension(0:4) :: PP1
           ! Four-momentum of the 2nd evaporating nucleus
           real (kind=c_double), dimension(0:4) :: PP2
           ! Direction of the incident neutron (normalized)
           real (kind=c_double), dimension(1:3) :: ndir
           integer (kind=c_int) :: iZ1,iA1,iZ2,iA2
           integer (kind=c_int) :: m
           ! Ejectile momenta of all photons and neutrons
           real (kind=c_double), dimension(4,3*mMax) :: p
           ! Ejectile type 
           integer (kind=c_int), dimension(3*mMax) :: id
         end subroutine msfreya_event_c
       end interface
!
       interface
         subroutine msfreya_getids_c(id0_c,id1_c,id2_c) &
           bind (C, name="msfreya_getids_c_")
           use, intrinsic :: iso_c_binding, only: C_INT
           use freyaparameters
           implicit none
           ! Ejectile type of prefission neutrons
           integer (kind=c_int), dimension(3*mMax) :: id0_c
           ! Ejectile type of 1st fission fragment
           integer (kind=c_int), dimension(3*mMax) :: id1_c
           ! Ejectile type of 2nd fission fragment
           integer (kind=c_int), dimension(3*mMax) :: id2_c
         end subroutine msfreya_getids_c
       end interface
!
       interface
         function msFREYA_SEPn_c (iK,iZ,iA)  &
         bind (C, name="msfreya_sepn_c_")
           use, intrinsic :: iso_c_binding
           implicit none
           real (kind=c_double) :: msFREYA_SEPn_c
           integer (kind=c_int), value :: iK,iZ,iA
         end function msFREYA_SEPn_c
       end interface
!
       interface
         subroutine msfreya_setpath_c (pathstring,strlength) &
         bind (C, name="msfreya_setpath_c_")
           use freyaPath, only: maxDATAPATH
           use, intrinsic :: iso_c_binding
           implicit none
           integer(kind=c_int), value :: strlength
           character (kind=c_char,len=1), &
             dimension(maxDATAPATH) :: pathstring
         end subroutine msfreya_setpath_c
       end interface
!
       interface
         subroutine msfreya_setup_c () &
         bind (C, name="msfreya_setup_c_")
           use, intrinsic :: iso_c_binding
           implicit none
         end subroutine msfreya_setup_c
       end interface
!
       interface
         subroutine msfreya_getniso_c(nisosf,nisoif) &
         bind (C, name="msfreya_getniso_c_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int) :: nisosf,nisoif
         end subroutine msfreya_getniso_c
       end interface
!
       interface
         subroutine msfreya_getzas_c(zas,fistypes) &
         bind (C, name="msfreya_getzas_c_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int), dimension(*) :: zas,fistypes
         end subroutine msfreya_getzas_c
       end interface
!
       interface
         function msFREYA_gsMassN_c (iZ,iA)  &
         bind (C, name="msfreya_gsmassn_c_")
           use, intrinsic :: iso_c_binding
           implicit none
           real (kind=c_double) :: msFREYA_gsMassN_c
           integer (kind=c_int), value :: iZ,iA
         end function msFREYA_gsMassN_c
       end interface
!
       interface
         subroutine msFREYA_usehostrng_c() &
         bind(C, name="msfreya_usehostrng_c_")
           use, intrinsic :: iso_c_binding
           implicit none
         end subroutine msFREYA_usehostrng_c
       end interface
!
       interface
         subroutine msFREYA_getlasterror_c(message,counter) &
         bind(C, name="msfreya_getlasterror_c_")
           use, intrinsic :: iso_c_binding
           implicit none
           character (kind=c_char,len=1), dimension(256) :: message
           integer (kind=c_int) :: counter
         end subroutine msFREYA_getlasterror_c
       end interface
!
       interface
         subroutine msFREYA_geterrors_c(messages,length) &
         bind(C, name="msfreya_geterrors_c_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int) :: length
           character (kind=c_char,len=1), dimension(length) :: messages
         end subroutine msFREYA_geterrors_c
       end interface
!
       interface
         subroutine msFREYA_reseterrorflag_c() &
         bind(C, name="msfreya_reseterrorflag_c_")
           use, intrinsic :: iso_c_binding
           implicit none
         end subroutine msFREYA_reseterrorflag_c
       end interface
!
       interface
         function msFREYA_errorflagset_c ()  &
         bind (C, name="msfreya_errorflagset_c_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int) :: msFREYA_errorflagset_c
         end function msFREYA_errorflagset_c
       end interface
!
       interface
         function rngdptr ()  &
         bind (C, name="rngdptr_")
           use, intrinsic :: iso_c_binding
           implicit none
           real (kind=c_double) :: rngdptr
         end function rngdptr
       end interface
!
       end module freyainterfaces
