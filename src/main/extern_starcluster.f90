!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module extern_starcluster
!
! This module contains routines to model the potential of a stellar cluster
! with a massive object at the centre 
!
! :References: Bonnell et al., 2001, MNRAS, 324, 573–579
!
! :Owner: Cheryl Lau 
!
! :Runtime parameters:
!   - Mcore       : *Total mass within the cluster core*
!   - Mclust      : *Total mass within the outer cluster*
!   - Rcore       : *Radius of cluster core*
!   - Rclust      : *Radius of whole cluster*
!   - dMfrac      : *Step size as fraction of Mcore and Mclust*
!
! :Dependencies: infile_utils, io, part 
!
 implicit none

 real,    public :: Rcore  = 0.1 
 real,    public :: Rclust = 0.5
 real,    public :: Mcore  = 1.2e2
 real,    public :: Mclust = 1e3 

 real,    public :: dMfrac = 5.d-3

 integer, public :: update_mass_freq = 100      ! update Mcore_phi & Mclust_phi every x-th call
 logical, public :: actual_mass_only = .false.  ! only account for mass present in the sim
 logical, public :: vary_potential   = .true.   ! vary potential to keep cluster virialized

 public :: starcluster_force,init_starcluster,update_Mcore_Mclust
 public :: write_options_starcluster,read_options_starcluster

 private

 integer :: icall = 1
 integer :: nzeromass = 0
 real    :: Mcore_phi,Mclust_phi,dMcore_phi,dMclust_phi
 real    :: ekin0,etherm0,epot0,Rcore2,Rclust2
 
contains
!-----------------------------------------------------------------
!+
!  Calculate Mcore_phi and Mclust_phi to be used in the potential 
!+
!-----------------------------------------------------------------
subroutine init_starcluster(ierr)
 use part,  only:npart,xyzh,vxyzu,isdead_or_accreted
 use part,  only:nptmass,xyzmh_ptmass,vxyz_ptmass
 use part,  only:massoftype,npartoftype,igas
 use dim,   only:gravity
 use io,    only:warning 
 integer, intent(out) :: ierr 
 integer :: ip,isink,io_potfile,io_energfile
 real    :: xi,yi,zi,r2
 real    :: Mgas_core,Mgas_clust,Msink_core,Msink_clust 

 ierr = 0

 Rcore2  = Rcore**2 
 Rclust2 = Rclust**2 

 !- Find actual mass contained within Rcore and Rclust 
 Mgas_core   = 0.d0 
 Mgas_clust  = 0.d0
 Msink_core  = 0.d0 
 Msink_clust = 0.d0 
 gas: do ip = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,ip))) then 
       xi = xyzh(1,ip)
       yi = xyzh(2,ip)
       zi = xyzh(3,ip)
       r2 = xi**2 + yi**2 + zi**2 
       if (r2 < Rcore2) then 
          Mgas_core = Mgas_core + massoftype(igas)
       elseif (r2 < Rclust2) then 
          Mgas_clust = Mgas_clust + massoftype(igas)
       else 
          call warning('extern_starcluster','Particle is beyond cluster radius')
          print*, 'particle id: ',ip
       endif 
    endif 
 enddo gas 
 sinks: do isink = 1,nptmass
    if (xyzmh_ptmass(4,isink) > 0) then  !- not merged 
       xi = xyzmh_ptmass(1,isink)
       yi = xyzmh_ptmass(2,isink)
       zi = xyzmh_ptmass(3,isink)
       r2 = xi**2 + yi**2 + zi**2 
       if (r2 < Rcore2) then 
          Msink_core = Msink_core + xyzmh_ptmass(4,isink)
       elseif (r2 < Rclust2) then 
          Msink_clust = Msink_clust + xyzmh_ptmass(4,isink) 
       else 
          call warning('extern_starcluster','Sink is beyond cluster radius')
          print*,'sink id',isink 
       endif 
    endif 
 enddo sinks 
 print*, 'Mass contained in Rcore: ',Mgas_core+Msink_core
 print*, 'Mass contained between Rcore and Rclust: ',Mgas_clust+Msink_clust 

 if (actual_mass_only) then 
    Mcore_phi  = Mgas_core + Msink_core
    Mclust_phi = Mgas_clust + Msink_clust 
 else 
    Mcore_phi  = Mcore 
    Mclust_phi = Mclust 
    if (gravity) then 
       !- If particles are already subjected to self-gravity, the mass terms in the
       !  potential only represent the extra 'hidden' mass which are not explicityly
       !  modelled in the simulation, serve to control the boundness of the gas. 
       print*,'Mcore and Mclust subtracted to account for self-grav'
       Mcore_phi  = Mcore_phi - Mgas_core - Msink_core 
       Mclust_phi = Mclust_phi - Mgas_clust - Msink_clust 
       Mcore_phi  = max(Mcore_phi,0.)
       Mclust_phi = max(Mclust_phi,0.)
    endif    
 endif 
 
 !-File to write Mcore_phi and Mclust_phi 
 print*,'Mcore_phi: ',Mcore_phi,'; Mclust_phi: ',Mclust_phi
 open(2020,file='Mcore_Mclust_evol.dat',status='replace',iostat=io_potfile)
 if (io_potfile /= 0) ierr = 1 
 write(2020,'(3A20)') 'time','Mcore_phi','Mclust_phi'
 write(2020,'(3E20.10)') 0.d0, Mcore_phi, Mclust_phi
 
 !-Init energy mem 
 ekin0   = 0.d0 
 etherm0 = 0.d0 
 epot0   = 0.d0

 !-File to record the computed energies 
 open(2030,file='extern_starcluster_energies.dat',status='replace',iostat=io_energfile)
 if (io_energfile /= 0) ierr = 1 
 write(2030,'(4A20)') 'time','ekin','etherm','epot'

 !-Step in Mcore and Mclust 
 if (vary_potential) then 
    dMcore_phi  = Mcore_phi*dMfrac
    dMclust_phi = Mclust_phi*dMfrac
 endif 

end subroutine init_starcluster 

!-----------------------------------------------------------------
!+
!  Check whether or not Mcore and Mclust needs to be varied
!  to keep the cluster core virialized
!  Called from subroutine update_externalforce every substep
!+
!-----------------------------------------------------------------
subroutine update_Mcore_Mclust(time)
 use io, only:fatal 
 real,    intent(in) :: time
 real    :: ekin,etherm,epot
 real    :: tol = 1.d4
 real    :: alpha_uppthresh = 1.2
 real    :: alpha_lowthresh = 0.8
 logical :: check_now,recalc_Mcore_Mclust

 check_now = .false. 
 if (vary_potential .and. mod(icall,update_mass_freq) == 0) check_now = .true. 

 !- Check if energies have changed since last read 
 recalc_Mcore_Mclust = .false.
 if (check_now) then 
    call compute_energies(ekin,etherm,epot)
    write(2030,'(4E20.10)') time, ekin, etherm, epot 
 
    if (abs(ekin-ekin0) > tol .or. abs(etherm-etherm0) > tol .or. abs(epot-epot0) > tol) then 
       recalc_Mcore_Mclust = .true. 
       ekin0   = ekin
       etherm0 = etherm 
       epot0   = epot 
    endif 
 endif 

 if (recalc_Mcore_Mclust) then 
    if (2.d0*ekin/abs(epot) > alpha_uppthresh) then      !- unbound 
       Mcore_phi  = Mcore_phi + dMcore_phi
       Mclust_phi = Mclust_phi + dMclust_phi
    elseif (2.d0*ekin/abs(epot) < alpha_lowthresh) then  !- bound 
       Mcore_phi  = Mcore_phi - dMcore_phi
       Mclust_phi = Mclust_phi - dMclust_phi
       Mcore_phi  = max(Mcore_phi,0.)
       Mclust_phi = max(Mclust_phi,0.)

       if (Mclust_phi <= 0 .or. Mcore_phi <= 0) then 
          nzeromass = nzeromass + 1 
          if (nzeromass > 50) call fatal('extern_starcluster','Particles are doomed to collapse')
       endif 
    endif 
    write(2020,'(3E20.10)') time, Mcore_phi, Mclust_phi
 endif 

 icall = icall + 1 

end subroutine update_Mcore_Mclust

!
! Compute kinetic, thermal, potential energy 
!
subroutine compute_energies(ekin,etherm,epot)
 use dim,    only:maxvxyzu,gravity
 use part,   only:npart,xyzh,vxyzu,massoftype,igas,poten,nptmass,xyzmh_ptmass
 use part,   only:epot_sinksink,isdead_or_accreted
 real,   intent(out) :: ekin,etherm,epot 
 integer :: ip 
 real    :: pmass,xi,yi,zi,hi,vxi,vyi,vzi,v2i,phi,epot_sinks,ekin_sinks 

 pmass = massoftype(igas)

 ekin = 0.d0 
 etherm = 0.d0 
 epot = 0.d0

 do ip = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,ip))) then 
       xi  = xyzh(1,ip)
       yi  = xyzh(2,ip)
       zi  = xyzh(3,ip)
       hi  = xyzh(4,ip)

       !--Kinetic energy 
       vxi = vxyzu(1,ip)
       vyi = vxyzu(2,ip)
       vzi = vxyzu(3,ip)
       v2i = vxi*vxi + vyi*vyi + vzi*vzi
       ekin = ekin + 5.d-1*pmass*v2i

       !--Thermal energy 
       if (maxvxyzu >= 4) etherm = etherm + vxyzu(4,ip)*pmass 

       !--Potential energy 
       call cluster_potential(xi,yi,zi,phi)
       epot = epot + phi*pmass 
       if (gravity) then 
          epot = epot + poten(ip)
       endif 
       if (nptmass > 0) then
          call get_accel_from_sinks(xi,yi,zi,hi,epot_sinks)
          epot = epot + pmass*epot_sinks 
       endif
    endif 
 enddo 

 !--Kinetic energy of sinks 
 if (nptmass > 0) then 
    call get_ekin_of_sinks(ekin_sinks)
    ekin = ekin + ekin_sinks 
 endif 

 !--Potential energy of sinks 
 if (nptmass > 1) then 
    epot = epot + epot_sinksink
 endif 

end subroutine compute_energies 

!
! Cluster potential that corresponds to density distribution 
! \rho \propto r^{-3/2};   r < Rcore 
!              r^{-2};     Rcore < r < Rclust 
!
subroutine cluster_potential(xi,yi,zi,phi)
 use io, only:fatal 
 real,   intent(in)  :: xi,yi,zi
 real,   intent(out) :: phi 
 real   :: r2,r

 r2  = xi**2 + yi**2 + zi**2 

 if (r2 < Rcore2) then 
    r   = sqrt(r2)
    phi = -3.*Mcore_phi*Rcore**(-1.5) * (Rcore**(0.5)- 2./3.*r**(0.5)) &
        & - Mclust_phi * (log(Rclust)-log(Rcore))/(Rclust-Rcore)
 elseif (r2 < Rclust2) then 
    r   = sqrt(r2)
    phi = -((Mcore_phi - Mclust_phi*Rcore/(Rclust-Rcore))*1./r + Mclust_phi/(Rclust-Rcore)) &
        & - Mclust_phi/(Rclust-Rcore) * (log(Rclust)-log(r))
 else 
    phi = -1.*tiny(phi)
 endif 

end subroutine cluster_potential


!
! Kinetic energy of sinks 
!
subroutine get_ekin_of_sinks(ekin_sinks)
 use part, only:nptmass,xyzmh_ptmass,vxyz_ptmass
 real, intent(out) :: ekin_sinks 
 integer :: isink 
 real    :: msink,vx,vy,vz,v2

 ekin_sinks = 0. 

 do isink = 1,nptmass 
    msink = xyzmh_ptmass(4,isink)
    if (msink < 0.) cycle 
    vx = vxyz_ptmass(1,isink)
    vy = vxyz_ptmass(2,isink)
    vz = vxyz_ptmass(3,isink)
    v2 = vx**2 + vy**2 + vz**2 
    ekin_sinks = ekin_sinks + 5.d-1*msink*v2 
 enddo 

end subroutine get_ekin_of_sinks 

!
! Copy of subroutine get_accel_sink_gas from module ptmass 
! but only for getting the potential exerted by sinks 
!
subroutine get_accel_from_sinks(xi,yi,zi,hi,epot_sinks)
 use part,     only:nptmass,xyzmh_ptmass
 use kernel,   only:kernel_softening,radkern
 real, intent(in)  :: xi,yi,zi,hi
 real, intent(out) :: epot_sinks 
 integer :: isink 
 real    :: dx,dy,dz,r,r2,msink,hsoft,qi,q2i,psoft,fsoft

 epot_sinks = 0. 

 do isink = 1,nptmass
    msink  = xyzmh_ptmass(4,isink)
    if (msink < 0.0) cycle

    dx = xi - xyzmh_ptmass(1,isink)
    dy = yi - xyzmh_ptmass(2,isink)
    dz = zi - xyzmh_ptmass(3,isink)
    r2 = dx**2 + dy**2 + dz**2 

    hsoft  = xyzmh_ptmass(6,isink)
    if (hsoft > 0.0) hsoft = max(hsoft,hi)

    if (r2 < (radkern*hsoft)**2) then  !-within softening length 
       q2i = r2/hsoft**2 
       qi  = sqrt(q2i)
       call kernel_softening(q2i,qi,psoft,fsoft) 
       epot_sinks = epot_sinks + msink*psoft/hsoft  ! potential (spline-softened)

    else !-no softening needed 
       r = sqrt(r2)
       epot_sinks = epot_sinks - msink/r   ! potential (GM/r)
    endif 
 enddo 

end subroutine get_accel_from_sinks

!-----------------------------------------------------------------
!+
!  compute the force on a given particle within a stellar cluster 
!+
!-----------------------------------------------------------------
subroutine starcluster_force(xi,yi,zi,fxi,fyi,fzi,phi)
 real,    intent(in)  :: xi,yi,zi
 real,    intent(out) :: fxi,fyi,fzi,phi
 real    :: r2,r
 real    :: eps2_soft = 1.d-1

 r2 = xi*xi + yi*yi + zi*zi + eps2_soft
 
 if (r2 < Rcore2) then 
    r = sqrt(r2)
    fxi = -Mcore_phi*Rcore**(-1.5d0) * xi*r**(-1.5d0)
    fyi = -Mcore_phi*Rcore**(-1.5d0) * yi*r**(-1.5d0)
    fzi = -Mcore_phi*Rcore**(-1.5d0) * zi*r**(-1.5d0)
 elseif (r2 < Rclust2) then 
    r = sqrt(r2)
    fxi = -(Mcore_phi + Mclust_phi*(r-Rcore)/(Rclust-Rcore)) * xi*r**(-3.d0)
    fyi = -(Mcore_phi + Mclust_phi*(r-Rcore)/(Rclust-Rcore)) * yi*r**(-3.d0)
    fzi = -(Mcore_phi + Mclust_phi*(r-Rcore)/(Rclust-Rcore)) * zi*r**(-3.d0)
 else 
    fxi = 0.d0
    fyi = 0.d0
    fzi = 0.d0
 endif 

 !--Get phi 
 call cluster_potential(xi,yi,zi,phi)

end subroutine starcluster_force


!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_starcluster(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(Mcore,'Mcore','Mass in cluster core',iunit)
 call write_inopt(Mclust,'Mclust','Mass in outer cluster',iunit)
 call write_inopt(Rcore,'Rcore','Radius of cluster core',iunit)
 call write_inopt(Rclust,'Rclust','Radius of whole cluster',iunit)
 call write_inopt(dMfrac,'dMfrac','Step size as fraction of Mcore/Mclust',iunit)
 call write_inopt(vary_potential,'vary_potential','Vary potential to keep cluster virialized',iunit)

end subroutine write_options_starcluster

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_starcluster(name,valstring,imatch,igotall,ierr)
 use physcon, only: pi
 use io,      only:fatal,error
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: where = 'read_options_externstarcluster'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('Mcore')
    read(valstring,*,iostat=ierr) Mcore
    ngot = ngot+1
 case('Mclust')
    read(valstring,*,iostat=ierr) Mclust
    ngot = ngot+1
 case('Rcore')
    read(valstring,*,iostat=ierr) Rcore
    ngot = ngot+1
 case('Rclust')
    read(valstring,*,iostat=ierr) Rclust
    ngot = ngot+1
 case('dMfrac')
    read(valstring,*,iostat=ierr) dMfrac
    ngot = ngot+1
 case('vary_potential')
    read(valstring,*,iostat=ierr) vary_potential
    ngot = ngot+1
 end select

 igotall = (ngot >= 6)

end subroutine read_options_starcluster

end module extern_starcluster
