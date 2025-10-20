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
! :References: 
!
! :Owner: Cheryl Lau 
!
! :Runtime parameters:
!   - Mclust            : *Total mass within the outer cluster*
!   - Rcore             : *Radius of cluster core*
!   - dMfrac            : *Step size as fraction of Mclust*
!   - vary_potential    : *Option to vary Mclust to keep cloud virialized*
!
! :Dependencies: infile_utils, io, part, dim, physcon, 
!
 implicit none

 real,    public :: Rcore  = 0.1 
 real,    public :: Mclust = 1e3 

 real,    public :: dMfrac = 5.d-3

 integer, public :: update_mass_freq = 100      ! update Mclust_phi every x-th call
 logical, public :: actual_mass_only = .false.  ! only account for mass present in the sim
 logical, public :: vary_potential   = .true.   ! vary potential to keep cluster virialized

 public :: starcluster_force,init_starcluster,update_Mclust
 public :: write_options_starcluster,read_options_starcluster

 private

 integer :: icall = 1
 integer :: nzeromass = 0
 real    :: Mclust_phi,dMclust_phi
 real    :: ekin0,etherm0,epot0
 
 !--Storage for cluster profiles 
 integer, parameter :: nRmax = 100000
 integer :: nR 
 real    :: r_profile(nRmax),rho_profile(nRmax)
 real    :: phi_profile(nRmax),force_profile(nRmax)

 !--Flags for checking and debugging 
 logical :: print_Mclust  = .true.  
 logical :: print_energy  = .true. 
 logical :: print_profile = .true. 
 logical :: use_current_sigma = .false. 

contains
!-----------------------------------------------------------------
!+
!  Calculate initial Mclust_phi to be used in the potential 
!+
!-----------------------------------------------------------------
subroutine init_starcluster(ierr)
 use part,  only:npart,nptmass,xyzmh_ptmass,vxyz_ptmass
 use part,  only:massoftype,npartoftype,igas
 use dim,   only:gravity
 use io,    only:warning,fatal 
 integer, intent(out) :: ierr 
 integer :: ip,isink,io_potfile,io_energfile
 real    :: sigma,j,j2,phi0,W0,Rclust
 real    :: Mgas_clust,Msink_clust 

 ierr = 0

 !--Find total mass present in the simulation 
 Mgas_clust  = npart * massoftype(igas)
 Msink_clust = 0.
 sinks: do isink = 1,nptmass
    if (xyzmh_ptmass(4,isink) > 0) then  !- not merged 
       Msink_clust = Msink_clust + xyzmh_ptmass(4,isink) 
    endif 
 enddo sinks 
 print*, 'Mass contained in simulation: ',Mgas_clust+Msink_clust 

 if (actual_mass_only) then 
    Mclust_phi = Mgas_clust + Msink_clust 
 else 
    Mclust_phi = Mclust 
    if (gravity) then 
       !- If particles are already subjected to self-gravity, the mass terms in the
       !  potential only represent the extra 'hidden' mass which are not explicityly
       !  modelled in the simulation, serve to control the boundness of the gas. 
       print*,'Mclust subtracted to account for self-grav'
       Mclust_phi = Mclust_phi - Mgas_clust - Msink_clust 
       Mclust_phi = max(Mclust_phi,0.)
    endif    
 endif 
 

 !-File to write Mclust_phi, phi0, W0, sigma 
 if (print_Mclust) then 
    open(2020,file='Mclust_phi0_evol.dat',status='replace',iostat=io_potfile)
    if (io_potfile /= 0) ierr = 1 
    write(2020,'(5A20)') 'time','Mclust_phi','phi0','W0','sigma'
 endif 

 !--Compute cluster potential for the first time 
 call cluster_profile(phi0,W0,Rclust,sigma)
 if (print_Mclust) write(2020,'(5E20.10)') 0.d0, Mclust_phi, phi0, W0, sigma 

! call fatal('extern_starcluster','force stop')

 !-Init energy mem 
 ekin0   = 0.d0 
 etherm0 = 0.d0 
 epot0   = 0.d0

 !-File to record the computed energies 
 if (print_energy) then 
    open(2030,file='extern_starcluster_energies.dat',status='replace',iostat=io_energfile)
    if (io_energfile /= 0) ierr = 1 
    write(2030,'(4A20)') 'time','ekin','etherm','epot'
 endif 

 !-Step in Mclust 
 if (vary_potential) then 
    dMclust_phi = Mclust_phi*dMfrac
 endif 

end subroutine init_starcluster 

!-----------------------------------------------------------------
!+
!  Check whether or not Mclust needs to be varied
!  to keep the cluster core virialized
!  Called from subroutine update_externalforce every substep
!+
!-----------------------------------------------------------------
subroutine update_Mclust(time)
 use io, only:fatal 
 real,    intent(in) :: time
 real    :: ekin,etherm,epot,phi0,W0,Rclust,sigma
 real    :: tol = 1.d4
 real    :: alpha_uppthresh = 1.2
 real    :: alpha_lowthresh = 0.8
 logical :: check_now,recalc_Mclust

 check_now = .false. 
 if (vary_potential .and. mod(icall,update_mass_freq) == 0) check_now = .true. 

 !- Check if energies have changed since last read 
 recalc_Mclust = .false.
 if (check_now) then 
    call compute_energies(ekin,etherm,epot)
    if (print_energy) write(2030,'(4E20.10)') time, ekin, etherm, epot 
 
    if (abs(ekin-ekin0) > tol .or. abs(etherm-etherm0) > tol .or. abs(epot-epot0) > tol) then 
       recalc_Mclust = .true. 
       ekin0   = ekin
       etherm0 = etherm 
       epot0   = epot 
    endif 
 endif 

 if (recalc_Mclust) then 
    if (2.d0*ekin/abs(epot) > alpha_uppthresh) then      !- unbound 
       Mclust_phi = Mclust_phi + dMclust_phi
    elseif (2.d0*ekin/abs(epot) < alpha_lowthresh) then  !- bound 
       Mclust_phi = Mclust_phi - dMclust_phi
       Mclust_phi = max(Mclust_phi,0.)

       if (Mclust_phi <= 0.) then 
          nzeromass = nzeromass + 1 
          if (nzeromass > 50) call fatal('extern_starcluster','Particles are doomed to collapse')
       endif 
    else                                                 !- virialized
       use_current_sigma = .true.  !!!! TESTING !!!!
    endif 

    !--Recompute cluster profile with current Mclust_phi & sigma  
    !  and tabulate for later access (stored globally)
    call cluster_profile(phi0,W0,Rclust,sigma)

    if (print_Mclust) write(2020,'(5E20.10)') time, Mclust_phi, phi0, W0, sigma 
 endif 

 icall = icall + 1 

end subroutine update_Mclust

!
! Compute kinetic, thermal, potential energy 
!
subroutine compute_energies(ekin,etherm,epot)
 use dim,    only:maxvxyzu,gravity
 use part,   only:npart,xyzh,vxyzu,massoftype,igas,poten,nptmass,xyzmh_ptmass
 use part,   only:epot_sinksink,isdead_or_accreted
 real,   intent(out) :: ekin,etherm,epot 
 integer :: ip 
 real    :: pmass,xi,yi,zi,hi,vxi,vyi,vzi,v2i,fxi,fyi,fzi,phi,epot_sinks,ekin_sinks 

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
       call starcluster_force(xi,yi,zi,fxi,fyi,fzi,phi)
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
! With current Mclust_phi, computes cluster profile with King models 
! Produces tabulated density, potential and force as functions of R
!+
!-----------------------------------------------------------------
subroutine cluster_profile(phi0,W0,Rclust,sigma)
 use part,   only:npart,xyzh,vxyzu
 use units,  only:unit_density,utime,udist,umass 
 use io,     only:fatal 
 real,   intent(out) :: phi0,W0,Rclust,sigma
 integer :: iR,io_clusterfile
 real    :: j,j2,k,ve2,R,dWdR,d2WdR2,W,dR,rhomin,phi,rho,dphidr,force
 real    :: dRmax_dW,dRmax_dR
 real    :: r_pc,rho_cgs,phi_cgs,force_cgs 

 !--Use current Mclust_phi to compute W0 with Plummer model
 call get_vel_dispersion(npart,xyzh,vxyzu,sigma,j,j2)
 phi0 = -Mclust_phi/Rcore 
 W0   = -2.d0*j2*phi0

 !--Estimate k 
 ve2  = -2.d0*phi0 
 k    = (1.d0 - exp(-1.d0*j2*ve2))**(-1)

 !--Initialize 
 R    = 1.d-3  ! close to 0 
 dWdR = 1.d-10
 W    = W0 
 dR   = 1.d-3

 rhomin = 1.d-27/unit_density
 phi = -1   ! dummy 
 rho = 10.  ! dummy 
 iR  = 0

 scan_over_R: do while (rho > rhomin .and. phi < 0.)

    rho = rho_as_func_of_W(W,W0,k,j)

    d2WdR2 = d2WdR2_poisson(dWdR,R,rho,j2)

    phi = W/(-2.d0*j2)
    phi = min(phi,0.d0)

    dphidr = dWdR / (-2.d0*j2*Rcore) 
    force  = -1.d0*dphidr 
    force  = min(force,0.d0)

    !--Constrain dR and update 
    dRmax_dR = 1.d-1 * R
    dRmax_dW = 1.d-3 * abs(W/dWdR)
    dR = min(dR,dRmax_dR,dRmax_dW)
    R = R + dR

    !--Recalc for next iteration 
    dWdR = dWdR + d2WdR2 * dR 
    W = W + dWdR * dR

    !--Store results in code units 
    iR = iR + 1 
    if (iR > nRmax) call fatal('extern_starcluster','number of R entries exceeded limit')
    r_profile(iR)       = R*Rcore
    rho_profile(iR)     = rho
    phi_profile(iR)     = phi
    force_profile(iR)   = force 

 enddo scan_over_R
 
 !--Current number of entries in profile 
 nR = iR 

 !--Truncation radius 
 Rclust = R*Rcore

 !--Write profile to file for checking 
 if (print_profile) then 
    open(2050,file='cluster_profile.dat',status='replace',iostat=io_clusterfile)
    if (io_clusterfile /= 0) call fatal('extern_starcluster','error opening cluster profile file')
    write(2050,'(4A20)') 'r [pc]','rho [g/cm3]','potential [cm2/s2]','force [cm/s2]'
    do iR = 1,nR
       r_pc       = r_profile(iR)
       rho_cgs    = rho_profile(iR)*unit_density
       phi_cgs    = phi_profile(iR)*udist**2/utime**2
       force_cgs  = force_profile(iR)*umass*udist/utime**2
       write(2050,'(4E20.10)') r_pc, rho_cgs, phi_cgs, force_cgs 
    enddo 
    close(2050)
 endif 

 !--Remove particles beyond truncation radius 
 call truncate_cloud(Rclust,npart,xyzh)

end subroutine cluster_profile


!
! Computes current sigma, j, and j^2 
!
subroutine get_vel_dispersion(npart,xyzh,vxyzu,sigma,j,j2)
 use part,  only:isdead_or_accreted
 use units, only:unit_velocity 
 integer, intent(in)  :: npart 
 real,    intent(in)  :: xyzh(:,:),vxyzu(:,:)
 real,    intent(out) :: sigma,j,j2
 integer :: np_alive,i
 real    :: vsum,v2sum,v2,v,var 

 vsum     = 0.
 v2sum    = 0. 
 np_alive = 0 

 do i = 1,npart 
    if (.not.isdead_or_accreted(xyzh(4,i))) then 
       v2 = vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2
       v  = sqrt(v2)
       vsum  = vsum  + v
       v2sum = v2sum + v2 
       np_alive = np_alive + 1 
    endif 
 enddo 
 var   = (v2sum - vsum**2/np_alive) / (np_alive-1)
 sigma = sqrt(var)

 !--testing 
 if (.not.use_current_sigma) sigma = 2.3e5/unit_velocity

 j2 = 1.d0/(2*sigma**2)
 j  = sqrt(j2)

end subroutine get_vel_dispersion

!
! Expression for rho(W) derived from distribution function 
!
real function rho_as_func_of_W(W,W0,k,j)
 use physcon, only:pi
 real, intent(in) :: W,W0,k,j
 real :: integral 

 !print*,'W,W0,k,j in func',W,W0,k,j
 integral = simpson_rule_equi(rhoW_integrand,0.d0,W,1000)
 !print*,'integral',integral 
 rho_as_func_of_W = 4.d0/3.d0*pi*k*j**(-3)*exp(W-W0)*integral

end function rho_as_func_of_W

!
! Integrand in the expression of rho(W)
!
real function rhoW_integrand(eta)
 real, intent(in) :: eta

 rhoW_integrand = exp(-1.d0*eta) * eta**(3.d0/2.d0)

end function rhoW_integrand

!
! Simpson-rule with equidistant spacing
!
real function simpson_rule_equi(func,a,b,ninterval)
 integer, intent(in) :: ninterval    ! number of sub-intervals
 real,    intent(in) :: a,b          ! boundary values
 real,    external   :: func         ! the function to be integrated
 integer :: i
 real    :: dx,x1,x2,xm,f1,f2,fm,intsum

 dx  = (b-a)/dble(ninterval)
 x1  = a							! left
 f1  = func(a)
 intsum = 0.d0
 do i = 1,ninterval
    x2  = a + dble(i)*dx      ! right
    xm  = 0.5d0*(x1+x2)       ! midpoint
    f2  = func(x2)
    fm  = func(xm)
    intsum = intsum + (f1+4.d0*fm+f2)/6.d0*(x2-x1)  ! Simpson rule
    x1  = x2
    f1  = f2                                 ! save for next subinterval
 enddo
 simpson_rule_equi = intsum

end function simpson_rule_equi

!
! Expression for d^2W/dR^2 derived from Poisson equation 
!
real function d2WdR2_poisson(dWdR,R,rhoW,j2)
 use physcon, only:pi
 real, intent(in) :: dWdR,R,rhoW,j2

 d2WdR2_poisson = -8.d0*pi*j2*Rcore**2*rhoW - 2.d0/R*dWdR 

end function d2WdR2_poisson


!-----------------------------------------------------------------
!+
!  compute the force/potential on a given particle 
!+
!-----------------------------------------------------------------
subroutine starcluster_force(xi,yi,zi,fxi,fyi,fzi,phi)
 real,    intent(in)  :: xi,yi,zi
 real,    intent(out) :: fxi,fyi,fzi,phi
 real    :: r2i,ri,rho,fr,theta_angle,phi_angle 

 r2i = xi**2 + yi**2 + zi**2 
 ri  = sqrt(r2i)

 rho = interp_from_profile(ri,nR,r_profile,rho_profile)
 phi = interp_from_profile(ri,nR,r_profile,phi_profile)
 fr  = interp_from_profile(ri,nR,r_profile,force_profile)

 theta_angle = acos(zi/ri)
 phi_angle   = atan(yi/xi)
 fxi = fr*sin(theta_angle)*cos(phi_angle)
 fyi = fr*sin(theta_angle)*sin(phi_angle)
 fzi = fr*cos(theta_angle)

end subroutine starcluster_force

!
! Function to extract/interpolate from a tabulated cluster profile
!
real function interp_from_profile(ri,nR,rad_profile,A_profile)
 integer, intent(in) :: nR 
 real,    intent(in) :: ri
 real,    intent(in) :: rad_profile(nR),A_profile(nR) 
 integer :: i,j

 !--Locate the closest entries and bracket ri 
 ! ([i] closest index; [j] lower bound; [j+1] upper bound around ri)
 i = minloc(abs(rad_profile(:)-ri),1)
 j = 0
 if (i == 1) then
    j = 1
 elseif (i == nR) then
    j = i-1
 elseif (rad_profile(i) >= ri .and. rad_profile(i-1) < ri) then
    j = i-1
 elseif (rad_profile(i) < ri .and. rad_profile(i+1) >= ri) then
    j = i
 endif
 
 !--Interpolate 
 interp_from_profile = A_profile(j) + (ri-rad_profile(j))*(A_profile(j+1)  &
                     & - A_profile(j))/(rad_profile(j+1)-rad_profile(j))

end function interp_from_profile


!-----------------------------------------------------------------
!+
!  Remove particles beyond the cluster truncated radius 
!  We assume they are lost to tidal forces from galactic centre 
!+
!-----------------------------------------------------------------
subroutine truncate_cloud(Rclust,npart,xyzh)
 use part, only:kill_particle 
 integer, intent(inout) :: npart 
 real,    intent(inout) :: Rclust 
 real,    intent(inout) :: xyzh(:,:)
 integer :: ip
 real    :: Rclust2,xi,yi,zi,r2i
 
 Rclust2 = Rclust**2

 do ip = 1,npart 
    xi = xyzh(1,ip)
    yi = xyzh(2,ip)
    zi = xyzh(3,ip)
    r2i = xi**2 + yi**2 + zi**2 
    if (r2i > Rclust2) then 
       print*,'particle ',ip,' gone beyond cluster truncation radius'
       xyzh(4,ip) = 0.d0 
    endif 
 enddo 

end subroutine truncate_cloud


!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_starcluster(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(Mclust,'Mclust','Mass in outer cluster',iunit)
 call write_inopt(Rcore,'Rcore','Radius of cluster core',iunit)
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
 case('Mclust')
    read(valstring,*,iostat=ierr) Mclust
    ngot = ngot+1
 case('Rcore')
    read(valstring,*,iostat=ierr) Rcore
    ngot = ngot+1
 case('dMfrac')
    read(valstring,*,iostat=ierr) dMfrac
    ngot = ngot+1
 case('vary_potential')
    read(valstring,*,iostat=ierr) vary_potential
    ngot = ngot+1
 end select

 igotall = (ngot >= 4)

end subroutine read_options_starcluster

end module extern_starcluster
