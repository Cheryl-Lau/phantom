!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which gets the mean of physical properties in a given box
! i.e. spreading everything over the box. 
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies:
!

 implicit none
 character(len=20), parameter, public :: analysistype = 'mean_box'

 public :: do_analysis

 private
 real    :: centre(3) = (/ 0., 0., 0. /)
 real    :: radius = 10. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,pmass,npart,time,iunit)
 use units,   only:unit_velocity,unit_density,umass,unit_ergg 
 use physcon, only:kboltz,mass_proton_cgs,solarm
 use eos,     only:gmw,gamma 
 use io,      only:fatal
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: pmass,time
 integer :: i,ipart_inbox
 real    :: totmass,vol_box,vel_mean,u_mean,rho_mean,ekin,etherm
 real    :: m_cgs,rho_cgs,vel_cgs,u_cgs,temp 
 character(len=70) :: filename

 vol_box = (2.*radius)**3 
 
 !- Init counters 
 totmass  = 0.
 ekin = 0.
 etherm   = 0. 
 ipart_inbox = 0

 do i = 1,npart
    !
    ! Filter those beyond box
    !
    if (xyzh(1,i) < centre(1)-radius .or. xyzh(1,i) > centre(1)+radius .or. &
        xyzh(2,i) < centre(2)-radius .or. xyzh(2,i) > centre(2)+radius .or. &
        xyzh(3,i) < centre(3)-radius .or. xyzh(3,i) > centre(3)+radius) cycle 
    !
    ! Calculate ekin
    !
    ekin = ekin + 5d-1*pmass*(vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2)
    !
    ! Calculate etherm
    !
    etherm = etherm + vxyzu(4,i)*pmass
    !
    ! Calculate mass
    !
    totmass = totmass + pmass
    !
    ! Count number 
    !
    ipart_inbox = ipart_inbox + 1 
 enddo 

 print*,'number of particles in box: ', ipart_inbox,'/',npart 

 !
 ! Average density
 !
 rho_mean = totmass/vol_box
 !
 ! Average velocity 
 !
 vel_mean = sqrt(2.d0*ekin/totmass)
 !
 ! Average u 
 !
 u_mean = etherm/totmass 

 !
 ! Convert to cgs units 
 ! 
 m_cgs   = totmass*umass 
 rho_cgs = rho_mean*unit_density
 vel_cgs = vel_mean*unit_velocity 
 u_cgs   = u_mean*unit_ergg
 temp    = u_cgs/kboltz*(gmw*mass_proton_cgs*(gamma-1.))


 filename = 'box_mean_'//TRIM(dumpfile)//'.dat'
 open(unit=2024,file=filename)
 write(2024,'(5a20)') 'totmass (m_sun)','rho_mean','vel_mean','u_mean','temp_mean' 
 write(2024,'(5es20.10)') m_cgs/solarm, rho_cgs, vel_cgs, u_cgs, temp 
 close(2024)


end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
