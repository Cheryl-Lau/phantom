!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which bins the particles by velocity
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
 character(len=20), parameter, public :: analysistype = 'bin_velocity'

 public :: do_analysis

 private
 integer, parameter :: num_velbin = 500
 integer :: binned_vel(num_velbin)
 real    :: velbin_min_cgs = 1.5E+3
 real    :: velbin_max_cgs = 4E+7
 real    :: logvelbin_min,logvelbin_max
 logical :: locate_range = .true. 
 real    :: centre(3) = (/ 0.,0.,0. /)
 real    :: radius = 70. 
 logical :: box_only = .true. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units, only:unit_velocity
 use io,    only:fatal
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,ivelbin
 real    :: vel,dvel,velbin(num_velbin),maxvel,minvel
 character(len=70) :: filename

 logvelbin_min = log10(velbin_min_cgs/unit_velocity)
 logvelbin_max = log10(velbin_max_cgs/unit_velocity)

 !
 ! Separation of rho in log-scale
 !
 dvel = (logvelbin_max-logvelbin_min)/num_velbin
 !
 ! Init binned array
 !
 do ivelbin = 1,num_velbin
    binned_vel(ivelbin) = 0.
 enddo
 !
 ! Bin all particles by density
 !
 minvel = huge(minvel)
 maxvel = tiny(maxvel)
 do i = 1,npart
    !
    ! Filter 
    !
    if (box_only) then 
       if (xyzh(1,i) < centre(1)-radius .or. xyzh(1,i) > centre(1)+radius .or. &
           xyzh(2,i) < centre(2)-radius .or. xyzh(2,i) > centre(2)+radius .or. &
           xyzh(3,i) < centre(3)-radius .or. xyzh(3,i) > centre(3)+radius) cycle 
    endif 
    !
    ! Calculate |v|
    !
    vel = sqrt(vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2)
    if (.not.locate_range) then 
       !
       ! Bin velocity
       !
       ivelbin = (log10(vel)-logvelbin_min+dvel)/dvel
       if (ivelbin < 0) then 
          print*,vel*unit_velocity
          call fatal('analysis_binvelocity','require smaller logvelbin_min')
       elseif (ivelbin > num_velbin) then 
          print*,vel*unit_velocity
          call fatal('analysis_binvelocity','require larger logvelbin_max')
       endif 
       binned_vel(ivelbin) = binned_vel(ivelbin) + 1
    else 
       minvel = min(vel,minvel)
       maxvel = max(vel,maxvel)
    endif 
 enddo
 if (locate_range) then 
    print*,'max',maxvel*unit_velocity,'min',minvel*unit_velocity 
 else 
    !
    ! Convert to physical units and store results
    !
    do ivelbin = 1,num_velbin
       velbin(ivelbin) = (10**(logvelbin_min + (ivelbin-1)*dvel))*unit_velocity
    enddo

    filename = 'velocity_binned_'//TRIM(dumpfile)//'.dat'
    open(unit=2024,file=filename)
    do ivelbin = 1,num_velbin
       write(2024,*) velbin(ivelbin), binned_vel(ivelbin)
    enddo
    close(2024)
 endif 

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
