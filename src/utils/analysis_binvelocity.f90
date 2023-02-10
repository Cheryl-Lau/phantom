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
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies:
!

 implicit none
 character(len=20), parameter, public :: analysistype = 'bin_velocity'

 public :: do_analysis

 private
 integer, parameter :: num_velbin = 100
 integer :: binned_vel(num_velbin)
 real    :: logvelbin_min = -1.
 real    :: logvelbin_max = 5.

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units, only:unit_velocity
 use io,    only:fatal
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,ivelbin
 real    :: vel,dvel,velbin(num_velbin)
 character(len=70) :: filename

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
 do i = 1,npart
    !
    ! Calculate |v|
    !
    vel = sqrt(vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2)
    !
    ! Bin velocity
    !
    ivelbin = (log10(vel)-logvelbin_min+dvel)/dvel
    if (ivelbin < 0) call fatal('analysis_binvelocity','require smaller logvelbin_min')
    if (ivelbin > num_velbin) call fatal('analysis_binvelocity','require larger logvelbin_max')
    binned_vel(ivelbin) = binned_vel(ivelbin) + 1
 enddo
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

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
