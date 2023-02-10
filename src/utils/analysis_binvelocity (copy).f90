!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which just prints all position and velocities
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
 character(len=20), parameter, public :: analysistype = 'print_posvel'

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units, only:unit_velocity
 use io,    only:fatal
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character(len=70) :: filename

 filename = 'posvel_'//TRIM(dumpfile)//'.dat'
 open(unit=2025,file=filename)
 do ip = 1,npart
    write(2025,*) xyzh(1:3,ip), vxyzu(1:3,ip)
 enddo
 close(2025)

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
