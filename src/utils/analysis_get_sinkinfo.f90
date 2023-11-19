!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which gets sink info from dumpfiles
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
 character(len=20), parameter, public :: analysistype = 'get_info'

 public :: do_analysis

 private


contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units, only:unit_velocity
 use io,    only:fatal
 use part,  only:xyzmh_ptmass,vxyz_ptmass,ihacc,nptmass
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i
 character(len=70) :: filename

 filename = 'sinkinfo_'//TRIM(dumpfile)//'.dat'
 open(unit=2024,file=filename)
 write(2024,'(2a10)') 'time','nsink'
 write(2024,*) time, nptmass
 write(2024,'(10a10)') 'id','x','y','z','m','h','vx','vy','vz','ihacc'
 do i = 1,nptmass
    write(2024,*) i,xyzmh_ptmass(1:5,i),vxyz_ptmass(1:3,i),ihacc
 enddo
 close(2024)

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
