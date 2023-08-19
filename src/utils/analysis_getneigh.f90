!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for getting neighbours
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, kdtree, linklist
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'neighbours'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use kdtree,   only:getneigh
 use linklist, only:node,ifirstincell,listneigh
 use dim,      only:maxneigh
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in)    :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: isrc,nneigh,ixyzcachesize,ineigh,ipart
 real    :: xyzcache(3,neighcachesize)
 real    :: target_pos(3),rad_thresh,x

 !- set target location and threshold radius, e.g.
 target_pos = (/ 0.,0.,0. /)
 rad_thresh = 0.3

 call getneigh(node,target_pos,0.,rad_thresh,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)

 !- retrieve neighbours
 do ineigh = 1,nneigh
    ipart = listneigh(ineigh)   ! particle index

    ! Get particle position
    x = xyzh(1,ipart)
    ! or get it from cache list for fast retrieval
    x = xyzcache(1,ineigh)
 enddo

end subroutine do_analysis

end module analysis
