!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for measuring the radial momentum contained within a sphere around a source 
! for a range of sphere radii 
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
 character(len=20), parameter, public :: analysistype = 'momentum_in_spheres'

 public :: do_analysis

 private

 integer, parameter :: nrad = 4
! real     :: rad_list(nrad) = (/ 20.,50.,80.,100.,150.,200. /) ! radii of spherical surfaces in code units 
 real     :: rad_list(nrad) = (/ 8.,15.,25.,35. /) 
 real     :: xyz_src(3) = (/ 0., 0., 0. /)            ! Position of feedback source in code units 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use kernel,   only:get_kernel,cnormk,radkern2
 use units,    only:udist,utime,umass,unit_velocity 
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma
 use physcon,  only:pi 
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: irad,ip
 real    :: time_cgs,pmass,rad2,rad_cgs,momen_tot,momen_tot_cgs,momen_part
 real    :: dist2,r_unitvec(3),radvel
 character(len=70) :: filename

 filename = 'momentum_insphere_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename,status='replace')

 time_cgs  = time*utime

 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs
 write(2206,'(2a25)') 'radius [cm]','radial momentum [g cm s^-1]'
 close(2206)

 !- Particle mass
 pmass = massoftype(igas)

 each_radius: do irad = 1,nrad 
    rad2 = (rad_list(irad))**2
    momen_tot = 0.
    do ip = 1,npart
       dist2 = mag2(xyzh(1:3,ip) - xyz_src(1:3)) 
       if (dist2 < rad2) then 
           r_unitvec = (xyzh(1:3,ip) - xyz_src(1:3)) / sqrt(dist2)
           radvel = dot_product(vxyzu(1:3,ip),r_unitvec)
           momen_part = pmass*radvel
           momen_tot = momen_tot + momen_part 
       endif 
    enddo 
    ! convert to cgs units 
    rad_cgs = rad_list(irad) *udist 
    momen_tot_cgs = momen_tot *umass*unit_velocity 


    open(unit=2206,file=filename,position='append')
    write(2206,'(2e25.10)') rad_cgs, momen_tot_cgs
    close(2206)
 enddo each_radius 

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
