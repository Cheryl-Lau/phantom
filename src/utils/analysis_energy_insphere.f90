!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for measuring the energy contained within a sphere around a source 
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
 character(len=20), parameter, public :: analysistype = 'energy_in_spheres'

 public :: do_analysis

 private

 integer, parameter :: nrad = 6
 real     :: rad_list(nrad) = (/ 20.,50.,80.,100.,150.,200. /) ! radii of spherical surfaces in code units 
 real     :: xyz_src(3) = (/ 0., 0., 0. /)            ! Position of feedback source in code units 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use kernel,   only:get_kernel,cnormk,radkern2
 use units,    only:udist,utime,unit_energ 
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma
 use physcon,  only:pi 
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: irad,ip
 real    :: time_cgs,pmass,rad2,rad_cgs,ek_tot,ek_tot_cgs,ek_part,et_tot,et_tot_cgs,et_part
 real    :: dist2,totenerg_cgs 
 character(len=70) :: filename

 filename = 'energy_insphere_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename,status='replace')

 time_cgs  = time*utime

 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs
 write(2206,'(4a25)') 'radius [cm]','kinetic energy [erg]','thermal energy [erg]','total energy [erg]'
 close(2206)

 !- Particle mass
 pmass = massoftype(igas)

 each_radius: do irad = 1,nrad 
    rad2 = (rad_list(irad))**2
    ek_tot = 0.
    et_tot = 0.
    do ip = 1,npart
       dist2 = mag2(xyzh(1:3,ip) - xyz_src(1:3)) 
       if (dist2 < rad2) then 
           ek_part = 0.5*pmass*mag2(vxyzu(1:3,ip))
           et_part = vxyzu(4,ip)*pmass 
           ek_tot = ek_tot + ek_part 
           et_tot = et_tot + et_part 
       endif 
    enddo 
    ! convert to cgs units 
    rad_cgs = rad_list(irad) *udist 
    ek_tot_cgs = ek_tot *unit_energ 
    et_tot_cgs = et_tot *unit_energ 
    totenerg_cgs = ek_tot_cgs + et_tot_cgs 

    open(unit=2206,file=filename,position='append')
    write(2206,'(4e25.10)') rad_cgs, ek_tot_cgs, et_tot_cgs, totenerg_cgs 
    close(2206)
 enddo each_radius 

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis