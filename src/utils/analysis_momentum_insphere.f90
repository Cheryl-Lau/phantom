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

 integer, parameter :: nrad = 3
 integer :: isink_src = 14
 real    :: rad_list(nrad) = (/ 4.,7.,10. /) ! radii of spherical surfaces in code units 
 real    :: xyz_src_in(3) = (/ 0., 0., 0. /)            ! Position of feedback source in code units 
 logical :: use_sink = .true. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,utime,umass,unit_velocity 
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use eos,      only:gamma
 use physcon,  only:pi,solarm,km,mass_proton_cgs
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: irad,ip,nback
 real    :: time_cgs,pmass,rad2,rad_cgs,momen_tot,momen_tot_cgs,momen_part,xyz_src(3)
 real    :: dist2,r_unitvec(3),radvel,momen_kimmcen14,momen_cioffi88,nH,rho_cgs
 character(len=70) :: filename

 if (use_sink) then 
    if (isink_src > nptmass) call fatal('analysis_momentum_insphere','requested sink not found')
    xyz_src = xyzmh_ptmass(1:3,isink_src)
 else
    xyz_src = xyz_src_in 
 endif 

 filename = 'momentum_insphere_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename,status='replace')

 time_cgs  = time*utime

 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs

 !- Predicted terminal momentum from literature 
 rho_cgs = 4.d-25
 nH = rho_cgs/mass_proton_cgs
 write(2206,'(a60)') 'terminal momentum [g cm s^-1]'
 momen_cioffi88 = 4.8d5*km * (1.d0)**(13.d0/14.d0) * (1.d0)**(-3.d0/14.d0) * (nH)**(-1.d0/7.d0) * solarm 
 write(2206,'(1e60.10)') momen_cioffi88 
 momen_kimmcen14 = 3.d5*km * (1.d0)**(16.d0/17.d0) * (nH)**(-2.d0/17.d0) * solarm * (1.d0)**(-0.14)
 write(2206,'(1e60.10)') momen_kimmcen14

 !- Header 
 write(2206,'(2a30)') 'radius [cm]','radial momentum [g cm s^-1]'
 close(2206)

 !- Particle mass
 pmass = massoftype(igas)

 each_radius: do irad = 1,nrad 
    rad2 = (rad_list(irad))**2
    momen_tot = 0.
    nback = 0
    do ip = 1,npart
       dist2 = mag2(xyzh(1:3,ip) - xyz_src(1:3)) 
       if (dist2 < rad2) then 
          r_unitvec = (xyzh(1:3,ip) - xyz_src(1:3)) / sqrt(dist2)
          radvel = dot_product(vxyzu(1:3,ip),r_unitvec)
          if (radvel < 0.) nback = nback + 1 
          momen_part = pmass*radvel
          momen_tot = momen_tot + momen_part 
       endif 
    enddo 
    if (nback > 0) print*,nint(real(nback)/real(npart)*100.d0),'% of the particles went backwards!'
    if (momen_tot < 0.d0) call warning('analysis_momentum_insphere','shock reversed')
    ! convert to cgs units 
    rad_cgs = rad_list(irad) *udist 
    momen_tot_cgs = momen_tot *umass*unit_velocity 

    open(unit=2206,file=filename,position='append')
    write(2206,'(2e30.10)') rad_cgs, momen_tot_cgs
    close(2206)
 enddo each_radius 


end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
