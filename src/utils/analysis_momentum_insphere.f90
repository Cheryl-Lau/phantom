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
 real    :: rad_list(nrad) = (/ 4., 7., 10. /)  ! radii of spherical surfaces in code units 
 real    :: xyz_src_in(3)  = (/ 0., 0., 0. /)   ! Position of feedback source in code units 
 logical :: use_sink = .false. 

 real    :: cone_omegafrac = 0.1                ! Size of solid angle 
 real    :: cone_vec(3) = (/ -4., 8., 3. /)     ! Vector of the solid angle cone wrt source 
 logical :: cone_only = .true.                  ! measure energy within a solid angle only 

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
 real    :: cone_unitvec(3),cone_mag,cone_omega,rc(3),rp(3),absrc,absrp,theta_pc,omega,xp,yp,zp 
 character(len=70) :: filename

 if (use_sink) then 
    if (isink_src > nptmass) call fatal('analysis_momentum_insphere','requested sink not found')
    xyz_src = xyzmh_ptmass(1:3,isink_src)
 else
    xyz_src = xyz_src_in 
 endif 

 if (cone_only) then 
    cone_omega = cone_omegafrac*4.*pi 
    cone_unitvec = cone_vec/sqrt(mag2(cone_vec))
    filename = 'momentum_incone_'//TRIM(dumpfile)//'.dat'
 else 
    filename = 'momentum_insphere_'//TRIM(dumpfile)//'.dat'
 endif 

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
       rp = xyzh(1:3,ip) - xyz_src(1:3)
       dist2 = mag2(rp) 
       if (dist2 < rad2) then 
          r_unitvec = rp / sqrt(dist2)

          if (cone_only) then 
             rc = cone_unitvec * sqrt(rad2)
             absrc = sqrt(mag2(rc))
             absrp = sqrt(mag2(rp))
             theta_pc = acos(dot_product(rp,rc)/(absrp*absrc))  ! angle between rp and rc 
             if (theta_pc < pi/2.d0) then                       ! within 90 deg 
                omega = 2.d0*pi*(1.d0-cos(theta_pc))            ! solid angle subtended by cone of 2.d0*theta_pc 
                if (omega < cone_omega) then
                   radvel = dot_product(vxyzu(1:3,ip),r_unitvec)
                   if (radvel < 0.) nback = nback + 1 
                   momen_part = pmass*radvel
                   momen_tot = momen_tot + momen_part 
                endif 
             endif 
          else
             radvel = dot_product(vxyzu(1:3,ip),r_unitvec)
             if (radvel < 0.) nback = nback + 1 
             momen_part = pmass*radvel
             momen_tot = momen_tot + momen_part 
          endif 
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
