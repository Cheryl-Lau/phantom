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
 character(len=20), parameter, public :: analysistype = 'mass_energ_momen_in_spheres'

 public :: do_analysis

 private

 integer, parameter :: nrad = 5
 real    :: rad_list(nrad) = (/ 5., 10., 15., 20., 30. /)  ! radii of spherical surfaces in code units 
 real    :: xyz_src_in(3)  = (/ 0., 0., 0. /)   ! Position of feedback source in code units 
 logical :: use_sink = .true. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,utime,umass,unit_energ,unit_velocity
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use eos,      only:gamma
 use physcon,  only:pi,solarm,pc
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: irad,ip,isink_src
 real    :: time_cgs,pmass,rad2,rad_pc,ek_tot,ek_tot_cgs,ek_part,et_tot,et_tot_cgs,et_part
 real    :: rp(3),r_unitvec(3),dist2,totenerg_cgs,xyz_src(3)
 real    :: mass_tot,momen_tot,momen_part,momen_tot_cgs,radvel 
 character(len=70) :: filename_mass,filename_energ,filename_momen


 if (use_sink) then 
    open(unit=2204,file='sink_src.txt',status='old')
    read(2204,*) isink_src
    close(2204) 
    print*, 'Centre: Sink ', isink_src
    if (isink_src > nptmass) call fatal('analysis_mass_energy_momentum_insphere','requested sink not found')
    xyz_src = xyzmh_ptmass(1:3,isink_src)
 else
    xyz_src = xyz_src_in 
 endif 

 time_cgs  = time*utime
 pmass     = massoftype(igas)*umass/solarm 



 filename_mass  = 'mass_insphere_'//TRIM(dumpfile)//'.dat'
 open(unit=2205,file=filename_mass,status='replace')
 write(2205,'(1a20)') 'time [s]'
 write(2205,'(1e20.10)') time_cgs
 write(2205,'(2a25)') 'radius [pc]','mass [Msun]'
 close(2205)

 filename_energ = 'energy_insphere_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename_energ,status='replace')
 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs
 write(2206,'(4a25)') 'radius [pc]','kinetic energy [erg]','thermal energy [erg]','total energy [erg]'
 close(2206)

 filename_momen = 'momen_insphere_'//TRIM(dumpfile)//'.dat'
 open(unit=2207,file=filename_momen,status='replace')
 write(2207,'(1a20)') 'time [s]'
 write(2207,'(1e20.10)') time_cgs 
 write(2207,'(2a25)') 'radius [pc]','momentum [g cm/s]'
 close(2207)
 



 each_radius: do irad = 1,nrad 

    rad2 = (rad_list(irad))**2
    mass_tot = 0. 
    ek_tot = 0.
    et_tot = 0.
    momen_tot = 0.

    over_particles: do ip = 1,npart
       rp = xyzh(1:3,ip) - xyz_src(1:3)
       dist2 = mag2(rp) 

       if (dist2 < rad2) then 

          !-Mass
          mass_tot = mass_tot + pmass 

          !-Energies 
          ek_part = 0.5*pmass*mag2(vxyzu(1:3,ip))
          et_part = vxyzu(4,ip)*pmass 
          ek_tot = ek_tot + ek_part 
          et_tot = et_tot + et_part    

          !-Momentum radial component 
          r_unitvec = rp / sqrt(dist2)
          radvel = dot_product(vxyzu(1:3,ip),r_unitvec)
          momen_part = pmass*radvel
          momen_tot = momen_tot + abs(momen_part)

       endif 
    enddo over_particles


    rad_pc = rad_list(irad)*udist/pc 

    ! convert energy and momentum to cgs units 
    ek_tot_cgs = ek_tot *unit_energ 
    et_tot_cgs = et_tot *unit_energ 
    totenerg_cgs = ek_tot_cgs + et_tot_cgs 
    momen_tot_cgs = momen_tot *umass*unit_velocity 

    ! Write to file 
    open(unit=2205,file=filename_mass,position='append')
    write(2205,'(2e25.10)') rad_pc, mass_tot
    close(2205)

    open(unit=2206,file=filename_energ,position='append')
    write(2206,'(4e25.10)') rad_pc, ek_tot_cgs, et_tot_cgs, totenerg_cgs 
    close(2206)

    open(unit=2207,file=filename_momen,position='append')
    write(2207,'(2e25.10)') rad_pc, momen_tot_cgs
    close(2207)


 enddo each_radius 

end subroutine do_analysis



real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
