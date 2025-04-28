!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the gas flow radial properties (v,rho,P) through
! a shell at certain radius around the supernova
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
 character(len=20), parameter, public :: analysistype = 'sn_output'

 public :: do_analysis

 private

 integer :: isink_sn  = 14
 real    :: xyz_sn_in(3) = (/ 0.d0, 0.d0, 0.d0 /) 
 logical :: use_sink  = .true. 
 
 real    :: radius_detect = 6.0
 real    :: width_detect  = 1.0 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,umass,utime,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihsoft,ihacc
 use eos,      only:gamma,gmw
 use physcon,  only:mass_proton_cgs,kboltz
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: ip
 real    :: pmass,xyz_sn(3),r,r_unitvec(3),radvel,radmomen,u,T,rho,pressure,rampress
 real    :: radvel_cgs,radmomen_cgs,rho_cgs,pressure_cgs,rampress_cgs
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Set centre 
 if (use_sink) then 
    if (isink_sn > nptmass) call fatal('analysis_detectorshell_around_sn','sink no found.')
    xyz_sn = xyzmh_ptmass(1:3,isink_sn)
 else
    xyz_sn = xyz_sn_in 
 endif 

 filename = 'gasflow_shell_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename)

 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time*utime 
 write(2206,'(1a10,6a25)') 'ip','T [K]','rho [g cm^-3]','v_r [cm s^-1]','mu [g cm s^-1]',&
                         &'therm pr [g cm^-1 s^-2]','ram pr [g cm^-1 s^-2]'

 do ip = 1,npart 
    r = sqrt(mag2(xyzh(1:3,ip)-xyz_sn(1:3))) 
    if (r > radius_detect-width_detect .and. r < radius_detect+width_detect) then 
       r_unitvec = xyzh(1:3,ip) - xyz_sn(1:3)
       r_unitvec = r_unitvec/r 

       radvel = dot_product(vxyzu(1:3,ip),r_unitvec)
       radmomen = pmass*radvel
       u = vxyzu(4,ip) 
       T = u/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
       rho = rhoh(xyzh(4,ip),pmass)
       pressure = rho*(gamma-1.d0)*u 
       rampress = rho*radvel**2 

       ! Convert units 
       rho_cgs = rho*unit_density 
       radvel_cgs = radvel*unit_velocity 
       radmomen_cgs = radmomen*umass*unit_velocity 
       pressure_cgs = pressure*unit_pressure 
       rampress_cgs = rampress*unit_pressure 

       write(2206,'(1i10,6e25.10)') ip, T, rho_cgs, radvel_cgs, radmomen_cgs, pressure_cgs, rampress_cgs
    endif 
 enddo 

 close(2206)

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
