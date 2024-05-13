!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for measuring the density of the HII region(s)
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
 character(len=20), parameter, public :: analysistype = 'hiiregion_density'

 public :: do_analysis

 private


contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,utime,unit_energ,unit_ergg,unit_density 
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma,gmw 
 use physcon,  only:kboltz,mass_proton_cgs
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: ip
 real    :: pmass,time_cgs,temp,etherm_cgs,rho_cgs,x_cgs,y_cgs,z_cgs,rho_mean
 character(len=70) :: filename

 filename = 'hii_dens_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename,status='replace')

 time_cgs  = time*utime

 write(2206,'(1a20)') 'time [s]'
 write(2206,'(1e20.10)') time_cgs
 write(2206,'(6a25)') 'x [cm]','y [cm]','z [cm]','temp [K]','etherm [erg]','density [g cm^-3]'
 close(2206)

 !- Particle mass
 pmass = massoftype(igas)

 rho_mean = 0. 
 do ip = 1,npart
    temp = vxyzu(4,ip)/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
    if (temp > 5E3) then 
       etherm_cgs = vxyzu(4,ip)*pmass *unit_energ 
       rho_cgs = rhoh(xyzh(4,ip),pmass) *unit_density 
       x_cgs = xyzh(1,ip) *udist 
       y_cgs = xyzh(2,ip) *udist
       z_cgs = xyzh(3,ip) *udist
       write(2206,'(6a25)') x_cgs, y_cgs, z_cgs, temp, etherm_cgs, rho_cgs 
       rho_mean = rho_mean + rho_cgs
    endif 
 enddo 
 rho_mean = rho_mean/npart 

 print*,'Average density of HII region:', rho_mean 
 
end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
