!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for measuring the total amount of ionized mass 
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
 character(len=20), parameter, public :: analysistype = 'ionized_mass'

 public :: do_analysis

 private

 real :: threshold_temp = 0.8E4  ! minimum temp to be considered as ionized/heated 


contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,utime,umass,unit_ergg 
 use physcon,  only:kboltz,mass_proton_cgs,solarm
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma,gmw
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: ip
 real    :: time_cgs,pmass,ionized_mass,u,temp,ionized_mass_solarm
 character(len=70) :: filename

 time_cgs = time*utime
 pmass = massoftype(igas)

 ionized_mass = 0. 
 do ip = 1,npart
    u = vxyzu(4,ip)
    temp = u/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
    if (temp > threshold_temp) then 
       ionized_mass = ionized_mass + pmass 
    endif 
 enddo 
 ionized_mass_solarm = ionized_mass*umass/solarm 


 filename = 'ionized_mass_'//TRIM(dumpfile)//'.dat'
 open(unit=2209,file=filename,status='replace')
 write(2209,'(1a20)') 'time [s]'
 write(2209,'(1e20.10)') time_cgs
 write(2209,'(1a20)') 'ionized mass [M_sun]'
 write(2209,'(1e20.10)') ionized_mass_solarm
 close(2209) 


end subroutine do_analysis


!--------------------------------------------------------------------------
end module analysis
