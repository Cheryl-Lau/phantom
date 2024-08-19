!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the temperature within given ionization front radius
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
 character(len=20), parameter, public :: analysistype = 'hii_temp'

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,utime,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use physcon,  only:kboltz,mass_proton_cgs
 use eos,      only:gamma,gmw 
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: nentry,ientry,io,ip,npart_inrad
 real    :: pmass,time_cgs,mindiff,time_fromfile,ionrad,ionrad2,meantemp,r2,temp 
 real,   allocatable :: time_Myr(:),ionrad_pc(:),spitz_pc(:),hosoinu_pc(:)
 character(len=70) :: in_filename,out_filename

 !- Particle mass
 pmass = massoftype(igas)
 !- Time 
 time_cgs = time *utime 

 !- Import the ionization fronts 
 in_filename = 'rho19_implicit_radevol.dat'
 open(2205,file=in_filename,status='old')
 nentry = 0
 do
 	 read(2205,*,iostat=io)
	 if (io < 0) then
	 	 exit
    else
	 	 nentry = nentry + 1
	 endif
 enddo
 print*,'number of entries in file',nentry

 allocate(time_Myr(nentry))
 allocate(ionrad_pc(nentry))
 allocate(spitz_pc(nentry))
 allocate(hosoinu_pc(nentry))

 rewind(2205)
 do ientry = 1,nentry
	 read(2205,*) time_Myr(ientry), ionrad_pc(ientry), spitz_pc(ientry), hosoinu_pc(ientry) 
 enddo
 close(2205) 

 !- Find current ionization front 
 mindiff = huge(mindiff)
 do ientry = 1,nentry
    time_fromfile = (time_Myr(ientry)*1e6*365*24*60*60) /utime 
    if (abs(time - time_fromfile) < mindiff) then 
       mindiff = abs(time - time_fromfile)
       ionrad = ionrad_pc(ientry)
    endif 
 enddo 
 deallocate(time_Myr)
 deallocate(ionrad_pc)
 deallocate(spitz_pc)
 deallocate(hosoinu_pc)

 print*,'current ionization front radius',ionrad 

 ionrad2 = (ionrad*0.8)**2  ! slightly inwards 
 
 !- Get mean temperature inside radius 
 meantemp = 0.
 npart_inrad = 0
 do ip = 1,npart 
    r2 = mag2(xyzh(1:3,ip))
    if (r2 < ionrad2) then 
       temp = vxyzu(4,ip)/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
       meantemp = meantemp + temp 
       npart_inrad = npart_inrad + 1 
    endif 
 enddo
 meantemp = meantemp/npart_inrad 

 out_filename = 'tempHII_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=out_filename)
 write(2206,'(2a20)') 'time [s]','temp [K]'
 write(2206,'(2e20.10)') time_cgs, meantemp
 close(2206)

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
