!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the total mass, kinetic energy, thermal energy and 
! momentum of particles in feedback outflows beyond the cloud
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
 character(len=20), parameter, public :: analysistype = 'cloud_extern'

 public :: do_analysis

 private

 real    :: rcloud = 8. 
 real    :: rmax   = 50. 
 real    :: xyz_src_in(3) = (/ 0.d0, 0.d0, 0.d0 /) 
 logical :: use_sink  = .true. 
 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,umass,utime,unit_velocity,unit_energ
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:xyzmh_ptmass,vxyz_ptmass,nptmass,ihsoft,ihacc
 use physcon,  only:solarm
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,io,isink_src
 real    :: rcloud2,rmax2,r2,rvec(3),r_unitvec(3)
 real    :: pmass,time_cgs,xyz_src(3),radvel,vel2,vel,momen
 real    :: totmass,totekin,totetherm,tote,totmomen
 character(len=70) :: filename

 rcloud2 = rcloud**2
 rmax2 = rmax**2 

 !- Set centre 
 if (use_sink) then 
    open(unit=2204,file='sink_src.txt',status='old',iostat=io)
    if (io /= 0) call fatal('analysis_cloud_extern','cannot open sink_src.txt')
    read(2204,*) isink_src
    close(2204) 
    print*, 'Centre: Sink ', isink_src
    if (isink_src > nptmass) call fatal('analysis_cloud_extern','sink no found.')
    xyz_src = xyzmh_ptmass(1:3,isink_src)
 else
    xyz_src = xyz_src_in 
 endif 


 filename = 'extprop_'//TRIM(dumpfile)//'.dat'
 open(unit=2207,file=filename,status='replace',iostat=io)
 if (io /= 0) call fatal('analysis_cloud_extern','cannot open file to write')
 close(2207)

 time_cgs  = time*utime

 !- Particle mass
 pmass = massoftype(igas)

 totmass = 0
 totekin = 0
 totetherm = 0
 totmomen = 0

 do i = 1,npart
    rvec = xyzh(1:3,i)-xyz_src(1:3)
    r2 = mag2(rvec)

    if (r2 > rcloud2 .and. r2 < rmax2) then 
       totmass = totmass + pmass
       totekin = totekin + 5d-1*pmass*mag2(vxyzu(1:3,i))
       totetherm = totetherm + vxyzu(4,i)*pmass 

       r_unitvec = rvec / sqrt(r2)
       radvel = dot_product(vxyzu(1:3,i),r_unitvec)
       vel2 = mag2(vxyzu(1:3,i))
       vel  = sqrt(vel2)
       momen = pmass*vel
       totmomen = totmomen + momen 
    endif 
 enddo 
 
 !- Convert to useful units 
 totmass = totmass*umass/solarm
 totekin = totekin*unit_energ
 totetherm = totetherm*unit_energ
 tote = totekin + totetherm
 totmomen = totmomen*umass*unit_velocity

 !- Write results 
 open(unit=2207,file=filename,status='old')
 write(2207,'(1a20)') 'time [s]'
 write(2207,'(1e20.10)') time_cgs
 write(2207,'(5a20)') 'mass [msun]','ekin [erg]','etherm [erg]','etot [erg]','momen [g cm s^-1]'
 write(2207,'(5e20.10)') totmass, totekin, totetherm, tote, totmomen
 close(2207)


end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
