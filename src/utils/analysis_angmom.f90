!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for calculating particle angular momentum
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies: dim, io
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'angular_momentum'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,    only:massoftype,igas
 use physcon, only:pc,solarm,years
 use units,   only:udist,umass,utime,unit_velocity
 character(len=*), intent(in)    :: dumpfile
 integer,          intent(in)    :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in)    :: particlemass,time
 integer :: icen,ifile,ncen,ip
 real    :: max_xcen,xcen,cen(3),v_vec(3),r_vec(3),angm(3)
 real    :: pmass_cgs,v_vec_cgs(3),r_vec_cgs(3),angm_cgs(3)
 character(len=50) :: filename
 character(len=5)  :: idump_char,ifile_char

 ifile = 3000
 ncen  = 10                ! number of centres to try
 max_xcen = 0.1*pc/udist   ! max distance on x-axis

 pmass_cgs = massoftype(igas)*umass

 each_centre: do icen = 1,ncen+1

    !- set centre
    xcen  = max_xcen/ncen*(icen-1)
    cen   = (/ xcen, 0., 0. /)

    write(ifile_char,'(i2.2)') icen-1
    write(idump_char,'(i5.5)') num
    filename = 'angm_'//trim(idump_char)//'_r'//trim(ifile_char)//'.txt'

    open(unit=ifile,file=trim(filename),status='replace')
    write(ifile,'(a22,e22.10)') 'particle mass [Msun]: ',pmass_cgs/solarm
    write(ifile,'(a22,e22.10)') 'time [Myr]: ',time*utime/(1E6*years)
    write(ifile,'(a22,a2,3e20.10,a2)') 'centre [pc]: ','( ',cen(1)*udist/pc, cen(2)*udist/pc, cen(3)*udist/pc,' )'
    write(ifile,*) ''

    write(ifile,'(a7,e20.10,a3)') 'udist: ',udist,'cm'
    write(ifile,'(a7,e20.10,a3)') 'utime: ',utime,'s'
    write(ifile,'(a7,e20.10,a3)') 'umass: ',umass,'g'
    write(ifile,*) ''

    write(ifile,'(a10,9a20)') 'index','r_x','r_y','r_z','v_x','v_y','v_z','L_x','L_y','L_z'
    do ip = 1,npart
       r_vec = xyzh(1:3,ip) - cen
       v_vec = vxyzu(1:3,ip)

       !- Calculate angular momentum
       r_vec_cgs = r_vec*udist
       v_vec_cgs = v_vec*unit_velocity
       angm_cgs  = pmass_cgs*cross_product(r_vec_cgs,v_vec_cgs)
       angm      = angm_cgs/(umass*udist**2/utime)
       write(ifile,'(i10,9e20.10)') ip, r_vec(1:3), v_vec(1:3), angm(1:3)
    enddo
    close(ifile)

    ifile = ifile + 1

 enddo each_centre

end subroutine do_analysis


function cross_product(vec1,vec2)
 real, intent(in)   :: vec1(3),vec2(3)
 real, dimension(3) :: cross_product

 cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
 cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
 cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

end function cross_product

end module analysis
