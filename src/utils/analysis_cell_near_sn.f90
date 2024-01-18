!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the gas flow properties within a cell of unit volume
! placed at some distance from the supernova
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

 real    :: xcell = 10.  ! Position of cell
 real    :: ycell = 0.
 real    :: zcell = 0.
 real    :: len_cell = 5.

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,  only:udist,utime,unit_velocity,unit_density,unit_pressure
 use io,     only:fatal
 use part,   only:hfact,massoftype,igas
 use eos,    only:gamma
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: ip,npart_incell
 real    :: xmin,xmax,ymin,ymax,zmin,zmax
 real    :: pmass,x,y,z,h,rho,vx,vy,vz,v,u,thermpr,rampr
 real    :: rho_cell,v_cell,thermpr_cell,rampr_cell
 real    :: time_si,xcell_si,ycell_si,zcell_si
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Define boundaries of cell
 xmin = xcell - len_cell/2.
 xmax = xcell + len_cell/2.
 ymin = ycell - len_cell/2.
 ymax = ycell + len_cell/2.
 zmin = zcell - len_cell/2.
 zmax = zcell + len_cell/2.

 !- Extract particle properties within cell and get mean
 rho_cell     = 0.
 v_cell       = 0.
 thermpr_cell = 0.
 rampr_cell   = 0.
 npart_incell = 0
 do ip = 1,npart
    x = xyzh(1,ip)
    y = xyzh(2,ip)
    z = xyzh(3,ip)
    if (x > xmin .and. x < xmax .and. y > ymin .and. y < ymax .and. &
        z > zmin .and. z < zmax) then
       h   = xyzh(4,ip)
       rho = pmass*(hfact/abs(h))**3
       vx  = vxyzu(1,ip)
       vy  = vxyzu(2,ip)
       vz  = vxyzu(3,ip)
       v   = sqrt(vx*vx + vy*vy + vz*vz)
       u   = vxyzu(4,ip)
       thermpr = rho*(gamma-1)*u
       rampr   = rho*v**2

       rho_cell     = rho_cell + rho
       v_cell       = v_cell + v
       thermpr_cell = thermpr_cell + thermpr
       rampr_cell   = rampr_cell + rampr
       npart_incell = npart_incell + 1
    endif
 enddo
 if (npart_incell == 0) call fatal('analysis_cell_near_sn','no particles within cell')
 !if (npart_incell <= 3) call fatal('analysis_cell_near_sn','too few particles')
 rho_cell     = rho_cell/npart *unit_density *1000.
 v_cell       = v_cell/npart *unit_velocity *0.01
 thermpr_cell = thermpr_cell/npart *unit_pressure *0.1
 rampr_cell   = rampr_cell/npart *unit_pressure *0.1
 xcell_si = xcell *udist *0.01
 ycell_si = ycell *udist *0.01
 zcell_si = zcell *udist *0.01
 time_si  = time  *utime             !- all converted to SI units

 filename = 'gasflow_'//TRIM(dumpfile)//'.dat'
 open(unit=2206,file=filename)
 write(2206,'(4a20)') 'time [s]','xcell [m]','ycell [m]','zcell [m]'
 write(2206,*) time_si, xcell_si, ycell_si, zcell_si
 write(2206,'(4a25)') 'rho [kg m^-3]','vel [m s^-1]','therm pr [kg m^-1 s^-2]','ram pr [kg m^-1 s^-2]'
 write(2206,*) rho_cell,v_cell,thermpr_cell,rampr_cell

 close(2206)

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
