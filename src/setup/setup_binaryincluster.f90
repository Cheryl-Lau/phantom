!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up a rotating sphere with randomly-positioned sinks and
! a massive sink at the centre (single star or binary)
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - totmass_req      : *mass of sphere in solarm*
!   - pmass            : *particle mass resolution in solarm*
!   - r_sphere         : *radius of sphere in pc*
!   - mach             : *turbulence mach number*
!   - angvel_cgs       : *angular velocity in rad/s*
!   - cs_cgs           : *sound speed in sphere in cm/s*
!   - nptmass_clust    : *number of sinks in the cluster*
!   - binary_cen_sink  : *option to set the central sink as a binary*
!   - pin_cen_sink     : *option to stop the central sink from moving*
!   - make_sinks       : *option to create sinks dynamically*
!
! :Dependencies: 
!
 implicit none
 public  :: setpart

 private
 integer :: nptmass_clust
 real    :: totmass_req,pmass,r_sphere
 real    :: mach,angvel_cgs,cs_cgs
 logical :: binary_cen_sink,pin_cen_sink,make_sinks
 logical :: isotherm    = .true.    ! isothermal; otherwise adiabatic 

 integer :: iseed = -12345

contains

!----------------------------------------------------------------
!+
!  setup for a sphere with sinks 
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,pc,au,kboltz,mass_proton_cgs,gg 
 use dim,          only:maxvxyzu
 use setup_params, only:rhozero,npart_total
 use io,           only:master,fatal,iverbose
 use spherical,    only:set_sphere
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,udist,unit_density,unit_velocity
 use eos,          only:ieos,gmw
 use part,         only:set_particle_type,igas,nptmass,xyzmh_ptmass,vxyz_ptmass
 use timestep,     only:dtmax,tmax,nout,dtwallmax
 use centreofmass, only:reset_centreofmass
 use options,      only:nfulldump,nmaxdumps,icooling,iexternalforce,ishock_heating,ipdv_heating
 use kernel,       only:hfact_default
 use domain,       only:i_belong
 use ptmass,       only:icreate_sinks,rho_crit,rho_crit_cgs,r_crit,h_acc,h_soft_sinksink,h_soft_sinkgas
 use ptmass,       only:pin_sink,pin_all,isink_to_pin 
 use options,      only:iexternalforce
 use cooling,      only:Tfloor
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use random,       only:ran2
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=20), parameter     :: filevx = 'cube_v1.dat'
 character(len=20), parameter     :: filevy = 'cube_v2.dat'
 character(len=20), parameter     :: filevz = 'cube_v3.dat'
 integer :: ierr,npart_req,i,isink 
 real    :: cs,angvel,totmass,t_ff,psep,rmsmach,v2i,turbfac,turbboxsize,temp,u
 real    :: sep,mbinary,angmom,h_acc0
 real    :: lowmass_innersink,lowmass_outersink,uppmass_innersink,uppmass_outersink
 real    :: r_thresh2,r2,mass,x,y,z 
 logical :: setexists,inexists
 character(len=120) :: filex,filey,filez
 character(len=100) :: filein,fileset,cwd

 if (isotherm .and. maxvxyzu >= 4) call fatal('setup_binaryincluster','set ISOTHERMAL=yes in Makefile')
 if (maxvxyzu == 3 .and. .not.isotherm) call fatal('setup_binaryincluster','set isotherm=.true. in setup')

 !--Set units
 call set_units(dist=pc,mass=solarm,G=1.)

 !--Default values for the input params 
 totmass_req      = 1d2
 pmass            = 1d-3
 r_sphere         = 0.1
 mach             = 20.0 
 angvel_cgs       = 1.006d-12 
 cs_cgs           = 2.19d4  ! 8K assuming mu = 2.31 & gamma = 5/3
 nptmass_clust    = 50
 binary_cen_sink  = .true. 
 pin_cen_sink     = .false.
 make_sinks       = .true. 

 !--Check for existence of the .in and .setup files
 filein = trim(fileprefix)//'.in'
 inquire(file=filein,exist=inexists)
 fileset = trim(fileprefix)//'.setup'
 inquire(file=fileset,exist=setexists)

 !--Read values from .setup
 if (setexists) then
    call read_setupfile(fileset,ierr)
    if (ierr /= 0) then
       if (id == master) call write_setupfile(fileset)
       stop
    endif
    !--Prompt to get inputs and write to file
 elseif (id == master) then
    print "(a,/)",trim(fileset)//' not found: using interactive setup'
    call get_input_from_prompts()
    call write_setupfile(fileset)
 endif

 !--Convert units 
 cs      = cs_cgs/unit_velocity 
 angvel  = angvel_cgs*utime

 !--Set sphere 
 npart_req = nint(totmass_req/pmass)
 if (npart_req > size(xyzh(1,:))) call fatal('setup_binaryincluster','npart_req exceeded limit')
 call set_sphere('closepacked',id,master,0.,r_sphere,psep,hfact_default,npart,xyzh,nptot=npart_total, &
                 exactN=.true.,np_requested=npart_req,mask=i_belong)

 !--Set particle properties
 npartoftype(:) = 0
 npartoftype(1) = npart
 massoftype(1)  = pmass 
 do i = 1,npart
    call set_particle_type(i,igas)
 enddo

 !--Recomptue cloud properties 
 totmass    = npart*pmass
 rhozero    = totmass/(4./3.*pi*r_sphere**3)
 t_ff       = sqrt(3.*pi/(32.*rhozero)) 
 if (abs(totmass-totmass_req)/totmass_req > 0.1) call fatal('setup_binaryincluster','mass does not match requested value')

 !--Impose turbulent velocity field 
 vxyzu = 0.
 turbulent: if (mach > 0.) then
    call getcwd(cwd)

    filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
    filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
    filez  = find_phantom_datafile(filevz,'velfield_sphng_small')

    turbboxsize = r_sphere
    call set_velfield_from_cubes(xyzh(:,1:npart),vxyzu(:,:npart),npart, &
                                 filex,filey,filez,1.,turbboxsize,.false.,ierr)
    if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

    rmsmach = 0.0
    do i = 1,npart
       v2i     = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       rmsmach = rmsmach + v2i/cs**2
    enddo
    rmsmach = sqrt(rmsmach/npart)
    if (rmsmach > 0.) then
       turbfac = mach/rmsmach ! normalise the energy to the desired mach number
     else
       turbfac = 0.
    endif
    do i = 1,npart
       vxyzu(1:3,i) = turbfac*vxyzu(1:3,i)
    enddo
 endif turbulent

 !--Impose uniform rotation velocity 
 rotating: if (angvel > 0.) then 
    do i = 1,npart
       vxyzu(1,i) = vxyzu(1,i) - angvel*xyzh(2,i)
       vxyzu(2,i) = vxyzu(2,i) + angvel*xyzh(1,i)
    enddo
 endif rotating 

 !--Setting the centre of mass of the cloud to be zero
 call reset_centreofmass(npart,xyzh,vxyzu)

 !--Thermodynamic properties  
 gamma   = 5./3.           
 Tfloor  = 3.
 temp    = 10.
 polyk   = kboltz*temp/(gmw*mass_proton_cgs)*(utime/udist)**2
 u       = polyk/(gamma-1.)
 if (maxvxyzu >= 4) vxyzu(4,1:npart) = u

 !--Set sinks 
 nptmass = nptmass_clust 

 !--Place central massive sink 
 if (nptmass > 0) then 
    isink = 1
    xyzmh_ptmass(1:3,isink) = (/ 0.,0.,0. /)
    if (binary_cen_sink) then 
       sep      = 100*au/udist                     ! binary separation 
       mbinary  = 8. * 2.                          ! total mass of the binary 
       angmom   = 2.5d-1 * sqrt(mbinary**3 * sep)  ! for equal mass pairs 
       xyzmh_ptmass(4,isink)  = mbinary ! mass 
       xyzmh_ptmass(5,isink)  = sep     ! h_acc 
       xyzmh_ptmass(6,isink)  = sep     ! h_soft 
       xyzmh_ptmass(8,isink)  = 0.      ! spinx 
       xyzmh_ptmass(9,isink)  = 0.      ! spiny 
       xyzmh_ptmass(10,isink) = angmom  ! spinz 
    else 
       h_acc0 = 5*au/udist 
       xyzmh_ptmass(4,isink)  = 8. * 2. ! mass
       xyzmh_ptmass(5,isink)  = h_acc0  ! h_acc
       xyzmh_ptmass(6,isink)  = h_acc0  ! h_soft 
       xyzmh_ptmass(8,isink)  = 0.      ! spinx 
       xyzmh_ptmass(9,isink)  = 0.      ! spiny 
       xyzmh_ptmass(10,isink) = 0.      ! spinz 
   endif 
 endif 

 !--Place the rest of the sinks 
 if (nptmass > 1) then 
    uppmass_innersink   = 8.                 ! range of stellar masses in inner regions 
    lowmass_innersink   = 4. 
    uppmass_outersink   = 5.                 ! range of stellar masses in outer regions 
    lowmass_outersink   = 2. 
    r_thresh2           = (0.5*r_sphere)**2  ! where to divide the regions 

    do isink = 2,nptmass_clust
       call gen_random_pos(r_sphere,x,y,z)
       r2 = x**2 + y**2 + z**2 
       if (r2 < r_thresh2) then 
          call gen_random_mass(lowmass_innersink,uppmass_innersink,mass)
       else 
          call gen_random_mass(lowmass_outersink,uppmass_outersink,mass)
       endif 
       h_acc0 = 5*au/udist
       xyzmh_ptmass(1,isink)  = x 
       xyzmh_ptmass(2,isink)  = y 
       xyzmh_ptmass(3,isink)  = z 
       xyzmh_ptmass(4,isink)  = mass 
       xyzmh_ptmass(5,isink)  = h_acc0
       xyzmh_ptmass(6,isink)  = h_acc0
       xyzmh_ptmass(8,isink)  = 0.      ! spinx 
       xyzmh_ptmass(9,isink)  = 0.      ! spiny 
       xyzmh_ptmass(10,isink) = 0.      ! spinz 
       vxyz_ptmass(1,isink)   = -1.d0*angvel*xyzmh_ptmass(2,isink)
       vxyz_ptmass(2,isink)   = angvel*xyzmh_ptmass(1,isink)
    enddo 
 endif 

 !--Set options for input file, if .in file does not exist
 if (.not.inexists) then
    tmax      = 1.*t_ff
    dtmax     = 1.d-3*t_ff
    nout      = 10
    nfulldump = 1
    nmaxdumps = 1000
    dtwallmax = 86400   ! s
    iverbose  = 1

    if (isotherm) then 
       ieos           = 1
       icooling       = 0
    else 
       ieos           = 2
       icooling       = 7
       ipdv_heating   = 1
       ishock_heating = 1
    endif 

    iexternalforce = 17  ! cluster potential

    if (make_sinks) then 
       icreate_sinks    = 1
       h_acc            = 1.d-3
       r_crit           = 5.0*h_acc
       rho_crit_cgs     = 1.d-14 
       rho_crit         = rho_crit_cgs/unit_density
       h_soft_sinkgas   = h_acc
       h_soft_sinksink  = h_acc
    else 
       icreate_sinks    = 0
    endif 

    if (pin_cen_sink) then 
       pin_sink         = .true. 
       isink_to_pin     = 1
    endif 
 endif 

end subroutine setpart


subroutine gen_random_pos(maxr,x,y,z)
 use random, only:ran2
 real, intent(in)  :: maxr 
 real, intent(out) :: x,y,z
 real :: r

 r = maxr + 1 
 do while (r > maxr)
    x = 2.*(ran2(iseed)-0.5) *maxr 
    y = 2.*(ran2(iseed)-0.5) *maxr 
    z = 2.*(ran2(iseed)-0.5) *maxr 
    r = sqrt(x**2 + y**2 + z**2) 
 enddo 

end subroutine gen_random_pos 


subroutine gen_random_mass(lowbound,uppbound,mass)
 use random, only:ran2
 real, intent(in)  :: lowbound,uppbound 
 real, intent(out) :: mass 

 mass = lowbound + ran2(iseed)*(uppbound-lowbound)

end subroutine gen_random_mass 

!----------------------------------------------------------------
!
!  Prompt user for inputs
!
!----------------------------------------------------------------
subroutine get_input_from_prompts()
 use prompting, only:prompt

 call prompt('Enter the requested mass of the sphere (in Msun)',totmass_req)
 call prompt('Enter the particle mass (in Msun)',pmass)
 call prompt('Enter the radius of the sphere (in pc)',r_sphere)
 call prompt('Enter the turbulent mach number',mach)
 call prompt('Enter the angular velocity in rad s^-1',angvel_cgs)
 call prompt('Enter the sound speed in cm s^-1',cs_cgs)
 call prompt('Enter the number of stars in the cluster',nptmass_clust)
 call prompt('Is the central star a binary',binary_cen_sink)
 call prompt('Do you wish to pin the central star',pin_cen_sink)
 call prompt('Do you wish to dynamically create sinks during runtime',make_sinks)

end subroutine get_input_from_prompts

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20
 integer                      :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(a)") '# input file for binary-in-cluster setup routines'
 call write_inopt(totmass_req,'totmass_req','requested total mass in sphere',iunit)
 call write_inopt(pmass,'pmass','particle mass',iunit)
 call write_inopt(r_sphere,'r_sphere','radius of sphere',iunit)
 call write_inopt(mach,'mach','turbulent mach number',iunit)
 call write_inopt(angvel_cgs,'angvel_cgs','angular velocity of sphere',iunit)
 call write_inopt(cs_cgs,'cs_cgs','sound speed in sphere',iunit)
 call write_inopt(nptmass_clust,'nptmass_clust','number of sinks in sphere',iunit)
 call write_inopt(binary_cen_sink,'binary_cen_sink','set the central star as a binary',iunit)
 call write_inopt(pin_cen_sink,'pin_cen_sink','pin the central star',iunit)
 call write_inopt(make_sinks,'make_sinks','dynamically create sink particles',iunit)

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: i,nerr,jerr,kerr
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(totmass_req,'totmass_req',db,ierr)
 call read_inopt(pmass,'pmass',db,ierr)
 call read_inopt(r_sphere,'r_sphere',db,ierr)
 call read_inopt(mach,'mach',db,ierr)
 call read_inopt(angvel_cgs,'angvel_cgs',db,ierr)
 call read_inopt(cs_cgs,'cs_cgs',db,ierr)
 call read_inopt(nptmass_clust,'nptmass_clust',db,ierr)
 call read_inopt(binary_cen_sink,'binary_cen_sink',db,ierr)
 call read_inopt(pin_cen_sink,'pin_cen_sink',db,ierr)
 call read_inopt(make_sinks,'make_sinks',db,ierr)

 call close_db(db)

 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_binaryincluster: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

end module setup
