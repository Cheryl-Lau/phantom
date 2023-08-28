!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up a sphere or ellipsoid with no surrounding medium
! Includes settings for injecting photoionization and supernovae
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - dist_unit        : *distance unit (e.g. pc)*
!   - mass_unit        : *mass unit (e.g. solarm)*
!   - is_sphere        : *set up a sphere; otherwise an ellipsoid*
!   - npart_sphere     : *requested number of particles in sphere*
!   - r_sphere         : *radius of sphere in code units*
!   - r_ellipsoid      : *semi-axes of ellipsoid a, b and c in code units*
!   - totmass_sphere   : *mass of sphere/ellipsoid in code units*
!   - rms_mach         : *turbulence mach number*
!   - gauss_density    : *flag to set a gaussian distribution for density within sphere*
!   - make_sinks       : *flag to dynamically create sinks*
!
! :Dependencies: dim, eos, infile_utils, io, options, part, physcon, prompting, ptmass,
!   rho_profile, units
!
!
 implicit none
 public :: setpart

 private

 integer      :: npart_sphere
 real         :: totmass_sphere,r_sphere,r_ellipsoid(3),temp_sphere,rms_mach
 real(kind=8) :: udist,umass
 logical      :: is_sphere,gauss_density,apply_cooling,make_sinks
 character(len=20) :: dist_unit,mass_unit

contains

!----------------------------------------------------------------
!+
!  setup for a sphere or ellipsoid
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,years,pc,kboltz,mass_proton_cgs,Rg,gg
 use dim,          only:maxvxyzu
 use setup_params, only:npart_total
 use io,           only:master,fatal,iverbose
 use spherical,    only:set_sphere,set_ellipse
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_velocity,unit_ergg
 use eos,          only:ieos,gmw
 use part,         only:igas,set_particle_type
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use ptmass,       only:icreate_sinks,rho_crit_cgs,r_crit,h_acc,h_soft_sinksink,h_soft_sinkgas
 use timestep,     only:dtmax,tmax,dtwallmax,nout
 use options,      only:nfulldump,nmaxdumps,icooling
 use kernel,       only:hfact_default
 use domain,       only:i_belong
 use cooling,      only:Tfloor,ufloor
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use options,      only:ipdv_heating,ishock_heating
 use stretchmap,   only:rho_func
 use photoionize_cmi,  only:monochrom_source,fix_temp_hii,treat_Rtype_phase
 use photoionize_cmi,  only:photoionize_tree,tree_accuracy_cmi,nHlimit_fac
 use photoionize_cmi,  only:rcut_opennode_cgs,rcut_leafpart_cgs,delta_rcut_cgs
 use photoionize_cmi,  only:sink_ionsrc,inject_rad
 use inject,           only:inject_sn,sink_progenitor
 procedure(rho_func),  pointer :: density_func
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
 integer            :: ip,npart_max,ierr
 real               :: psep,vol_sphere,dens_sphere,u_sphere,cs_sphere,cs_sphere_cgs,dens_sphere_cgs
 real               :: rmax,t_ff,rmsmach,v2i,turbfac,turbboxsize
 real               :: v2_sum,v_rms,v_rms_kms
 real               :: jeans_mass,jeans_mass_cgs
 real               :: h_acc_cgs,h_soft_sinksink_cgs,h_soft_sinkgas_cgs,rho_crit_cgs_recomm
 logical            :: iexist,in_iexist
 logical            :: place_sink_in_setup = .false.
 character(len=120) :: filex,filey,filez,filename,infilename,cwd,lattice
 character(len=10)  :: c_shape,proceed

 print "(a)", 'Setup for Sphere or Ellipsoid in which radiation and/or supernovae can be injected'

 filename = trim(fileprefix)//'.setup'
 inquire(file = filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    if (ierr /= 0) then
       if (id == master) call write_setupfile(filename)
       stop
    endif
 elseif (id == master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    dist_unit = 'pc'
    mass_unit = 'solarm'
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo
    !
    ! units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)
    !
    ! Limit in number of particles
    !
    npart_max = int(0.9*size(xyzh(1,:)))

    !
    ! prompt user for settings
    !
    c_shape = 'sphere'
    call prompt('Is the cloud initially a [sphere] or an [ellipsoid]?',c_shape)

    if (c_shape == 'sphere') then
       is_sphere = .true.

       npart_sphere = 1E5
       call prompt('Enter the approximate number of particles in the sphere',npart_sphere,0,npart_max)

       r_sphere = 0.5
       call prompt('Enter radius of sphere in units of '//dist_unit,r_sphere,0.)

       totmass_sphere = 1E3
       call prompt('Enter total mass in sphere in units of '//mass_unit,totmass_sphere,0.)

       gauss_density = .false.
       call prompt('Do you wish to centrally condense the sphere? ',gauss_density)

    elseif (c_shape == 'ellipsoid') then
       is_sphere = .false.

       npart_sphere = 1E5
       call prompt('Enter the approximate number of particles in the ellipsoid',npart_sphere,0,npart_max)

       r_ellipsoid(1) = 5.
       r_ellipsoid(2) = 2.
       r_ellipsoid(3) = 2.
       call prompt('Enter the semi-axis {a} of ellipsoid in units of '//dist_unit,r_ellipsoid(1),0.)
       call prompt('Enter the semi-axis {b} of ellipsoid in units of '//dist_unit,r_ellipsoid(2),0.)
       call prompt('Enter the semi-axis {c} of ellipsoid in units of '//dist_unit,r_ellipsoid(3),0.)

       totmass_sphere = 1E4
       call prompt('Enter total mass in ellipsoid in units of '//mass_unit,totmass_sphere,0.)

       !- density map cannot be used on ellipsoids
       gauss_density = .false.

    else
       print "(a)",' ERROR: shape not recognised - Cloud has to be either a sphere or an ellipsoid.'
       stop
    endif

    apply_cooling = .false.
    call prompt('Do you wish to apply cooling to the gas particles? ',apply_cooling)

    temp_sphere = 20.
    if (apply_cooling) then
       print*, 'Initital temperature set to ',temp_sphere,' K; will drift to thermal equilibrium.'
    else
       call prompt('Enter temperature of sphere in K',temp_sphere,0.)
    endif

    rms_mach = 5.5
    call prompt('Enter the Mach number of the cloud turbulence',rms_mach,0.)

    make_sinks = .true.
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)

    if (id == master) call write_setupfile(filename)
    stop 'please edit .setup file and rerun phantomsetup'
 else
    stop ! MPI, stop on other threads, interactive on master
 endif

 !
 ! units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 ! set gamma
 !
 if (maxvxyzu >= 4) then
    gamma = 5./3.
 else
    gamma = 1.
 endif

 !
 ! General parameters
 !
 time  = 0.
 hfact = hfact_default

 if (is_sphere) then
    rmax = r_sphere
    vol_sphere = 4./3.*pi*r_sphere**3
 else
    rmax = max(r_ellipsoid(1),r_ellipsoid(2),r_ellipsoid(3))
    vol_sphere = 4./3.*pi*r_ellipsoid(1)*r_ellipsoid(2)*r_ellipsoid(3)
 endif
 dens_sphere = totmass_sphere/vol_sphere
 dens_sphere_cgs = dens_sphere*unit_density
 t_ff = sqrt(3.*pi/(32.*dens_sphere))

 !
 ! setup particles in the sphere
 !
 lattice = 'closepacked'
 if (is_sphere) then
    if (gauss_density) then
       density_func => gauss_density_func
       call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh,nptot=npart_total,&
                       rhofunc=density_func,exactN=.true.,np_requested=npart_sphere,mask=i_belong)
    else
       call set_sphere(trim(lattice),id,master,0.,r_sphere,psep,hfact,npart,xyzh,nptot=npart_total,&
                       exactN=.true.,np_requested=npart_sphere,mask=i_belong)
    endif
    print "(a)",' Initialised sphere'
 else
    call set_ellipse(trim(lattice),id,master,r_ellipsoid,vol_sphere,psep,hfact,xyzh,npart,&
                     nptot=npart_total,exactN=.true.,np_requested=npart_sphere,mask=i_belong)
    print "(a)",' Initialised ellipsoid'
 endif

 print*,'requested number of particles:   ',npart_sphere
 print*,'final number of particles:       ',npart,npart_total

 !
 ! Particle mass
 !
 massoftype(igas) = totmass_sphere/npart

 !
 ! set particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 do ip = 1,npartoftype(igas)
    call set_particle_type(ip,igas)
 enddo

 !
 ! Temperature and sound speed
 !
 if (apply_cooling) then
    !- Estiamte temp_sphere based on mean density
    temp_sphere = get_eqtemp_from_rho(dens_sphere_cgs)
 else
    !- Set gas internal energy with input temp_sphere
    u_sphere = kboltz * temp_sphere / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
    do ip = 1,npart
       vxyzu(4,ip) = u_sphere
    enddo
 endif
 cs_sphere_cgs = sqrt(temp_sphere*(gamma*kboltz)/(gmw*mass_proton_cgs))
 cs_sphere     = cs_sphere_cgs/unit_velocity

 !
 ! Polytropic constant
 !
 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs_sphere**2
 else
    polyk = 0.
 endif

 !
 ! Turbulent velocity field
 !
 vxyzu = 0.
 if (rms_mach > 0.) then
    call getcwd(cwd)

    filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
    filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
    filez  = find_phantom_datafile(filevz,'velfield_sphng_small')
    !- To convert endian for different vfield files: setenv GFORTRAN_CONVERT_UNIT big/small_endian

    turbboxsize = 1.1*rmax
    call set_velfield_from_cubes(xyzh(:,1:npart),vxyzu(:,:npart),npart, &
                                 filex,filey,filez,1.,turbboxsize,.false.,ierr)
    if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

    rmsmach = 0.0
    print*, 'Turbulence being set by user'
    do ip = 1,npart
       v2i     = dot_product(vxyzu(1:3,ip),vxyzu(1:3,ip))
       rmsmach = rmsmach + v2i/cs_sphere**2
    enddo
    rmsmach = sqrt(rmsmach/npart)
    if (rmsmach > 0.) then
       turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
     else
       turbfac = 0.
    endif
    do ip = 1,npart
       vxyzu(1:3,ip) = turbfac*vxyzu(1:3,ip)
    enddo

    !- Calculate rms-velocity
    v2_sum = 0.
    do ip = 1,npart
       v2_sum = v2_sum + mag2(vxyzu(1:3,ip))
    enddo
    v_rms = sqrt(1./npart*v2_sum)
    v_rms_kms = v_rms*unit_velocity*1E-5

 endif

 !
 ! Manually place a sink as feedback source
 !
 if (place_sink_in_setup) then
    nptmass                      = 1
    xyzmh_ptmass(:,:)            = 0.
    xyzmh_ptmass(1:3,nptmass)    = 0.
    xyzmh_ptmass(4,nptmass)      = 40.*solarm/umass
    xyzmh_ptmass(ihacc,nptmass)  = 0.005*pc/udist
    xyzmh_ptmass(ihsoft,nptmass) = 0.005*pc/udist
    vxyz_ptmass                  = 0.
 endif

  !
  ! Set default runtime parameters if .in file does not exist
  !
  infilename=trim(fileprefix)//'.in'
  inquire(file=infilename,exist=in_iexist)
  if (.not.in_iexist) then
     tmax      = 3.15360E15/utime ! 1E2 Myr
     dtmax     = 3.15360E9/utime  ! 1E-4 Myr
     nout      = 100
     nfulldump = 1
     nmaxdumps = 1000
     dtwallmax = 1800.  ! s
     iverbose  = 1

     ieos      = 2    ! adiabatic eos with P = (gamma-1)*rho*u
     gmw       = 1.29
     if (apply_cooling) then
        icooling = 7
        Tfloor   = 3.
        ufloor   = kboltz*Tfloor/(gmw*mass_proton_cgs*(gamma-1.))/unit_ergg
     else
        icooling = 0
     endif
     ipdv_heating   = 1
     ishock_heating = 1

     !
     ! Sinks settings
     !
     if (make_sinks) then
        icreate_sinks = 1
        rho_crit_cgs  = 1.E-16          ! density above which sink particles are created
        h_acc_cgs     = 0.005*pc        ! accretion radius for new sink particles
        h_soft_sinksink_cgs = 0.005*pc  ! softening length between sink particles
        h_soft_sinkgas_cgs  = 0.        ! softening length for new sink particles

        !- Check if rho_crit is sensible [10^5 times the initial MC density (Bate et al 1995)]
        rho_crit_cgs_recomm = 1.d5*dens_sphere_cgs
        if (abs(log10(rho_crit_cgs_recomm)-log10(rho_crit_cgs))/log10(rho_crit_cgs_recomm) > 0.1) then
           print*,'Recommend setting rho_crit_cgs to ',rho_crit_cgs_recomm,' instead'
           proceed = 'n'
           call prompt('Do you wish to continue?',proceed)
           if (trim(adjustl(proceed)) == 'n') then
              stop 'Edit rho_crit_cgs and rerun'
           elseif (trim(adjustl(proceed)) /= 'y') then
              stop 'Invalid input'
           endif
        endif

        !- convert to code units
        h_acc           = h_acc_cgs/udist
        r_crit          = 2.*h_acc
        h_soft_sinksink = h_soft_sinksink_cgs/udist
        h_soft_sinkgas  = h_soft_sinkgas_cgs/udist

     else
        icreate_sinks = 0
     endif

     !
     ! Photoionization settings
     !
     inject_rad  = .false.
     sink_ionsrc = .false.

     monochrom_source  = .false.
     fix_temp_hii      = .false.
     treat_Rtype_phase = .true.

     photoionize_tree  = .true.
     tree_accuracy_cmi = 0.3
     nHlimit_fac       = 50
     rcut_opennode_cgs = 6.0*pc
     rcut_leafpart_cgs = 4.0*pc
     delta_rcut_cgs    = 0.1*pc

     !
     ! Supernova settings
     !
     inject_sn = .false.
     sink_progenitor = .false.

  endif

  !
  ! Calculate the approx number of stars that will form
  !
  jeans_mass_cgs = (5.*Rg*temp_sphere/(2.*gg*gmw))**(3./2.) * (4./3.*pi*dens_sphere_cgs)**(-1./2.)
  jeans_mass = jeans_mass_cgs/umass

  !
  ! Write summary
  !
  print*,'-Cloud-'
  print*,'total mass        ',totmass_sphere,mass_unit
  print*,'Jeans mass        ',jeans_mass,mass_unit
  if (is_sphere) then
     print*,'radius            ',r_sphere,dist_unit
  else
     print*,'semi-axis a       ',r_ellipsoid(1),dist_unit
     print*,'semi-axis b       ',r_ellipsoid(2),dist_unit
     print*,'semi-axis c       ',r_ellipsoid(3),dist_unit
  endif
  print*,'volume            ',vol_sphere*udist**3*3.4E-56,'pc^3'
  if (gauss_density) then
     print*,'mean density      ',dens_sphere_cgs,'g/cm^3'
     print*,'est. free-fall time ',t_ff*utime/(1E6*365*24*60*60),'Myr'
  else
     print*,'density           ',dens_sphere_cgs,'g/cm^3'
     print*,'free-fall time    ',t_ff*utime/(1E6*365*24*60*60),'Myr'
  endif
  if (apply_cooling) then
     print*,'est. temperature  ',temp_sphere,'K'
  else
     print*,'temperature       ',temp_sphere,'K'
     print*,'internal energy   ',u_sphere*unit_ergg,'erg/g'
  endif
  print*,'sound speed       ',cs_sphere_cgs*1E-5,'km/s'
  print*,'v_rms             ',v_rms_kms,'km/s'
  print*,''
  print*,'-Particles-'
  print*,'total number      ',npart
  print*,'particle mass     ',massoftype(igas),mass_unit

  proceed = 'y'
  call prompt('Do you wish to continue?',proceed)
  if (trim(adjustl(proceed)) == 'n') then
     stop
  elseif (trim(adjustl(proceed)) /= 'y') then
     stop 'Invalid input'
  else
     open(2022,file='cloud_particles_info.txt')
     write(2022,*) '-Cloud-'
     write(2022,*) 'total mass        ',totmass_sphere,mass_unit
     if (gauss_density) then
        write(2022,*) 'est. Jeans mass   ',jeans_mass,mass_unit
     else
        write(2022,*) 'Jeans mass        ',jeans_mass,mass_unit
     endif
     if (is_sphere) then
        write(2022,*) 'radius            ',r_sphere,dist_unit
     else
        write(2022,*) 'semi-axis a       ',r_ellipsoid(1),dist_unit
        write(2022,*) 'semi-axis b       ',r_ellipsoid(2),dist_unit
        write(2022,*) 'semi-axis c       ',r_ellipsoid(3),dist_unit
     endif
     write(2022,*) 'volume            ',vol_sphere*udist**3*3.4E-56,'pc^3'
     if (gauss_density) then
        write(2022,*) 'mean density      ',dens_sphere_cgs,'g/cm^3'
        write(2022,*) 'est. free-fall time ',t_ff*utime/(1E6*365*24*60*60),'Myr'
     else
        write(2022,*) 'density           ',dens_sphere_cgs,'g/cm^3'
        write(2022,*) 'free-fall time    ',t_ff*utime/(1E6*365*24*60*60),'Myr'
     endif
     if (apply_cooling) then
        write(2022,*) 'est. temperature  ',temp_sphere,'K'
        write(2022,*) 'est. sound speed  ',cs_sphere_cgs*1E-5,'km/s'
     else
        write(2022,*) 'temperature       ',temp_sphere,'K'
        write(2022,*) 'sound speed       ',cs_sphere_cgs*1E-5,'km/s'
     endif

     write(2022,*) 'Mach number       ',rms_mach
     write(2022,*) 'v_rms             ',v_rms_kms,'km/s'
     write(2022,*) ''
     write(2022,*) '-Particles-'
     write(2022,*) 'total number      ',npart
     write(2022,*) 'particle mass     ',massoftype(igas),mass_unit

     close(2022)
  endif

end subroutine setpart

!
! Set Gaussian density profile
!
real function gauss_density_func(r)
 use physcon,  only:pi
 real, intent(in) :: r
 real :: sigma

 sigma = 0.67*r_sphere  !- rho at center is 3x of edge
 gauss_density_func = exp(-r**2/(2*sigma**2))

end function gauss_density_func


real function mag2(vec)
 real, intent(in) :: vec(3)

 mag2 = vec(1)**2 + vec(2)**2 + vec(3)**2

end function mag2

!
! Estimate temperature with given rho from cooling curve
!
real function get_eqtemp_from_rho(rho_cgs)
 use datafiles, only:find_phantom_datafile
 real,   intent(in) :: rho_cgs
 integer, parameter :: maxrow = 1000
 integer, parameter :: maxcol = 4
 integer :: irow,icol,iclosest
 real    :: rhoteq(maxcol,maxrow)

 open(3000,file=find_phantom_datafile('rhoteq_table.txt','coolcurve_rho_temp'))
 do irow = 1,maxrow
    read(3000,*) (rhoteq(icol,irow), icol=1,maxcol)
 enddo
 close(3000)

 iclosest = minloc(abs(rhoteq(1,:)-rho_cgs),1)
 get_eqtemp_from_rho = rhoteq(3,iclosest)

end function get_eqtemp_from_rho


!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for sphere setup routines'
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)

 write(iunit,"(/,a)") '# options for shape'
 call write_inopt(is_sphere,'is_sphere','the intent is to set up a sphere',iunit)

 if (is_sphere) then
    write(iunit,"(/,a)") '# options for sphere'
    call write_inopt(npart_sphere,'npart_sphere','requested number of particles in sphere',iunit)
    call write_inopt(r_sphere,'r_sphere','radius of sphere in code units',iunit)
    call write_inopt(totmass_sphere,'totmass_sphere','mass of sphere in code units',iunit)
    call write_inopt(gauss_density,'gauss_density','the intent is to centrally condense the sphere',iunit)
 else
    write(iunit,"(/,a)") '# options for ellipsoid'
    call write_inopt(npart_sphere,'npart_sphere','requested number of particles in ellipsoid',iunit)
    call write_inopt(r_ellipsoid(1),'r_ellipsoid(1)','semi-axis {a} of ellipsoid in code units',iunit)
    call write_inopt(r_ellipsoid(2),'r_ellipsoid(2)','semi-axis {b} of ellipsoid in code units',iunit)
    call write_inopt(r_ellipsoid(3),'r_ellipsoid(3)','semi-axis {c} of ellipsoid in code units',iunit)
    call write_inopt(totmass_sphere,'totmass_sphere','mass of ellipsoid in code units',iunit)
 endif

 call write_inopt(apply_cooling,'apply_cooling','the intent is to apply cooling',iunit)
 call write_inopt(temp_sphere,'temp_sphere','temperature of cloud in K',iunit)

 call write_inopt(rms_mach,'rms_mach','turbulent rms mach number',iunit)

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
 use unifdis,      only: is_valid_lattice
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)

 call read_inopt(is_sphere,'is_sphere',db,ierr)
 call read_inopt(npart_sphere,'npart_sphere',db,ierr)
 call read_inopt(totmass_sphere,'totmass_sphere',db,ierr)

 if (is_sphere) then
   call read_inopt(r_sphere,'r_sphere',db,ierr)
   call read_inopt(gauss_density,'gauss_density',db,ierr)
 else
    call read_inopt(r_ellipsoid(1),'r_ellipsoid(1)',db,ierr)
    call read_inopt(r_ellipsoid(2),'r_ellipsoid(2)',db,ierr)
    call read_inopt(r_ellipsoid(3),'r_ellipsoid(3)',db,ierr)
 endif

 call read_inopt(apply_cooling,'apply_cooling',db,ierr)
 call read_inopt(temp_sphere,'temp_sphere',db,ierr)

 call read_inopt(rms_mach,'rms_mach',db,ierr)

 call read_inopt(make_sinks,'make_sinks',db,ierr)

 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_sphere_sphngsne','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_sphere_sphngsne','length unit not recognised')
    ierr = ierr + 1
 endif

 if (ierr > 0) then
    print "(1x,a,i2,a)",'setup_sphere_rhdsne: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif

end subroutine read_setupfile

end module setup
