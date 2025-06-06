!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module readwrite_infile
!
! This module contains all routines required for
!  reading and writing of input file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - C_cour             : *Courant number*
!   - C_force            : *dt_force number*
!   - alpha              : *art. viscosity parameter*
!   - alphaB             : *art. resistivity parameter*
!   - alphamax           : *MAXIMUM art. viscosity parameter*
!   - alphau             : *art. conductivity parameter*
!   - avdecayconst       : *decay time constant for viscosity switches*
!   - beta               : *beta viscosity*
!   - bulkvisc           : *magnitude of bulk viscosity*
!   - calc_erot          : *include E_rot in the ev_file*
!   - dtmax              : *time between dumps*
!   - dtmax_dratio       : *dynamic dtmax: density ratio controlling decrease (<=0 to ignore)*
!   - dtmax_max          : *dynamic dtmax: maximum allowed dtmax (=dtmax if <= 0)*
!   - dtmax_min          : *dynamic dtmax: minimum allowed dtmax*
!   - dtwallmax          : *maximum wall time between dumps (hhh:mm, 000:00=ignore)*
!   - dumpfile           : *dump file to start from*
!   - flux_limiter       : *limit radiation flux*
!   - hdivbbmax_max      : *max factor to decrease cleaning timestep propto B/(h|divB|)*
!   - hfact              : *h in units of particle spacing [h = hfact(m/rho)^(1/3)]*
!   - iopacity_type      : *opacity method (0=inf,1=mesa)*
!   - ipdv_heating       : *heating from PdV work (0=off, 1=on)*
!   - irealvisc          : *physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)*
!   - iresistive_heating : *resistive heating (0=off, 1=on)*
!   - ishock_heating     : *shock heating (0=off, 1=on)*
!   - iverbose           : *verboseness of log (-1=quiet 0=default 1=allsteps 2=debug 5=max)*
!   - logfile            : *file to which output is directed*
!   - nfulldump          : *full dump every n dumps*
!   - nmax               : *maximum number of timesteps (0=just get derivs and stop)*
!   - nmaxdumps          : *stop after n full dumps (-ve=ignore)*
!   - nout               : *write dumpfile every n dtmax (-ve=ignore)*
!   - overcleanfac       : *factor to increase cleaning speed (decreases time step)*
!   - psidecayfac        : *div B diffusion parameter*
!   - ptol               : *tolerance on pmom iterations*
!   - rhofinal_cgs       : *maximum allowed density (cgs) (<=0 to ignore)*
!   - rkill              : *deactivate particles outside this radius (<0 is off)*
!   - shearparam         : *magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)*
!   - tmax               : *end time*
!   - tolh               : *tolerance on h-rho iterations*
!   - tolv               : *tolerance on v iterations in timestepping*
!   - twallmax           : *maximum wall time (hhh:mm, 000:00=ignore)*
!   - use_mcfost         : *use the mcfost library*
!   - xtol               : *tolerance on xyz iterations*
!
! :Dependencies: cooling, damping, dim, dust, dust_formation, eos,
!   externalforces, forcing, growth, infile_utils, inject, io, linklist,
!   metric, nicil_sup, options, part, photoevap, ptmass, ptmass_radiation,
!   timestep, viscosity, photoionize_cmi
!
 use timestep,  only:dtmax_dratio,dtmax_max,dtmax_min
 use options,   only:nfulldump,nmaxdumps,twallmax,iexternalforce,idamp,tolh, &
                     alpha,alphau,alphaB,beta,avdecayconst,damp,rkill, &
                     ipdv_heating,ishock_heating,iresistive_heating, &
                     icooling,psidecayfac,overcleanfac,hdivbbmax_max,alphamax,calc_erot,rhofinal_cgs, &
                     use_mcfost, use_Voronoi_limits_file, Voronoi_limits_file, use_mcfost_stellar_parameters,&
                     exchange_radiation_energy,limit_radiation_flux,iopacity_type, mcfost_computes_Lacc
 use timestep,  only:dtwallmax,tolv,xtol,ptol
 use viscosity, only:irealvisc,shearparam,bulkvisc
 use part,      only:hfact
 use io,        only:iverbose
 use dim,       only:do_radiation
 implicit none
 logical :: incl_runtime2 = .false.
 character(len=80), parameter, public :: &
    modid="$Id$"

contains

!-----------------------------------------------------------------
!+
!  writes an input file
!+
!-----------------------------------------------------------------
subroutine write_infile(infile,logfile,evfile,dumpfile,iwritein,iprint)
 use timestep,        only:tmax,dtmax,nmax,nout,C_cour,C_force
 use io,              only:fatal
 use infile_utils,    only:write_inopt
#ifdef DRIVING
 use forcing,         only:write_options_forcing
#endif
 use externalforces,  only:write_options_externalforces
 use damping,         only:write_options_damping
 use linklist,        only:write_inopts_link
#ifdef DUST
 use dust,            only:write_options_dust
#ifdef DUSTGROWTH
 use growth,          only:write_options_growth
#endif
#endif
#ifdef PHOTO
 use photoevap,       only:write_options_photoevap
#endif
#ifdef PHOTOION
 use photoionize_cmi, only:write_options_photoionize
#endif
#ifdef INJECT_PARTICLES
 use inject,          only:write_options_inject
 use dust_formation,  only:write_options_dust_formation
#endif
#ifdef NONIDEALMHD
 use nicil_sup,       only:write_options_nicil
#endif
#ifdef GR
 use metric,          only:write_options_metric
#endif
 use eos,             only:write_options_eos,ieos
 use ptmass,          only:write_options_ptmass
 use ptmass_radiation,only:write_options_ptmass_radiation
 use cooling,         only:write_options_cooling
 use dim,             only:maxvxyzu,maxptmass,gravity,sink_radiation,gr
 use part,            only:h2chemistry,maxp,mhd,maxalpha,nptmass
 character(len=*), intent(in) :: infile,logfile,evfile,dumpfile
 integer,          intent(in) :: iwritein,iprint
 integer                      :: ierr
 character(len=10)            :: startdate,starttime

 if (iwritein /= iprint) then
    open(unit=iwritein,iostat=ierr,file=infile,status='replace',form='formatted')
    if (ierr /= 0) call fatal('write_infile','error creating input file, exiting...')
 endif

 if (iwritein /= iprint) then
!
! get date and time to timestamp the header
!
    call date_and_time(startdate,starttime)
    startdate = startdate(7:8)//'/'//startdate(5:6)//'/'//startdate(1:4)
    starttime = starttime(1:2)//':'//starttime(3:4)//':'//starttime(5:)
!
! write header to the file
!
    write(iwritein,"(a)") '# Runtime options file for Phantom, written '//startdate//' '//starttime
    write(iwritein,"(a)") '# Options not present assume their default values'
    write(iwritein,"(a)") '# This file is updated automatically after a full dump'
 endif

 write(iwritein,"(/,a)") '# job name'
 call write_inopt(trim(logfile),'logfile','file to which output is directed',iwritein)
 call write_inopt(trim(dumpfile),'dumpfile','dump file to start from',iwritein)

 write(iwritein,"(/,a)") '# options controlling run time and input/output'
 call write_inopt(tmax,'tmax','end time',iwritein)
 call write_inopt(dtmax,'dtmax','time between dumps',iwritein)
 call write_inopt(nmax,'nmax','maximum number of timesteps (0=just get derivs and stop)',iwritein)
 call write_inopt(nout,'nout','write dumpfile every n dtmax (-ve=ignore)',iwritein)
 call write_inopt(nmaxdumps,'nmaxdumps','stop after n full dumps (-ve=ignore)',iwritein)
 call write_inopt(real(twallmax),'twallmax','maximum wall time (hhh:mm, 000:00=ignore)',iwritein,time=.true.)
 call write_inopt(real(dtwallmax),'dtwallmax','maximum wall time between dumps (hhh:mm, 000:00=ignore)',iwritein,time=.true.)
 call write_inopt(nfulldump,'nfulldump','full dump every n dumps',iwritein)
 call write_inopt(iverbose,'iverbose','verboseness of log (-1=quiet 0=default 1=allsteps 2=debug 5=max)',iwritein)

 if (incl_runtime2 .or. rhofinal_cgs > 0.0 .or. dtmax_dratio > 1.0 .or. calc_erot) then
    write(iwritein,"(/,a)") '# options controlling run time and input/output: supplementary features'
    call write_inopt(rhofinal_cgs,'rhofinal_cgs','maximum allowed density (cgs) (<=0 to ignore)',iwritein)
    call write_inopt(dtmax_dratio,'dtmax_dratio','dynamic dtmax: density ratio controlling decrease (<=0 to ignore)',iwritein)
    call write_inopt(dtmax_max,'dtmax_max','dynamic dtmax: maximum allowed dtmax (=dtmax if <= 0)',iwritein)
    call write_inopt(dtmax_min,'dtmax_min','dynamic dtmax: minimum allowed dtmax',iwritein)
    call write_inopt(calc_erot,'calc_erot','include E_rot in the ev_file',iwritein)
 endif

 write(iwritein,"(/,a)") '# options controlling accuracy'
 call write_inopt(C_cour,'C_cour','Courant number',iwritein)
 call write_inopt(C_force,'C_force','dt_force number',iwritein)
 call write_inopt(tolv,'tolv','tolerance on v iterations in timestepping',iwritein,exp=.true.)
 call write_inopt(hfact,'hfact','h in units of particle spacing [h = hfact(m/rho)^(1/3)]',iwritein)
 call write_inopt(tolh,'tolh','tolerance on h-rho iterations',iwritein,exp=.true.)
 if (gr) then
    call write_inopt(xtol,'xtol','tolerance on xyz iterations',iwritein)
    call write_inopt(ptol,'ptol','tolerance on pmom iterations',iwritein)
 endif

 call write_inopts_link(iwritein)

 write(iwritein,"(/,a)") '# options controlling hydrodynamics, artificial dissipation'
 if (maxalpha==maxp) then
    call write_inopt(alpha,'alpha','MINIMUM art. viscosity parameter',iwritein)
    call write_inopt(alphamax,'alphamax','MAXIMUM art. viscosity parameter',iwritein)
 else
    call write_inopt(alpha,'alpha','art. viscosity parameter',iwritein)
 endif
 if (maxvxyzu >= 4) then
    call write_inopt(alphau,'alphau','art. conductivity parameter',iwritein)
 endif
 if (mhd) then
    call write_inopt(alphaB,'alphaB','art. resistivity parameter',iwritein)
    call write_inopt(psidecayfac,'psidecayfac','div B diffusion parameter',iwritein)
    call write_inopt(overcleanfac,'overcleanfac','factor to increase cleaning speed (decreases time step)',iwritein)
    call write_inopt(hdivbbmax_max,'hdivbbmax_max','max factor to decrease cleaning timestep propto B/(h|divB|)',iwritein)
 endif
 call write_inopt(beta,'beta','beta viscosity',iwritein)
 call write_inopt(avdecayconst,'avdecayconst','decay time constant for viscosity switches',iwritein)

 call write_options_damping(iwritein,idamp)

 !
 ! thermodynamics
 !
 call write_options_eos(iwritein)
 if (maxvxyzu >= 4 .and. (ieos==2 .or. ieos==10 .or. ieos==15 .or. ieos==12 .or. ieos==16 .or. ieos==19) ) then
    call write_inopt(ipdv_heating,'ipdv_heating','heating from PdV work (0=off, 1=on)',iwritein)
    call write_inopt(ishock_heating,'ishock_heating','shock heating (0=off, 1=on)',iwritein)
    if (mhd) then
       call write_inopt(iresistive_heating,'iresistive_heating','resistive heating (0=off, 1=on)',iwritein)
    endif
 endif

 if (maxvxyzu >= 4) call write_options_cooling(iwritein)

#ifdef MCFOST
 call write_inopt(use_mcfost,'use_mcfost','use the mcfost library',iwritein)
 if (use_Voronoi_limits_file) call write_inopt(Voronoi_limits_file,'Voronoi_limits_file',&
      'Limit file for the Voronoi tesselation',iwritein)
 call write_inopt(use_mcfost_stellar_parameters,'use_mcfost_stars',&
      'Fix the stellar parameters to mcfost values or update using sink mass',iwritein)
 call write_inopt(mcfost_computes_Lacc,'mcfost_computes_Lacc',&
      'Should mcfost compute the accretion luminosity',iwritein)
#endif

 ! only write sink options if they are used, or if self-gravity is on
 if (nptmass > 0 .or. gravity) call write_options_ptmass(iwritein)

 call write_options_externalforces(iwritein,iexternalforce)

 write(iwritein,"(/,a)") '# options controlling physical viscosity'
 call write_inopt(irealvisc,'irealvisc','physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)',iwritein)
 call write_inopt(shearparam,'shearparam','magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)',iwritein)
 call write_inopt(bulkvisc,'bulkvisc','magnitude of bulk viscosity',iwritein)

#ifdef DRIVING
 call write_options_forcing(iwritein)
#endif

#ifdef DUST
 call write_options_dust(iwritein)
#ifdef DUSTGROWTH
 call write_options_growth(iwritein)
#endif
#endif

#ifdef PHOTO
 call write_options_photoevap(iwritein)
#endif

#ifdef PHOTOION
 call write_options_photoionize(iwritein)
#endif

#ifdef INJECT_PARTICLES
 call write_options_inject(iwritein)
 call write_options_dust_formation(iwritein)
#endif

 write(iwritein,"(/,a)") '# options for injecting/removing particles'
 call write_inopt(rkill,'rkill','deactivate particles outside this radius (<0 is off)',iwritein)

 if (sink_radiation) then
    write(iwritein,"(/,a)") '# options controling radiation pressure from sink particles'
    call write_options_ptmass_radiation(iwritein)
 endif
#ifdef NONIDEALMHD
 call write_options_nicil(iwritein)
#endif

 if (do_radiation) then
    write(iwritein,"(/,a)") '# options for radiation'
    call write_inopt(exchange_radiation_energy,'gas-rad_exchange','exchange energy between gas and radiation',iwritein)
    call write_inopt(limit_radiation_flux,'flux_limiter','limit radiation flux',iwritein)
    call write_inopt(iopacity_type,'iopacity_type','opacity method (0=inf,1=mesa)',iwritein)
 endif
#ifdef GR
 call write_options_metric(iwritein)
#endif

 if (iwritein /= iprint) close(unit=iwritein)
 if (iwritein /= iprint) write(iprint,"(/,a)") ' input file '//trim(infile)//' written successfully.'

 return
end subroutine write_infile

!-----------------------------------------------------------------
!+
!  reads parameters for the run from the input file
!+
!-----------------------------------------------------------------
subroutine read_infile(infile,logfile,evfile,dumpfile)
 use dim,             only:maxvxyzu,maxptmass,gravity,sink_radiation
 use timestep,        only:tmax,dtmax,nmax,nout,C_cour,C_force
 use eos,             only:use_entropy,read_options_eos,ieos
 use io,              only:ireadin,iwritein,iprint,warn,die,error,fatal,id,master
 use infile_utils,    only:read_next_inopt,contains_loop,write_infile_series
#ifdef DRIVING
 use forcing,         only:read_options_forcing,write_options_forcing
#endif
 use externalforces,  only:read_options_externalforces
 use linklist,        only:read_inopts_link
#ifdef DUST
 use dust,            only:read_options_dust
#ifdef DUSTGROWTH
 use growth,          only:read_options_growth
#endif
#endif
#ifdef GR
 use metric,        only:read_options_metric
#endif
#ifdef PHOTO
 use photoevap,       only:read_options_photoevap
#endif
#ifdef PHOTOION
 use photoionize_cmi, only:read_options_photoionize
#endif
#ifdef INJECT_PARTICLES
 use inject,          only:read_options_inject
 use dust_formation,  only:read_options_dust_formation
#ifdef wind
 use dust_formation,  only:idust_opacity
#endif
#endif
#ifdef NONIDEALMHD
 use nicil_sup,       only:read_options_nicil
#endif
 use part,            only:mhd,nptmass
 use cooling,         only:read_options_cooling
 use ptmass,          only:read_options_ptmass
 use ptmass_radiation,only:read_options_ptmass_radiation
#ifdef WIND
 use ptmass_radiation,only:isink_radiation,alpha_rad
#endif
 use damping,         only:read_options_damping
 character(len=*), parameter   :: label = 'read_infile'
 character(len=*), intent(in)  :: infile
 character(len=*), intent(out) :: logfile,evfile,dumpfile
 character(len=len(infile)+4)  :: infilenew
 character(len=10) :: cline
 character(len=20) :: name
 character(len=120) :: valstring
 integer :: ierr,ireaderr,line,idot,ngot,nlinesread
 real    :: ratio
 logical :: imatch,igotallrequired,igotallturb,igotalllink,igotloops
 logical :: igotallbowen,igotallcooling,igotalldust,igotallextern,igotallinject,igotallgrowth
 logical :: igotallionise,igotallnonideal,igotalleos,igotallptmass,igotallphoto,igotalldamping
 logical :: igotallprad,igotalldustform,igotallphotoion
 integer, parameter :: nrequired = 1

 ireaderr = 0
 ierr     = 0
 line     = 0
 idot     = index(infile,'.in') - 1
 if (idot <= 1) idot = len_trim(infile)
 logfile    = infile(1:idot)//'01.log'
 dumpfile   = infile(1:idot)//'_00000.tmp'
 ngot            = 0
 igotallturb     = .true.
 igotalldust     = .true.
 igotallgrowth   = .true.
 igotallphoto    = .true.
 igotallphotoion = .true.
 igotalllink     = .true.
 igotallextern   = .true.
 igotallinject   = .true.
 igotalleos      = .true.
 igotallcooling  = .true.
 igotalldamping  = .true.
 igotloops       = .false.
 igotallionise   = .true.
 igotallnonideal = .true.
 igotallbowen    = .true.
 igotallptmass   = .true.
 igotallprad     = .true.
 igotalldustform = .true.
 use_Voronoi_limits_file = .false.

 open(unit=ireadin,err=999,file=infile,status='old',form='formatted')
 do while (ireaderr == 0)
    call read_next_inopt(name,valstring,ireadin,ireaderr,nlinesread)
    if (contains_loop(valstring)) igotloops = .true.
    line = line + nlinesread

    !print*,'name: '//trim(name),' value: '//trim(valstring)
    select case(trim(name))
    case('logfile')
       logfile = trim(valstring(1:min(len(logfile),len(valstring))))
       idot = index(logfile,'.log') - 1
       if (idot <= 1) idot = len_trim(logfile)
       evfile  = logfile(1:idot)//'.ev'
    case('dumpfile')
       dumpfile = trim(valstring(1:min(len(dumpfile),len(valstring))))
       ngot = ngot + 1
    case('tmax')
       read(valstring,*,iostat=ierr) tmax
    case('dtmax')
       read(valstring,*,iostat=ierr) dtmax
    case('nmax')
       read(valstring,*,iostat=ierr) nmax
    case('nout')
       read(valstring,*,iostat=ierr) nout
    case('nmaxdumps')
       read(valstring,*,iostat=ierr) nmaxdumps
    case('twallmax')
       read(valstring,*,iostat=ierr) twallmax
    case('dtwallmax')
       read(valstring,*,iostat=ierr) dtwallmax
    case('iverbose')
       read(valstring,*,iostat=ierr) iverbose
    case('rhofinal_cgs')
       read(valstring,*,iostat=ierr) rhofinal_cgs
       incl_runtime2 = .true.
    case('calc_erot')
       read(valstring,*,iostat=ierr) calc_erot
       incl_runtime2 = .true.
    case('dtmax_dratio')
       read(valstring,*,iostat=ierr) dtmax_dratio
       incl_runtime2 = .true.
    case('dtmax_max')
       read(valstring,*,iostat=ierr) dtmax_max
       if (dtmax_max <= 0.0) dtmax_max = dtmax
       ! to prevent comparison errors from round-off
       ratio = dtmax_max/dtmax
       ratio = int(ratio+0.5)+0.0001
       dtmax_max = dtmax*ratio
    case('dtmax_min')
       read(valstring,*,iostat=ierr) dtmax_min
       ! to prevent comparison errors from round-off
       if (dtmax_min > epsilon(dtmax_min)) then
          ratio = dtmax/dtmax_min
          ratio = int(ratio+0.5)+0.0001
          dtmax_min = dtmax/ratio
       endif
    case('C_cour')
       read(valstring,*,iostat=ierr) C_cour
    case('C_force')
       read(valstring,*,iostat=ierr) C_force
    case('tolv')
       read(valstring,*,iostat=ierr) tolv
    case('xtol')
       read(valstring,*,iostat=ierr) xtol
    case('ptol')
       read(valstring,*,iostat=ierr) ptol
    case('hfact')
       read(valstring,*,iostat=ierr) hfact
    case('tolh')
       read(valstring,*,iostat=ierr) tolh
    case('rkill')
       read(valstring,*,iostat=ierr) rkill
    case('nfulldump')
       read(valstring,*,iostat=ierr) nfulldump
    case('alpha')
       read(valstring,*,iostat=ierr) alpha
    case('alphamax')
       read(valstring,*,iostat=ierr) alphamax
    case('alphau')
       read(valstring,*,iostat=ierr) alphau
    case('alphaB')
       read(valstring,*,iostat=ierr) alphaB
    case('psidecayfac')
       read(valstring,*,iostat=ierr) psidecayfac
    case('overcleanfac')
       read(valstring,*,iostat=ierr) overcleanfac
    case('hdivbbmax_max')
       read(valstring,*,iostat=ierr) hdivbbmax_max
    case('beta')
       read(valstring,*,iostat=ierr) beta
    case('avdecayconst')
       read(valstring,*,iostat=ierr) avdecayconst
    case('ipdv_heating')
       read(valstring,*,iostat=ierr) ipdv_heating
    case('ishock_heating')
       read(valstring,*,iostat=ierr) ishock_heating
    case('iresistive_heating')
       read(valstring,*,iostat=ierr) iresistive_heating
    case('irealvisc')
       read(valstring,*,iostat=ierr) irealvisc
    case('shearparam')
       read(valstring,*,iostat=ierr) shearparam
    case('bulkvisc')
       read(valstring,*,iostat=ierr) bulkvisc
#ifdef MCFOST
    case('use_mcfost')
       read(valstring,*,iostat=ierr) use_mcfost
    case('Voronoi_limits_file')
       read(valstring,*,iostat=ierr) Voronoi_limits_file
       use_Voronoi_limits_file = .true.
    case('use_mcfost_stars')
       read(valstring,*,iostat=ierr) use_mcfost_stellar_parameters
    case('mcfost_computes_Lacc')
       read(valstring,*,iostat=ierr) mcfost_computes_Lacc
#endif
    case('gas-rad_exchange')
       read(valstring,*,iostat=ierr) exchange_radiation_energy
    case('flux_limiter')
       read(valstring,*,iostat=ierr) limit_radiation_flux
    case('iopacity_type')
       read(valstring,*,iostat=ierr) iopacity_type
    case default
       imatch = .false.
       if (.not.imatch) call read_options_externalforces(name,valstring,imatch,igotallextern,ierr,iexternalforce)
#ifdef DRIVING
       if (.not.imatch) call read_options_forcing(name,valstring,imatch,igotallturb,ierr)
#endif
       if (.not.imatch) call read_inopts_link(name,valstring,imatch,igotalllink,ierr)
#ifdef DUST
       !--Extract if one-fluid dust is used from the fileid
       if (.not.imatch) call read_options_dust(name,valstring,imatch,igotalldust,ierr)
#ifdef DUSTGROWTH
       if (.not.imatch) call read_options_growth(name,valstring,imatch,igotallgrowth,ierr)
#endif
#endif
#ifdef GR
       if (.not.imatch) call read_options_metric(name,valstring,imatch,igotalldust,ierr)
#endif
#ifdef PHOTO
       if (.not.imatch) call read_options_photoevap(name,valstring,imatch,igotallphoto,ierr)
#endif
#ifdef PHOTOION
       if (.not.imatch) call read_options_photoionize(name,valstring,imatch,igotallphotoion,ierr)
#endif
#ifdef INJECT_PARTICLES
       if (.not.imatch) call read_options_inject(name,valstring,imatch,igotallinject,ierr)
       if (.not.imatch) call read_options_dust_formation(name,valstring,imatch,igotalldustform,ierr)
#endif
       if (.not.imatch .and. sink_radiation) then
          call read_options_ptmass_radiation(name,valstring,imatch,igotallprad,ierr)
       endif
#ifdef NONIDEALMHD
       if (.not.imatch) call read_options_nicil(name,valstring,imatch,igotallnonideal,ierr)
#endif
       if (.not.imatch) call read_options_eos(name,valstring,imatch,igotalleos,ierr)
       if (.not.imatch .and. maxvxyzu >= 4) call read_options_cooling(name,valstring,imatch,igotallcooling,ierr)
       if (.not.imatch) call read_options_damping(name,valstring,imatch,igotalldamping,ierr,idamp)
       if (maxptmass > 0) then
          if (.not.imatch) call read_options_ptmass(name,valstring,imatch,igotallptmass,ierr)
          !
          ! read whatever sink options are present, but make them not compulsory
          ! if there are no sinks used and there is no self-gravity
          !
          if (nptmass==0 .and. .not.gravity) igotallptmass = .true.
       endif

       if (len_trim(name) /= 0 .and. .not.imatch) then
          call warn('read_infile','unknown variable '//trim(adjustl(name))// &
                     ' in input file, value = '//trim(adjustl(valstring)))
       endif
    end select
    if (ierr /= 0 .and. len_trim(name) /= 0) &
       call fatal('read_infile','error extracting '//trim(adjustl(name))//' from input file')
 enddo
 close(unit=ireadin)

 igotallrequired = (ngot  >=  nrequired) .and. igotalllink   .and. igotallbowen   .and. igotalldust &
                    .and. igotalleos    .and. igotallcooling .and. igotallextern  .and. igotallturb &
                    .and. igotallptmass .and. igotallinject  .and. igotallionise  .and. igotallnonideal &
                    .and. igotallphoto  .and. igotallgrowth  .and. igotalldamping .and. igotallprad &
                    .and. igotalldustform .and. igotallphotoion

 if (ierr /= 0 .or. ireaderr > 0 .or. .not.igotallrequired) then
    ierr = 1
    if (id==master) then
       if (igotallrequired) then
          write(cline,"(i10)") line
          call error('read_infile','error reading '//trim(infile)//' at line '//trim(adjustl(cline)))
          infilenew = trim(infile)//'.new'
       else
          call error('read_infile','input file '//trim(infile)//' is incomplete for current compilation')
          if (.not.igotalleos) write(*,*) 'missing equation of state options'
          if (.not.igotallcooling) write(*,*) 'missing cooling options'
          if (.not.igotalldamping) write(*,*) 'missing damping options'
          if (.not.igotalllink) write(*,*) 'missing link options'
          if (.not.igotallbowen) write(*,*) 'missing Bowen dust options'
          if (.not.igotalldust) write(*,*) 'missing dust options'
          if (.not.igotallgrowth) write(*,*) 'missing growth options'
          if (.not.igotallphoto) write(*,*) 'missing photoevaporation options'
          if (.not.igotallphotoion) write(*,*) 'missing photoionize_cmi options'
          if (.not.igotallextern) write(*,*) 'missing external force options'
          if (.not.igotallinject) write(*,*) 'missing inject-particle options'
          if (.not.igotallionise) write(*,*) 'missing ionisation options'
          if (.not.igotallnonideal) write(*,*) 'missing non-ideal MHD options'
          if (.not.igotallturb) write(*,*) 'missing turbulence-driving options'
          if (.not.igotallprad) write(*,*) 'missing sink particle radiation options'
          if (.not.igotallptmass) write(*,*) 'missing sink particle options'
          if (.not.igotalldustform) write(*,*) 'missing dusty wind options'
          infilenew = trim(infile)
       endif
       write(*,"(a)") ' REWRITING '//trim(infilenew)//' with all current and available options...'
       call write_infile(infilenew,logfile,evfile,dumpfile,iwritein,iprint)

       if (trim(infilenew) /= trim(infile)) then
          write(*,"(/,a)") ' useful commands to cut and paste:'
          write(*,*) ' diff '//trim(infilenew)//' '//trim(infile)
          write(*,*) ' mv '//trim(infilenew)//' '//trim(infile)
       else
          write(*,"(/,a)") ' (try again using revised input file)'
       endif
    endif
    call die
 elseif (igotloops) then
    if (id==master) then
       write(iprint,"(1x,a)") 'loops detected in input file:'
       call write_infile_series(ireadin,iwritein,infile,line,ierr)
       if (ierr /= 0) call fatal(label,'error writing input file series')
    endif
    call die
 endif
!
!--check options for possible errors
!
 if (id==master) then
    if (dtmax > tmax) call warn(label,'no output dtmax > tmax',1)
    if (nout > nmax)  call warn(label,'no output nout > nmax',1)
    if (nout==0)     call fatal(label,'nout = 0')
    if (C_cour <= 0.)  call fatal(label,'Courant number < 0')
    if (C_cour > 1.)  call fatal(label,'ridiculously big courant number!!')
    if (C_force <= 0.) call fatal(label,'bad choice for force timestep control')
    if (tolv <= 0.)    call fatal(label,'silly choice for tolv (< 0)')
    if (tolv > 1.e-1) call warn(label,'dangerously large tolerance on v iterations')
    if (xtol <= 0.)    call fatal(label,'silly choice for xtol (< 0)')
    if (xtol > 1.e-1) call warn(label,'dangerously large tolerance on xyz iterations')
    if (ptol <= 0.)    call fatal(label,'silly choice for ptol (< 0)')
    if (ptol > 1.e-1) call warn(label,'dangerously large tolerance on pmom iterations')
    if (nfulldump==0 .or. nfulldump > 10000) call fatal(label,'nfulldump = 0')
    if (nfulldump >= 50) call warn(label,'no full dumps for a long time...',1)
    if (twallmax < 0.)  call fatal(label,'invalid twallmax (use 000:00 to ignore)')
    if (dtwallmax < 0.) call fatal(label,'invalid dtwallmax (use 000:00 to ignore)')
    if (hfact < 1. .or. hfact > 5.) &
                         call warn(label,'ridiculous choice of hfact',4)
    if (tolh > 1.e-3)   call warn(label,'tolh is quite large!',2)
    if (tolh < epsilon(tolh)) call fatal(label,'tolh too small to ever converge')
    !if (damp < 0.)     call fatal(label,'damping < 0')
    !if (damp > 1.)     call warn(label,'damping ridiculously big')
    if (alpha < 0.)    call fatal(label,'stupid choice of alpha')
    if (alphau < 0.)   call fatal(label,'stupid choice of alphau')
    if (alphau > tiny(alphau) .and. use_entropy) &
                        call fatal(label,'cannot use thermal conductivity if evolving entropy')
    if (alpha > 10.)   call warn(label,'very large alpha, need to change timestep',2)
    if (alphau > 10.)  call warn(label,'very large alphau, check timestep',3)
    if (alphamax < tiny(alphamax)) call warn(label,'alphamax = 0 means no shock viscosity',2)
    if (alphamax < 1.) call warn(label,'alphamax < 1 is dangerous if there are shocks: don''t publish crap',2)
    if (alphamax < 0. .or. alphamax > 100.) call fatal(label,'stupid value for alphamax (generally 0.0-1.0)')
    if (mhd) then
       if (alphaB < 0.)   call warn(label,'stupid choice of alphaB',4)
       if (alphaB > 10.)  call warn(label,'very large alphaB, check timestep',3)
       if (alphaB < 1.)   call warn(label,'alphaB < 1 is not recommended, please don''t publish rubbish',2)
       if (psidecayfac < 0.) call fatal(label,'stupid value for psidecayfac')
       if (psidecayfac > 2.) call warn(label,'psidecayfac set outside recommended range (0.1-2.0)')
       if (overcleanfac < 1.0) call warn(label,'overcleanfac less than 1')
       if (hdivbbmax_max < overcleanfac) then
          call warn(label,'Resetting hdivbbmax_max = overcleanfac')
          hdivbbmax_max = overcleanfac
       endif
    endif
    if (beta < 0.)     call fatal(label,'beta < 0')
    if (beta > 4.)     call warn(label,'very high beta viscosity set')
#ifndef MCFOST
    if (maxvxyzu >= 4 .and. (ieos /= 2 .and. ieos /= 4 .and. ieos /= 10 .and. ieos /=11 .and. &
                             ieos /=12 .and. ieos /= 15 .and. ieos /= 16 .and. ieos /= 19)) &
       call fatal(label,'only ieos=2 makes sense if storing thermal energy')
#endif
    if (irealvisc < 0 .or. irealvisc > 12)  call fatal(label,'invalid setting for physical viscosity')
    if (shearparam < 0.)                     call fatal(label,'stupid value for shear parameter (< 0)')
    if (irealvisc==2 .and. shearparam > 1) call error(label,'alpha > 1 for shakura-sunyaev viscosity')
    if (iverbose > 99 .or. iverbose < -9)   call fatal(label,'invalid verboseness setting (two digits only)')
    if (icooling > 0 .and. ieos /= 2) call fatal(label,'cooling requires adiabatic eos (ieos=2)')
    if (icooling > 0 .and. (ipdv_heating <= 0 .or. ishock_heating <= 0)) &
         call fatal(label,'cooling requires shock and work contributions')
#ifdef WIND
    if (isink_radiation == 1 .and. idust_opacity == 0 .and. alpha_rad < 1.d-10) &
         call fatal(label,'no radiation pressure force! adapt isink_radiation/idust_opacity/alpha_rad')
#endif
 endif
 return

999 continue
 if (id == master) then
    call error('Phantom','input file '//trim(infile)//' not found')
    if (adjustl(infile(1:1)) /= ' ') then
       write(*,*) ' creating one with default options...'
       infilenew = trim(infile)
       call write_infile(infilenew,logfile,evfile,dumpfile,iwritein,iprint)
    endif
 endif
 call die

end subroutine read_infile

end module readwrite_infile
