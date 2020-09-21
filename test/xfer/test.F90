program test
  use phot_kind_mod, only: r8 => kind_phot

  use tuv_radiation_transfer, only: tuv_radiation_transfer_init
  use tuv_radiation_transfer, only: tuv_radiation_transfer_run

  use molec_ox_xsect, only: molec_ox_xsect_init
  use molec_ox_xsect, only: molec_ox_xsect_run
  use params_mod, only: input_data_root
  use wavelength_grid, only: wavelength_grid_init, nwave

  use environ_conditions_mod, only: environ_conditions_create, environ_conditions

  implicit none

  character(len=*), parameter :: wlgridfile = &
       'data/wavelength_grid.nc'

  integer :: errflg
  character(len=444) :: errmsg

  real(r8) :: zenith
  real(r8) :: albedo
  real(r8), allocatable :: alt(:)
  real(r8), allocatable :: press_mid(:)
  real(r8), allocatable :: press_int(:)
  real(r8), allocatable :: temp(:)

  real(r8), allocatable :: o2vmrcol(:)
  real(r8), allocatable :: o3vmrcol(:)
  real(r8), allocatable :: so2vmrcol(:)
  real(r8), allocatable :: no2vmrcol(:)

  real(r8), allocatable :: cldfrac(:)
  real(r8), allocatable :: cldwat(:)
  real(r8), allocatable :: srb_o2_xs(:,:)
  real(r8), allocatable :: dto2(:,:)
  real(r8), allocatable :: radfld(:,:)

  integer :: unitn
  character(len=*), parameter :: nml_file = 'xfer/test_nml'
  character(len=128) :: env_conds_file='NONE'
  real :: env_conds_lat=45.
  real :: env_conds_lon=180.
  type(environ_conditions),pointer :: colEnvConds => null()

  integer :: nlevels,k

  namelist /drv_opts/ env_conds_file, env_conds_lat, env_conds_lon

  write(*,*) 'BEGIN TEST'

  open(newunit=unitn, file=trim(nml_file), status='old')
  read(unit=unitn, nml=drv_opts)
  close(unitn)

  input_data_root = 'data'
  call wavelength_grid_init( wlgridfile, errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if
  
  colEnvConds => environ_conditions_create( env_conds_file, lat=env_conds_lat, lon=env_conds_lon )

  nlevels = colEnvConds%nlayers()

  allocate(srb_o2_xs(nwave,nlevels), dto2(nlevels,nwave))
  
  allocate(radfld(nwave,nlevels))



  allocate(alt(nlevels))
  allocate(press_mid(nlevels))
  allocate(press_int(nlevels+1))
  allocate(temp(nlevels))

  allocate(o2vmrcol(nlevels))
  allocate(o3vmrcol(nlevels))
  allocate(so2vmrcol(nlevels))
  allocate(no2vmrcol(nlevels))

  allocate(cldfrac(nlevels))
  allocate(cldwat(nlevels))

  call molec_ox_xsect_init( errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call tuv_radiation_transfer_init( r8, errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  zenith = colEnvConds%getsrf('SZA')

  albedo = colEnvConds%getsrf('ASDIR')
  press_mid(:nlevels) = colEnvConds%press_mid(nlevels)
  press_int(:nlevels+1) = colEnvConds%press_int(nlevels+1)
  alt(:nlevels) = colEnvConds%getcol('Z3',nlevels) ! meters
  temp(:nlevels) = colEnvConds%getcol('T',nlevels)

  o2vmrcol(:nlevels) = colEnvConds%getcol('O2',nlevels)
  o3vmrcol(:nlevels) = colEnvConds%getcol('O3',nlevels)
  so2vmrcol(:nlevels) = colEnvConds%getcol('SO2',nlevels)
  no2vmrcol(:nlevels) = colEnvConds%getcol('NO2',nlevels)

  cldfrac = 0._r8
  cldwat = 0._r8
  
  call molec_ox_xsect_run( nlevels, zenith, alt, temp, press_mid, press_int(1), o2vmrcol, dto2, srb_o2_xs, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call tuv_radiation_transfer_run( nlevels, nwave, zenith, albedo, press_mid, press_int(1), alt, temp, o3vmrcol, &
                                   so2vmrcol, no2vmrcol, cldfrac, cldwat, dto2, radfld, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if


  write(*,*) radfld(:,1)
  write(*,*) 'nwave=', nwave
  write(*,*) radfld(:,nlevels)
  
  write(*,*) 'END TEST'

end program test
