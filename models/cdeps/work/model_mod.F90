! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! model_mod module to querry metadata from ESMF_field object

module model_mod


use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength

use time_manager_mod, only : time_type, set_time, set_date, set_calendar_type

use     location_mod, only : location_type, get_close_type, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, &
                             get_location, query_location, VERTISLEVEL, &
                             VERTISHEIGHT, set_vertical

use    utilities_mod, only : error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, &
                             to_upper

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 nc_open_file_readonly, nc_close_file, &
                                 nc_get_variable, nc_get_variable_size, &
                                 NF90_MAX_NAME

use        quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
                                   set_quad_coords, quad_lon_lat_locate, &
                                   quad_lon_lat_evaluate, quad_interp_handle, &
                                   GRID_QUAD_FULLY_IRREGULAR, &
                                   QUAD_LOCATED_CELL_CENTERS

use state_structure_mod, only : add_domain, get_domain_size, &
                                get_model_variable_indices, &
                                get_kind_string, get_varid_from_kind, &
                                get_dart_vector_index

use distributed_state_mod, only : get_state, get_state_array

use obs_kind_mod, only : get_index_for_quantity, QTY_U_CURRENT_COMPONENT, &
                         QTY_V_CURRENT_COMPONENT, QTY_LAYER_THICKNESS, &
                         QTY_DRY_LAND, QTY_SALINITY

use ensemble_manager_mod, only : ensemble_type

use ESMF  !  HK might be best to have use ESMF, use NUOPC while developing

use ESMF             , only : ESMF_VM, ESMF_VMBroadcast, ESMF_FIELDGET
use ESMF             , only : ESMF_Mesh, ESMF_GridComp, ESMF_SUCCESS, ESMF_Grid, ESMF_GridGet
use ESMF             , only : ESMF_DistGrid, ESMF_DistGridConnection, ESMF_GridGetCoord
use ESMF             , only : ESMF_GridCompSetEntryPoint, ESMF_METHOD_INITIALIZE, ESMF_DistGridCreate, ESMF_GridCompSet
use ESMF             , only : ESMF_MethodRemove, ESMF_State, ESMF_Clock, ESMF_TimeInterval
use ESMF             , only : ESMF_Field, ESMF_LOGMSG_INFO, ESMF_ClockGet
use ESMF             , only : ESMF_Time, ESMF_Alarm
! use ESMF             , only : operator(+), ESMF_TimeIntervalGet
use ESMF             , only : ESMF_TimeIntervalGet, ESMF_ClockGetAlarm, ESMF_TimeGet
use ESMF             , only : ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_StateGet
use ESMF             , only : ESMF_FieldGet, ESMF_MAXSTR, ESMF_VMBroadcast, ESMF_Array, ESMF_FieldWriteVTK
use ESMF             , only : ESMF_TraceRegionEnter, ESMF_TraceRegionExit, ESMF_GridCompGet
use ESMF             , only : ESMF_KIND_R8, ESMF_LogFoundError
use ESMF             , only : ESMF_LOGERR_PASSTHRU, ESMF_LOGWRITE


use NUOPC            , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_CompCheckSetClock
use NUOPC            , only : NUOPC_CompFilterPhaseMap
use NUOPC, only : NUOPC_Write ! HK what is the write for?
use NUOPC, only : NUOPC_SetAttribute, NUOPC_CompAttributeSet, NUOPC_Realize, NUOPC_GetAttribute

!HK three of these rename on import did not seem to work, using the model_ version for now
use NUOPC_Model, only : NUOPC_ModelGet, setVM

use cop_comp_shr, only: ChkErr
use cop_comp_shr, only: StringListGetName
use cop_comp_shr, only: StringListGetNum
use cop_comp_shr, only: StateWrite


! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, adv_1step

implicit none
private

! routines required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
         !  nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          adv_1step,              &
          init_time,              &
          init_conditions,        &
          shortest_time_between_assimilations, &
          write_model_time




character(len=256), parameter :: modName   = "model_mod"
! character(len=*), parameter :: u_FILE_u = &
      !   __FILE__
logical :: module_initialized = .false.
type(time_type) :: assimilation_time_step
integer :: dom_id ! used to access the state structure
! Suman: What does state structure means?
integer(i8) :: model_size
integer :: nfields ! number of fields in the state vector
! Grid parameters, nz is number of layers
integer :: nx=-1, ny=-1, nz=-1     ! grid counts for each field
! Suman: Explain me this!
real(r8), allocatable    :: geolon(:,:),   geolat(:,:),     & ! T
                            geolon_u(:,:), geolat_u(:,:),   & ! U
                            geolon_v(:,:), geolat_v(:,:)      ! V
type(quad_interp_handle) :: interp_t_grid,  &
                            interp_u_grid,  &
                            interp_v_grid

! Ocean vs land
real(r8), allocatable :: wet(:,:), basin_depth(:,:)

! DART state vector contents are specified in the input.nml:&model_nml namelist.
! Suman: This is something I don't have to change it, but what about model_nml?
integer, parameter :: MAX_STATE_VARIABLES = 10
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 3

! model_interpolate failure codes
integer, parameter :: NOT_IN_STATE = 12
integer, parameter :: THICKNESS_NOT_IN_STATE = 13
integer, parameter :: QUAD_LOCATE_FAILED = 14
integer, parameter :: THICKNESS_QUAD_EVALUTATE_FAILED = 15
integer, parameter :: QUAD_EVALUTATE_FAILED = 16
integer, parameter :: QUAD_ON_LAND = 17
integer, parameter :: QUAD_ON_BASIN_EDGE = 18
integer, parameter :: OBS_ABOVE_SURFACE = 20
integer, parameter :: OBS_TOO_DEEP = 22
! Suman: What is Quad?


! namelist
character(len=256) :: template_file = 'mom6.r.nc'
character(len=256) :: static_file = 'c.e22.GMOM.T62_g16.nuopc.001.mom6.static.nc'
! Suman: We need to change this file name?
character(len=256) :: ocean_geometry = 'ocean_geometry.nc'
integer  :: assimilation_period_days      = -1
integer  :: assimilation_period_seconds   = -1
character(len=vtablenamelength) :: model_state_variables(MAX_STATE_VARIABLES * NUM_STATE_TABLE_COLUMNS ) = ' '

namelist /model_nml/ template_file, static_file, ocean_geometry, assimilation_period_days, &
                     assimilation_period_seconds, model_state_variables

! Suman: namelist file 'model_nml' contain those variables as listed above, and these can be read
! from or written to using the namelist in an external file. This allows for flexible runtime configuration.
! This approach is useful for runtime configuration. Instead of recompiling the code every time a variable 
! needs to be changed, such as file paths or assimilation periods, users can modify a namelist file.
! It is used for setting up input parameters for Fortran programs. This method allows users to modify
! input variables without directly editing the source code.

interface on_land
   module procedure on_land_point
   module procedure on_land_quad
end interface on_land

! 
! Suman: We basically need to make changes where model_mod is calling functions into this namelist nedCDF file


contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.


 


subroutine static_init_model(dgcomp, rc)

type(ESMF_GridComp)  :: dgcomp
integer, intent(out) :: rc
integer  :: iunit, io
character(len=vtablenamelength) :: variable_table(MAX_STATE_VARIABLES, NUM_STATE_TABLE_COLUMNS)

integer :: state_qty_list(MAX_STATE_VARIABLES)
logical :: update_var_list(MAX_STATE_VARIABLES)

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX = 1
integer, parameter :: VAR_QTY_INDEX = 2
integer, parameter :: VAR_UPDATE_INDEX = 3
character(len=*),parameter :: subname='(static_init_model)'


! Suman: we can create one namelist for the cdeps model

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('gregorian')

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime 
! feel free to hardcode it and remove from the namelist.
assimilation_time_step = set_time(assimilation_period_seconds, &
                                  assimilation_period_days)

! verify that the model_state_variables namelist was filled in correctly.
! returns variable_table which has variable names, kinds and update strings.
call verify_state_variables(model_state_variables, nfields, variable_table, state_qty_list, update_var_list)

! Define which variables are in the model state
dom_id = add_domain(template_file, nfields, &
                    var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                    kind_list = state_qty_list(1:nfields), &
                    update_list = update_var_list(1:nfields))

model_size = get_domain_size(dom_id)  !Suman: Why didn't you use get_model_size() fuction here? 

call read_horizontal_esmf_grid(dgcomp, rc)
call setup_interpolation()

call read_num_layers(dgcomp, rc)
call read_ocean_geometry(dgcomp, rc) ! ocean vs. land and basin depth

end subroutine static_init_model



!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
function get_model_size(dgcomp, rc)

integer(i8) :: get_model_size
type(ESMF_GridComp)  :: dgcomp
integer, intent(out) :: rc

rc = ESMF_SUCCESS

if ( .not. module_initialized ) call static_init_model(dgcomp, rc)

get_model_size = get_domain_size(dom_id)

end function get_model_size

!-------------------------------------------------------------------
! Returns if the coordinate is on land or ocean






!------------------------------------------------------------------
! Given a state handle, a location, and a state quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case a positive istatus should be returned.

subroutine model_interpolate(dgcomp, rc, state_handle, ens_size, location, qty, expected_obs, istatus)

type(ESMF_GridComp)  :: dgcomp
integer, intent(out) :: rc
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer  :: which_vert, four_ilons(4), four_ilats(4), lev(ens_size,2)
integer  :: locate_status, quad_status
real(r8) :: lev_fract(ens_size)
real(r8) :: lon_lat_vert(3)
real(r8) :: quad_vals(4, ens_size)
real(r8) :: expected(ens_size, 2) ! level below and above obs
type(quad_interp_handle) :: interp
integer :: varid, i, e, thick_id
integer(i8) :: th_indx, indx(ens_size)
real(r8) :: depth_at_x(ens_size), thick_at_x(ens_size) ! depth, layer thickness at obs lat lon
logical :: found(ens_size)

rc = ESMF_SUCCESS

if ( .not. module_initialized ) call static_init_model(dgcomp, rc)

expected_obs(:) = MISSING_R8
istatus(:) = 1

varid = get_varid_from_kind(dom_id, qty)
if (varid < 0) then ! not in state
   istatus = NOT_IN_STATE
   return
endif

thick_id = get_varid_from_kind(dom_id, QTY_LAYER_THICKNESS)
if (thick_id < 0) then ! thickness not in state
   istatus = THICKNESS_NOT_IN_STATE
   return ! HK else use pseudo depth?
endif

! find which grid the qty is on
! Suman: Does it want which processor the grid is on?
interp = get_interp_handle(qty)

! unpack the location type into lon, lat, vert, vert_type
lon_lat_vert = get_location(location)
which_vert   = nint(query_location(location))

! get the indices for the 4 corners of the quad in the horizontal
call quad_lon_lat_locate(interp, lon_lat_vert(1), lon_lat_vert(2), &
                         four_ilons, four_ilats, locate_status)
if (locate_status /= 0) then
  istatus(:) = QUAD_LOCATE_FAILED
  return
endif

! check if all four corners are in the ocean
! Suman: so in order to locate these coordinates we just have to provide the corner coordinates?
if (on_land(four_ilons, four_ilats)) then
   istatus(:) = QUAD_ON_LAND
   return
endif

! find which layer the observation is in. Layer thickness is a state variable.
! HK @todo Do you need to use t_grid interp for thickess four_ilons, four_ilats?
found(:) = .false.
depth_at_x(:) = 0
FIND_LAYER: do i = 2, nz

   ! corner1
   th_indx = get_dart_vector_index(four_ilons(1), four_ilats(1), i, dom_id, thick_id)
   quad_vals(1, :) = get_state(th_indx, state_handle)
   
   ! corner2
   th_indx = get_dart_vector_index(four_ilons(1), four_ilats(2), i, dom_id, thick_id)
   quad_vals(2, :) = get_state(th_indx, state_handle)
   
   ! corner3
   th_indx = get_dart_vector_index(four_ilons(2), four_ilats(1), i, dom_id, thick_id)
   quad_vals(3, :) = get_state(th_indx, state_handle)
   
   ! corner4
   th_indx = get_dart_vector_index(four_ilons(2), four_ilats(2), i, dom_id, thick_id)
   quad_vals(4, :) = get_state(th_indx, state_handle)
   
   call quad_lon_lat_evaluate(interp, &
                              lon_lat_vert(1), lon_lat_vert(2), & ! lon, lat of obs
                              four_ilons, four_ilats, &
                              ens_size, &
                              quad_vals, & ! 4 corners x ens_size
                              thick_at_x, &
                              quad_status)
   if (quad_status /= 0) then
      istatus(:) = THICKNESS_QUAD_EVALUTATE_FAILED
      return
   endif

   depth_at_x = depth_at_x + thick_at_x

   do e = 1, ens_size
      if (lon_lat_vert(3) < depth_at_x(e)) then
         lev(e,1) = i ! layer_below
         lev(e,2) = i-1 ! layer_above
         lev_fract(e) = (depth_at_x(e) - lon_lat_vert(3)) / thick_at_x(e)
         found(e) = .true.
         if (all(found)) exit FIND_LAYER
      endif
   enddo

enddo FIND_LAYER

if (any(found .eqv. .false.)) then
   istatus(:) = OBS_TOO_DEEP
   return
endif

if (on_basin_edge(four_ilons, four_ilats, ens_size, depth_at_x)) then
   istatus(:) = QUAD_ON_BASIN_EDGE
   return
endif

do i = 1, 2 
   !HK which corner of the quad is which?
   ! corner1
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(1), four_ilats(1), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(1, :), indx, state_handle)

   ! corner2
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(1), four_ilats(2), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(2, :), indx, state_handle)

   ! corner3
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(2), four_ilats(1), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(3, :), indx, state_handle)

   ! corner4
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(2), four_ilats(2), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(4, :), indx, state_handle)

   call quad_lon_lat_evaluate(interp, &
                              lon_lat_vert(1), lon_lat_vert(2), & ! lon, lat of obs
                              four_ilons, four_ilats, &
                              ens_size, &
                              quad_vals, & ! 4 corners x ens_size
                              expected(:,i), &
                              quad_status)
   if (quad_status /= 0) then
      istatus(:) = QUAD_EVALUTATE_FAILED
      return
   else
      istatus = 0
   endif

enddo

! Interpolate between levels
! expected_obs = bot_val + lev_fract * (top_val - bot_val)
expected_obs = expected(:,1) + lev_fract(:) * (expected(:,2) - expected(:,1))

if (qty == QTY_SALINITY) then  ! convert from PSU (model) to MSU (obersvation)
   expected_obs = expected_obs/1000.0_r8
endif


end subroutine model_interpolate



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations(dgcomp, rc)

type(time_type)     :: shortest_time_between_assimilations
type(ESMF_GridComp) :: dgcomp
integer, intent(out) :: rc

rc = ESMF_SUCCESS

if ( .not. module_initialized ) call static_init_model(dgcomp, rc)

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations

! Suman: Why do we want to advance the model in time? I thought it remain stop when assimilation is happening!

!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(dgcomp, rc, index_in, location, qty)

type(ESMF_GridComp)                        :: dgcomp
integer, intent(out)                       :: rc
integer(i8),         intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: qty

real(r8) :: lat, lon
integer :: lon_index, lat_index, level, local_qty

rc = ESMF_SUCCESS

if ( .not. module_initialized ) call static_init_model(dgcomp, rc)

call get_model_variable_indices(index_in, lon_index, lat_index, level, kind_index=local_qty)

call get_lon_lat(lon_index, lat_index, local_qty, lon, lat)

location = set_location(lon, lat, real(level,r8), VERTISLEVEL)

if (present(qty)) then
   qty = local_qty
   if (on_land(lon_index, lat_index)) qty = QTY_DRY_LAND
endif


end subroutine get_state_meta_data

!------------------------------------------------------------------
subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(num) !locations
integer,             intent(in)    :: loc_qtys(num) !qty at location
integer(i8),         intent(in)    :: loc_indx(num) !state index
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

integer :: i,j,k
integer :: ii, layer ! loop variables
integer :: thick_id
integer(i8) :: indx
real(r8) :: depth(1)

! assert(which_vert == VERTISHEIGHT)

thick_id = get_varid_from_kind(dom_id, QTY_LAYER_THICKNESS)
if (thick_id < 0) then
   istatus = THICKNESS_NOT_IN_STATE
   return
endif

do ii = 1, num

   if (loc_qtys(ii) == QTY_DRY_LAND) then
      call set_vertical(locs(ii), 0.0_r8, VERTISHEIGHT)
   else

      call get_model_variable_indices(loc_indx(ii), i, j, k)

      depth = 0.0_r8
      do layer = 1, k
         indx = get_dart_vector_index(i, j, layer, dom_id, thick_id)
         depth = depth + get_state(indx, state_handle)
      enddo

      call set_vertical(locs(ii), depth(1), VERTISHEIGHT)

   endif

enddo

istatus = 0

end subroutine convert_vertical_state
!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     ! location of interest
integer,                       intent(in)    :: base_type    ! observation TYPE
type(location_type),           intent(inout) :: locs(:)      ! state locations
integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
integer,                       intent(out)   :: num_close    ! how many are close
integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'

integer :: ii ! loop index
integer :: i, j, k
real(r8) :: lon_lat_vert(3)

call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)

if (.not. present(dist)) return

! Put any land or sea floor points very far away
! so they are not updated by assimilation
do ii = 1, num_close

  if(loc_qtys(close_ind(ii)) == QTY_DRY_LAND) dist = 1.0e9_r8

  lon_lat_vert = get_location(locs(close_ind(ii))) ! assuming VERTISHEIGHT
  call get_model_variable_indices(loc_indx(ii), i, j, k)
  if ( below_sea_floor(i,j,lon_lat_vert(3)) ) dist = 1.0e9_r8

enddo

end subroutine get_close_state


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()


end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

! subroutine nc_write_model_atts(ncid, domain_id)

! integer, intent(in) :: ncid      ! netCDF file identifier
! integer, intent(in) :: domain_id

! if ( .not. module_initialized ) call static_init_model(dgcomp, rc)

! ! put file into define mode.

! call nc_begin_define_mode(ncid)

! call nc_add_global_creation_time(ncid)

! call nc_add_global_attribute(ncid, "model_source", source )
! call nc_add_global_attribute(ncid, "model", "MOM6")

! call nc_end_define_mode(ncid)

! ! Flush the buffer and leave netCDF file open
! call nc_synchronize_file(ncid)

! end subroutine nc_write_model_atts

!------------------------------------------------------------
! Read lon, lat for T,U,V grids from mom6 static file
subroutine read_horizontal_esmf_grid(dgcomp, rc)


character(len=*), parameter :: routine = 'read_horizontal_esmf_grid'

! input/output variable
type(ESMF_GridComp)         :: dgcomp
integer, intent(out)        :: rc

! local variables
integer                             :: itemCount
type(ESMF_State)                    :: importState, exportState
type(ESMF_Field)                    :: field
type(ESMF_Grid)                     :: grid
type(ESMF_Clock)                    :: clock
real(ESMF_KIND_R8), pointer         :: mask2D(:,:)
CHARACTER(ESMF_MAXSTR), allocatable :: itemNameList(:)
character(len=*),parameter          :: subname='(static_init_model)'

type(ESMF_Field)                    :: field_lon, field_lat
real(ESMF_KIND_R8), pointer         :: geolon(:,:), geolat(:,:)
real(ESMF_KIND_R8), pointer         :: geolon_u(:,:), geolat_u(:,:)
real(ESMF_KIND_R8), pointer         :: geolon_v(:,:), geolat_v(:,:)
real(ESMF_KIND_R8), pointer         :: sst_field(:,:)
integer                             :: i, j, localDE, tileCount, dimCount
integer                             :: compLBnd(2), compUBnd(2), exclLBnd(2), exclUBnd(2)


rc = ESMF_SUCCESS

! Since model_mod interacts with the DART grid component, the first step would be to get the grid component
! and then extract the export state from it.

call NUOPC_ModelGet(dgcomp, exportState=exportState, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return ! bail out

if (.not. allocated(itemNameList)) allocate(itemNameList)

! get number of fields from the state
call ESMF_StateGet(exportState, itemCount=itemCount, itemNameList=itemNameList, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return ! bail out

! loop over the field and find the sea surface temperature
do i = 1, itemCount
   if (trim(itemNameList(i))=="So_t") then
   call ESMF_LogWrite("Getting the field: "//trim(itemNameList(i)), rc=rc)

   call ESMF_StateGet(exportState, itemName=trim(itemNameList(i)), field=field, rc=rc)
   ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
   !    line=__LINE__, &
   !    file=__FILE__)) &
   !    return ! bail out

   endif
end do


! retrieve the Fortran data pointer from the Field and bounds
! Suman: I have doubt in computationalbound and exclusive bound.
call ESMF_FieldGet(field=field, localDe=0, farrayPtr=sst_field, &
   computationalLBound=compLBnd, computationalUBound=compUBnd, &
   exclusiveLBound=exclLBnd, exclusiveUBound=exclUBnd, &
   rc=rc)

! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return ! bail out


! Accessing grid information from the field
call ESMF_GridGet(grid, dimCount=dimCount, tileCount=tileCount, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return ! bail out



!-----------------------------------------------------------------------------------------
! Accessing the coordinate array for the physical position and put into an Array
call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, coordDim=1, &
   localDE=0, array=geolon, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!    line=__LINE__, &
!    file=__FILE__)) &
!    return ! bail out

call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, coordDim=2, &
   localDE=0, array=geolat, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! line=__LINE__, &
! file=__FILE__)) &
! return ! bail out

call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE1, coordDim=1, &
   localDE=0, array=geolon_u, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! line=__LINE__, &
! file=__FILE__)) &
! return ! bail out

call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE1, coordDim=2, &
   localDE=0, array=geolat_u, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! line=__LINE__, &
! file=__FILE__)) &
! return ! bail out

call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE2, coordDim=1, &
   localDE=0, array=geolon_v, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! line=__LINE__, &
! file=__FILE__)) &
! return ! bail out

call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE2, coordDim=2, &
   localDE=0, array=geolat_v, rc=rc)
! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! line=__LINE__, &
! file=__FILE__)) &
! return ! bail out

!----------------------------------------------------------------------------------


! Step 2: Adjust the longitude grid if the value exceeds 360
where(geolon > 360.0_r8)   geolon   = geolon   - 360.0_r8
where(geolon_u > 360.0_r8) geolon_u = geolon_u - 360.0_r8
where(geolon_v > 360.0_r8) geolon_v = geolon_v - 360.0_r8

! No need to close any files since we are working with ESMF_Field data
call ESMF_LogWrite("read_horizontal_grid_esmf: Done", ESMF_LOGMSG_INFO)


end subroutine read_horizontal_esmf_grid





!------------------------------------------------------------
! Read number of vertical layers from mom6 template file
subroutine read_num_layers(dgcomp, rc)

type(ESMF_GridComp)                 :: dgcomp
integer, INTENT(OUT)                :: rc
type(ESMF_State)                    :: importState, exportState
type(ESMF_Field)                    :: field
type(ESMF_Grid)                     :: grid
type(ESMF_Clock)                    :: clock
real(ESMF_KIND_R8), pointer         :: mask2D(:,:)
CHARACTER(ESMF_MAXSTR), allocatable :: itemNameList(:)
character(len=*),parameter          :: subname='(static_init_model)'
real(ESMF_KIND_R8), pointer         :: verticalCoord(:)
integer                             :: i, j, localDE, tileCount, verticalDim
integer                             :: compLBnd(3), compUBnd(3), exclLBnd(3), exclUBnd(3)

character(len=*), parameter :: routine = 'read_num_layers'
rc = ESMF_SUCCESS




call NUOPC_ModelGet(dgcomp, exportState=exportState, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
line=__LINE__, &
file=__FILE__)) &
return ! bail out


call ESMF_StateGet(exportState, itemName="So_t", field=field, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
line=__LINE__, &
file=__FILE__)) &
return ! bail out

! retrieve the Fortran data pointer from the Field and bounds
call ESMF_FieldGet(field=field, localDe=0, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
line=__LINE__, &
file=__FILE__)) &
return ! bail out

! Step 2: Get the bounds of the vertical dimension (assuming 3D grid), if it's 2d give an error message
! compLBnd and compUBnd are arrays of bounds for each dimension
call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
    computationalLBound=compLBnd, computationalUBound=compUBnd, rc=rc)
if (rc /= ESMF_SUCCESS) then
    print*, "Error getting grid bounds"
    return
endif

! Step 3: Get the vertical dimension (usually the third dimension)
verticalDim = compUBnd(3) - compLBnd(3) + 1
print*, "Number of vertical layers: ", verticalDim

! Step 4: Retrieve the vertical coordinates (e.g., heights or depths)
! Assuming the vertical coordinates are stored at the center stagger location
call ESMF_GridGetCoord(grid, coordDim=3, staggerloc=ESMF_STAGGERLOC_CENTER, &
    farrayPtr=verticalCoord, rc=rc)
if (rc /= ESMF_SUCCESS) then
    print*, "Error getting vertical coordinates"
    return
endif

! Print the vertical coordinates
do i = 1, verticalDim
    print*, "Vertical layer ", i, " coordinate: ", verticalCoord(i)
end do

! Suman: We got the vertical coordinates, but still have to calculate total number of vertical layer
! and layer thickness!

end subroutine read_num_layers





!------------------------------------------------------------
! ocean_geom are 2D state sized static data
! HK Do these arrays become too big in high res cases?
subroutine read_ocean_geometry(dgcomp, rc)

! variable declations
type(ESMF_GridComp)         :: dgcomp
type(ESMF_State)            :: exportState
integer, INTENT(OUT)        :: rc
integer                     :: i, j, nx, ny
type(ESMF_Grid)             :: grid
type(ESMF_Field)            :: field
real(ESMF_KIND_R8), pointer :: basin_depth(:,:)
real(ESMF_KIND_R8), pointer :: mask2D(:,:)
character(len=*), parameter :: routine = 'read_ocean_geometry'
integer                     :: compLBnd(2), compUBnd(2)  ! Computational bounds

rc = ESMF_SUCCESS




call NUOPC_ModelGet(dgcomp, exportState=exportState, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
line=__LINE__, &
file=__FILE__)) &
return ! bail out


call ESMF_StateGet(exportState, itemName="basin_depth", field=field, rc=rc)
if (ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
line=__LINE__, &
file=__FILE__)) &
return ! bail out

! retrieve the Fortran data pointer from the Field and bounds
call ESMF_FieldGet(field=field, localDe=0, farrayPtr=basin_depth, &
                     computationalLBound=compLBnd, computationalUBound=compUBnd, rc=rc)
if (ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
line=__LINE__, &
file=__FILE__)) &
return ! bail out

allocate(mask2D(compLBnd(1):compUBnd(1), compLBnd(2):compUBnd(2)))

call ESMF_GridGetItem(grid, localDE=0, staggerloc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, farrayPtr=mask2D, rc=rc)
if (ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
line=__LINE__, &
file=__FILE__)) &
return ! bail out

! ! Now you can access mask2D as a regular Fortran array
! do i = 1, size(mask2D, 1)
!    do j = 1, size(mask2D, 2)

!       if (mask2D(i,j) ==0) then
!          Print *, "Point (", i, j, ") is masked (inactive)"
!       endif
!    enddo
! enddo

!nx and ny are the dimensions of the data
nx = compUBnd(1) - compLBnd(1) + 1
ny = compUBnd(2) - compLBnd(2) + 1

end subroutine read_ocean_geometry


 
!------------------------------------------------------------
! wet is a 2D array of ones and zeros
! 1 is ocean
! 0 is land

! Suman: functions to corsscheck if the point falls on land or ocean

function on_land_quad(ilon, ilat)

integer :: ilon(4), ilat(4) ! these are indices into lon, lat
logical ::  on_land_quad

real(r8), pointer :: mask2D(:,:) ! Add declaration for mask2D

if ( mask2D(ilon(1), ilat(1)) + &
     mask2D(ilon(1), ilat(2)) + &
     mask2D(ilon(2), ilat(1)) + &
     mask2D(ilon(2), ilat(2))  < 4) then
   on_land_quad = .true.
else
   on_land_quad = .false.
endif

end function on_land_quad

!------------------------------------------------------------
function on_land_point(ilon, ilat)

integer :: ilon, ilat ! these are indices into lon, lat
logical :: on_land_point

real(r8), pointer :: mask2D(:,:) ! Add declaration for mask2D

if ( mask2D(ilon, ilat) + &
     mask2D(ilon, ilat) + &
     mask2D(ilon, ilat) + &
     mask2D(ilon, ilat)  < 4) then
   on_land_point = .true.
else
   on_land_point = .false.
endif


end function on_land_point

!------------------------------------------------------------
! basin_depth is a 2D array with the basin depth
function on_basin_edge(ilon, ilat, ens_size, depth)

! indices into lon, lat lev
integer, intent(in)  :: ilon(4), ilat(4)
integer, intent(in)  :: ens_size
real(r8), intent(in) :: depth(ens_size)
logical :: on_basin_edge

integer  :: i, e
real(r8) :: d(4) ! basin depth at each corner

d(1) = basin_depth(ilon(1), ilat(1))
d(2) = basin_depth(ilon(1), ilat(2))
d(3) = basin_depth(ilon(2), ilat(1))
d(4) = basin_depth(ilon(2), ilat(2))

do e = 1, ens_size
   do i = 1, 4
      if (d(i) < depth(e)) then
         on_basin_edge = .true.
         return
      endif
   enddo
enddo

! four points are in the ocean
on_basin_edge = .false.

end function on_basin_edge

!------------------------------------------------------------
function below_sea_floor(ilon, ilat, depth)

! indices into lon, lat lev
integer,  intent(in) :: ilon, ilat
real(r8), intent(in) :: depth
logical :: below_sea_floor

if (basin_depth(ilon, ilat) < depth) then
   below_sea_floor = .true.
else
   below_sea_floor = .false.
endif

end function below_sea_floor

!------------------------------------------------------------
! longitude and latitide values from indices
subroutine get_lon_lat(lon_indx, lat_indx, qty, lon, lat)

integer, intent(in) :: lon_indx, lat_indx, qty
real(r8) :: lon, lat

if (on_u_grid(qty)) then
   lon = geolon_u(lon_indx, lat_indx)
   lat = geolat_u(lon_indx, lat_indx)
elseif (on_v_grid(qty)) then
   lon = geolon_v(lon_indx, lat_indx)
   lat = geolat_v(lon_indx, lat_indx)
else ! T grid
   lon = geolon(lon_indx, lat_indx)
   lat = geolat(lon_indx, lat_indx)
endif

end subroutine get_lon_lat

!------------------------------------------------------------
function on_v_grid(qty)

integer, intent(in)  :: qty
logical :: on_v_grid

if (qty == QTY_V_CURRENT_COMPONENT) then
  on_v_grid = .true.
else
  on_v_grid = .false.
endif

end function on_v_grid

!----------------------------------------------------------
function on_u_grid(qty)

integer, intent(in)  :: qty
logical :: on_u_grid

if (qty == QTY_U_CURRENT_COMPONENT) then
  on_u_grid = .true.
else
  on_u_grid = .false.
endif

end function on_u_grid

!----------------------------------------------------------
function on_t_grid(qty)

integer, intent(in)  :: qty
logical :: on_t_grid

if (qty == QTY_U_CURRENT_COMPONENT .or. qty == QTY_V_CURRENT_COMPONENT) then
  on_t_grid = .false.
else
  on_t_grid = .true.
endif

end function on_t_grid

!------------------------------------------------------------
function on_layer(qty)

logical :: on_layer
integer :: qty

! Salt, Temp, u, v all on layer 
on_layer = .true.

end function on_layer

!------------------------------------------------------------
subroutine setup_interpolation()

! T
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_t_grid)

call set_quad_coords(interp_t_grid, geolon, geolat)

! U
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_u_grid)

call set_quad_coords(interp_u_grid, geolon_u, geolat_u)


! V
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_v_grid)

call set_quad_coords(interp_v_grid, geolon_v, geolat_v)

end subroutine setup_interpolation

!------------------------------------------------------------
! return the appropriate quad_interp handle
function get_interp_handle(qty)

type(quad_interp_handle) :: get_interp_handle
integer, intent(in) :: qty

if (on_v_grid(qty)) then
  get_interp_handle = interp_v_grid
elseif (on_v_grid(qty)) then
  get_interp_handle = interp_u_grid
else
  get_interp_handle = interp_t_grid
endif

end function


!------------------------------------------------------------------
! Verify that the namelist was filled in correctly, and check
! that there are valid entries for the dart_kind.
! Returns a table with columns:
!
! netcdf_variable_name ; dart_qty_string ; update_string

subroutine verify_state_variables(state_variables, ngood, table, qty_list, update_var)

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: qty_list(:)   ! kind number
logical,           intent(out) :: update_var(:) ! logical update

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, update
character(len=256) :: string1, string2
type(ESMF_GridComp)  :: dgcomp
rc = ESMF_SUCCESS

if ( .not. module_initialized ) call static_init_model(dgcomp, rc)

nrows = size(table,1)

ngood = 0

MyLoop : do i = 1, nrows

   varname = trim(state_variables(3*i -2))
   dartstr = trim(state_variables(3*i -1))
   update  = trim(state_variables(3*i   ))
   
   call to_upper(update)

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(update)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' .and. table(i,3) == ' ') exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' .or. table(i,3) == ' ' ) then
      string1 = 'model_nml:model_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1)
   endif

   ! Make sure DART qty is valid

   qty_list(i) = get_index_for_quantity(dartstr)
   if( qty_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1)
   endif
   
   ! Make sure the update variable has a valid name

   select case (update)
      case ('UPDATE')
         update_var(i) = .true.
      case ('NO_COPY_BACK')
         update_var(i) = .false.
      case default
         write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
         write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update)
         call error_handler(E_ERR,'verify_state_variables',string1, text2=string2)
   end select

   ngood = ngood + 1
enddo MyLoop


end subroutine verify_state_variables

!--------------------------------------------------------------------
function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type) :: read_model_time

integer :: ncid
character(len=*), parameter :: routine = 'read_model_time'
real(r8) :: days
type(time_type) :: mom6_time
integer :: dart_base_date_in_days, dart_days

dart_base_date_in_days = 584388 ! 1601 1 1 0 0
ncid = nc_open_file_readonly(filename, routine)

call nc_get_variable(ncid, 'Time', days, routine)

call nc_close_file(ncid, routine)

! MOM6 counts days from year 1
! DART counts days from 1601 
dart_days = int(days) - dart_base_date_in_days

read_model_time = set_time(0,dart_days)

end function read_model_time


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

