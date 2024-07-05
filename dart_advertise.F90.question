module OceanVarAdvertiseMod
    use esmf
    use netcdf

    implicit none

    ! Public interfaces
    public :: AdvertiseOceanVariables

contains

    subroutine AdvertiseOceanVariables(importState, rc)
        type(ESMF_State) :: importState
        integer, intent(out) :: rc

        ! Local variables
        integer :: ncid, status

        rc = ESMF_SUCCESS

        ! Open the NetCDF file
        status = openNetCDF(ncid)
        if (status /= ESMF_SUCCESS) then
            rc = ESMF_FAILURE
            return
        endif

        ! Advertise variables
        call advertiseVariables(ncid, importState, rc)
        if (rc /= ESMF_SUCCESS) return

        ! Close the NetCDF file
        call closeNetCDF(ncid)
    end subroutine AdvertiseOceanVariables

    function openNetCDF(ncid) result(status)
        integer, intent(out) :: ncid
        integer :: status

        status = nf90_open('ocean_state.nc', nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            status = ESMF_FAILURE
        else
            status = ESMF_SUCCESS
        endif
    end function openNetCDF

    subroutine advertiseVariables(ncid, importState, rc)
        integer, intent(in) :: ncid
        type(ESMF_State) :: importState
        integer, intent(out) :: rc
        integer :: varid, status
        character(len=128) :: varName

        varName = 'temperature'  ! Example variable
        status = nf90_inq_varid(ncid, varName, varid)
        if (status == nf90_noerr) then
            call ESMF_StateAdd(importState, fieldName=varName, rc=rc)
            if (ESMF_LogFoundError(rc, msg='NUOPC Error:', line=__LINE__, file=__FILE__)) return
        else
            rc = ESMF_FAILURE
        endif
    end subroutine advertiseVariables

    subroutine closeNetCDF(ncid)
        integer, intent(in) :: ncid
        integer :: status

        status = nf90_close(ncid)
        if (status /= nf90_noerr) then
            print *, 'Error closing NetCDF file.'
        endif
    end subroutine closeNetCDF

end module OceanVarAdvertiseMod
