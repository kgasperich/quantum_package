program read_integrals
  PROVIDE ezfio_filename
  call run_save_mo_df_to_disk
end

subroutine run_save_mo_df_to_disk
  if (read_df_mo_integral_array) then
    print*,'integrals_bielec/disk_access_df_mo_integral_array == Read'
    print*,'MO df integrals already stored on disk'
  else
    call ezfio_set_integrals_bielec_disk_access_df_mo_integral_array('Write')
    PROVIDE df_mo_integral_array
    ! this should already be set by provider
!    call ezfio_set_integrals_bielec_disk_access_df_mo_integral_array('Read')
  endif
end
