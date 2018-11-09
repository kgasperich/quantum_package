program read_integrals
  PROVIDE ezfio_filename
  call run_save_mo_bielec_to_disk
end

subroutine run_save_mo_bielec_to_disk
  if (read_mo_integrals) then
    print*,'integrals_bielec/disk_access_mo_integrals == Read'
    print*,'MO bielec integrals already stored on disk'
  else
    call ezfio_set_integrals_bielec_disk_access_mo_integrals('Write')
    TOUCH read_mo_integrals read_ao_integrals write_mo_integrals write_ao_integrals
    PROVIDE mo_bielec_integrals_in_map
    call ezfio_set_integrals_bielec_disk_access_mo_integrals('Read')
    TOUCH read_mo_integrals read_ao_integrals write_mo_integrals write_ao_integrals
  endif
end
