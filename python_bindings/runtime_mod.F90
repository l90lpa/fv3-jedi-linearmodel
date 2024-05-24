module runtime_mod
    use fms_mod,  only: fms_init, fms_end
    use mpp_mod,  only: mpp_pe, mpp_npes
    use field_manager_mod,  only: field_manager_init
    implicit none
contains

    subroutine runtime_init(mpi_comm, fms_nml_file, field_table_file)

        integer         , intent(in) :: mpi_comm
        character(len=*), intent(in) :: fms_nml_file, field_table_file
        
        ! Initialize fms, mpp, etc.
        call fms_init(localcomm = mpi_comm, alt_input_nml_path = fms_nml_file)
        
        ! Initialize the tracers
        call field_manager_init(table_name = field_table_file)
        
    end subroutine runtime_init

    subroutine runtime_exit()
        call fms_end()
    end subroutine runtime_exit

    integer function runtime_pe()
        runtime_pe = mpp_pe()
    end function runtime_pe

    integer function runtime_npes()
        runtime_npes = mpp_npes()
    end function runtime_npes
    
end module runtime_mod