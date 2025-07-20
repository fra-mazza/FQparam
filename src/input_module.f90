!> @brief Module for reading input data.
!> @details This module is responsible for reading all the necessary input for the FQparam
!>          calculation. This includes the static molecular topology, the QM reference data
!>          points, initial FQ parameters, and optimization settings. Currently, most of
!>          this data is hard-coded, but it reads the core QM data from the external file
!>          'FQpar_data.dat'.
module input_module
    use types_module
    implicit none

contains

    !> @brief Reads all input data required for the optimization.
    !> @details This subroutine populates the main data structures of the program.
    !>          It defines the molecular structure (water molecule), sets the optimization
    !>          weights and initial parameters, and reads the QM reference data from an
    !>          external file. It also handles the conversion of coordinates from Angstroms
    !>          to atomic units.
    !> @param[out] molecule The molecular system definition.
    !> @param[out] qm_database The database of reference QM calculations.
    !> @param[out] opt_settings The settings for the optimization run.
    !> @param[out] opt_params The initial optimization parameters.
    subroutine read_input_data(molecule, qm_database, opt_settings, opt_params)
        type(molecule_t), intent(out):: molecule
        type(qm_datapoint_t), allocatable, intent(out):: qm_database(:)
        type(optimization_settings_t), intent(out):: opt_settings
        type(optimization_params_t), intent(out):: opt_params

        integer:: i, iostat, nQMcalc
        !> Conversion factor from Angstrom to Bohr (atomic units of length).
        real(wp), parameter:: Angstrom = 0.52917721090299996

        ! --- Hard-coded system definition (1 water molecule) ---
        molecule%nMol = 1
        molecule%nAtoms = 3
        molecule%nAtomTypes = 2  ! HW and OW

        ! --- Hard-coded optimization settings---
        opt_settings%w_eint = 10.0_wp
        opt_settings%w_dip = 2.0_wp
        opt_settings%w_quad = 0.0_wp
        opt_settings%max_iter = 10000000
        opt_settings%tolerance = 1.0e-9_wp  ! A reasonable default tolerance

        ! --- Read QM data from external file---
        open(unit = 16, file='FQpar_data.dat' , iostat = iostat, action='read')
        if (iostat /= 0) then
            print *, "Error opening FQpar_data.dat. Make sure the file exists."
            stop
        endif
        read(16, *) nQMcalc
        read(16, '(a)')  ! Skip header line

        ! --- Allocate all data structures---
        allocate(molecule%atom_types(molecule%nAtoms), &
                 molecule%nAtomsInMol(molecule%nMol), &
                 molecule%MolID(molecule%nAtoms), &
                 opt_params%chi(molecule%nAtomTypes), &
                 opt_params%eta(molecule%nAtomTypes), &
                 molecule%coordinates(3, molecule%nAtoms), &
                 molecule%symbols(molecule%nAtoms), &
                 molecule%type_symbols(molecule%nAtomTypes), &
                 molecule%MolCharges(molecule%nMol), &
                 qm_database(nQMcalc))

        ! --- Hard-coded molecular topology and initial parameters---
        molecule%atom_types = (/ 1, 1, 2 /)  ! HW, HW, OW
        molecule%nAtomsInMol = (/ 3 /)
        molecule%MolCharges = (/ 0.0_wp /)
        molecule%MolID =(/ 1, 1, 1 /)
        opt_params%chi = (/ 0.01_wp, 0.1_wp/)
        opt_params%eta = (/ 0.6_wp, 0.6_wp/)
        molecule%coordinates(:, 1) = (/    0.75308062_wp,     0.60025412_wp,     0.00000000_wp /)
        molecule%coordinates(:, 2) = (/   -0.75308062_wp,     0.60025412_wp,     0.00000000_wp /)
        molecule%coordinates(:, 3) = (/    0.00000000_wp,    -0.00025412_wp,     0.00000000_wp /)
        molecule%symbols = (/ 'HW ', 'HW ', 'OW ' /)
        molecule%type_symbols = (/ 'HW ', 'OW ' /)

        opt_settings%dipole_origin = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
        opt_settings%quad_origin = (/ 0.0_wp, 0.12652035_wp, 0.0_wp /)  ! Center of mass

        ! Read QM data points from file
        do i = 1, nQMcalc
          read(16, *) qm_database(i)%q_coords(:), qm_database(i)%q_int, qm_database(i)%ext_dipole(:), &
                      qm_database(i)%eint_qm, qm_database(i)%dip_qm(:), &
                      qm_database(i)%quad_qm(:)
        enddo

        close(16)

        ! Convert coordinates from Angstrom to atomic units (Bohr)
        molecule%coordinates = molecule%coordinates/Angstrom
        do i = 1, nQMcalc
            qm_database(i)%q_coords = qm_database(i)%q_coords/Angstrom
        enddo

    end subroutine read_input_data

end module input_module
