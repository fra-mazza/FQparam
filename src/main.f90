!> @brief Main program for the FQ parameter optimizer.
!> @details This program drives the optimization of Fluctuating Charge (FQ) model parameters.
!>          It reads the molecular system definition and QM reference data, optionally checks
!>          the analytical gradients against numerical ones, and then runs the optimization
!>          to find the best-fit chi and eta parameters.
program optimizer
    use types_module
    use input_module
    use optimizer_module
    use debug_module
    use fq_util

    implicit none

    type(molecule_t):: molecule
    type(qm_datapoint_t), allocatable:: qm_database(:)
    type(optimization_settings_t):: opt_settings
    type(optimization_params_t):: opt_params
    !> The primary data structure holding all molecule-specific FQ data.
    type(FQMol_t), dimension(:), allocatable:: FQMolecules
    real(wp), dimension(:,:), allocatable:: FQJinv

    ! 1. Read all input data from files and hard-coded values.
    call read_input_data(molecule, qm_database, opt_settings, opt_params)

    ! 2. Create and initialize the persistent FQ molecules data structure.
    !    This is allocated once and passed through the program to avoid
    !    repeated memory allocation/deallocation cycles.
    call CreateFQMolecules(molecule%nMol, molecule%nAtoms, molecule%nAtomsInMol, &
                         molecule%MolID, molecule%MolCharges, molecule%coordinates, &
                         molecule%symbols, FQMolecules)

    ! Calculate and print initial unperturbed charges
    call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)
    call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
    call calculate_and_print_unperturbed_charges(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
    if (allocated(FQJinv)) deallocate(FQJinv)

    ! 3. (Optional) Perform gradient checks for debugging.
    !    This is highly recommended when changing the cost function or its derivatives.
    call check_gradients(opt_params, molecule, qm_database, opt_settings)

    ! 4. Run the main optimization routine.
    call optimize_parameters(opt_params, molecule, qm_database, opt_settings, FQMolecules, molecule)

    ! 5. Print the final, optimized parameters.
    print *, 'Final Chi values: ', opt_params%chi
    print *, 'Final Eta values: ', opt_params%eta

    ! Calculate and print final unperturbed charges
    call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)
    call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
    call calculate_and_print_unperturbed_charges(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
    if (allocated(FQJinv)) deallocate(FQJinv)

    ! 6. Clean up memory.
    call DestroyFQMolecules(FQMolecules)
    if (allocated(qm_database)) deallocate(qm_database)
    if (allocated(molecule%atom_types)) deallocate(molecule%atom_types)
    if (allocated(molecule%symbols)) deallocate(molecule%symbols)
    if (allocated(molecule%type_symbols)) deallocate(molecule%type_symbols)
    if (allocated(molecule%coordinates)) deallocate(molecule%coordinates)
    if (allocated(molecule%nAtomsInMol)) deallocate(molecule%nAtomsInMol)
    if (allocated(molecule%MolID)) deallocate(molecule%MolID)
    if (allocated(molecule%MolCharges)) deallocate(molecule%MolCharges)
    if (allocated(opt_params%chi)) deallocate(opt_params%chi)
    if (allocated(opt_params%eta)) deallocate(opt_params%eta)

end program optimizer
