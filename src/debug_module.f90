!> @brief Module for debugging and verifying analytical gradients.
!> @details This module provides subroutines to check the correctness of the analytical
!>          derivatives (d(Charge)/d(parameter) and d(Cost)/d(parameter)) by comparing
!>          them against derivatives calculated numerically using finite differences.
!>          This is a crucial step to ensure the gradient-based optimization will work correctly.
module debug_module
    use types_module
    use cost_function_module
    use fq_util
    implicit none

contains

    !> @brief Numerically checks the analytical derivatives of fluctuating charges (dFQ/dp).
    !> @details This subroutine perturbs each optimization parameter (chi and eta) by a small
    !>          amount `h`, recalculates the charges, and computes a numerical derivative using
    !>          the central difference formula. It then prints a comparison with the analytical
    !>          derivative calculated by `GetFQsDerivatives`.
    !> @param[inout] opt_params The optimization parameters.
    !> @param[in] molecule The molecular system definition.
    !> @param[in] qm_database The database of reference QM calculations.
    !> @param[in] opt_settings The settings for the optimization run.
    !> @param[inout] FQMolecules The FQ molecules data structure.
    subroutine check_dq_dp_derivatives(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        type(optimization_params_t), intent(inout):: opt_params
        type(molecule_t), intent(in):: molecule
        type(qm_datapoint_t), allocatable, intent(in):: qm_database(:)
        type(optimization_settings_t), intent(in):: opt_settings
        type(FQMol_t), dimension(:), intent(inout):: FQMolecules

        real(wp), parameter:: h = 1.0e-2
        real(wp), dimension(molecule%nAtoms):: FQs_plus, FQs_minus
        real(wp), dimension(molecule%nAtoms, molecule%nAtomTypes):: dqdchi_analytical, dqdeta_analytical
        real(wp), dimension(molecule%nAtoms, molecule%nAtomTypes):: dqdchi_numerical, dqdeta_numerical
        real(wp), dimension(molecule%nAtoms):: V
        real(wp), dimension(:,:), allocatable:: FQJinv
        integer:: i, j, k

        ! Update parameters and calculate analytical derivatives for the first QM point
        call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)
        call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
        call GetPotentialPointCharge(molecule%nMol, molecule%nAtoms, FQMolecules, &
            V, qm_database(1)%q_int, qm_database(1)%q_coords, qm_database(1)%ext_dipole)
        call UpdateFQs(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv, V)
        call GetFQsDerivatives(molecule%nMol, molecule%nAtoms, molecule%nAtomTypes, &
            molecule%atom_types, FQMolecules, FQJinv, V, dqdchi_analytical, dqdeta_analytical)

        print *, '--------------------------------------------------'
        print *, 'd(FQ)/d(param) Derivatives Check'
        print *, '--------------------------------------------------'

        ! Check dq/dchi
        print *, 'Checking d(FQ)/d(chi)'
        do j = 1, molecule%nAtomTypes  ! Loop over chi parameters
            do i = 1, molecule%nAtoms  ! Loop over FQ charges
                ! Perturb chi(j)
                opt_params%chi(j) = opt_params%chi(j) + h
                call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)
                call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
                call UpdateFQs(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv, V)
                FQs_plus = FQMolecules(1)%FQs

                opt_params%chi(j) = opt_params%chi(j) - 2.0_wp*h
                call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)
                call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
                call UpdateFQs(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv, V)
                FQs_minus = FQMolecules(1)%FQs

                ! Reset parameter
                opt_params%chi(j) = opt_params%chi(j) + h

                dqdchi_numerical(i, j) = (FQs_plus(i) - FQs_minus(i)) / (2.0_wp*h)

                print '(a, i2, a, i2, a, 4(2x, es12.4e2))', 'dFQ(', i, ')/dchi(', j, ')', &
                    dqdchi_analytical(i, j), dqdchi_numerical(i, j), &
                    abs(dqdchi_analytical(i, j) - dqdchi_numerical(i, j)), &
                    abs(dqdchi_analytical(i, j) - dqdchi_numerical(i, j)) / abs(dqdchi_analytical(i, j))
            enddo
        enddo

        ! Check dq/deta
        print *, 'Checking d(FQ)/d(eta)'
        do j = 1, molecule%nAtomTypes  ! Loop over eta parameters
            do i = 1, molecule%nAtoms  ! Loop over FQ charges
                ! Perturb eta(j)
                opt_params%eta(j) = opt_params%eta(j) + h
                call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)
                call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
                call UpdateFQs(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv, V)
                FQs_plus = FQMolecules(1)%FQs

                opt_params%eta(j) = opt_params%eta(j) - 2.0_wp*h
                call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)
                call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)
                call UpdateFQs(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv, V)
                FQs_minus = FQMolecules(1)%FQs

                ! Reset parameter
                opt_params%eta(j) = opt_params%eta(j) + h

                dqdeta_numerical(i, j) = (FQs_plus(i) - FQs_minus(i)) / (2.0_wp*h)

                print '(a, i2, a, i2, a, 4(2x, es12.4e2))', 'dFQ(', i, ')/deta(', j, ')', &
                    dqdeta_analytical(i, j), dqdeta_numerical(i, j), &
                    abs(dqdeta_analytical(i, j) - dqdeta_numerical(i, j)), &
                    abs(dqdeta_analytical(i, j) - dqdeta_numerical(i, j)) / abs(dqdeta_analytical(i, j))
            enddo
        enddo

        if (allocated(FQJinv)) deallocate(FQJinv)
        print *, '--------------------------------------------------'

    end subroutine check_dq_dp_derivatives

    !> @brief Numerically checks the analytical gradient of the interaction energy cost.
    !> @details This subroutine isolates the interaction energy component of the cost function
    !>          and compares its analytical gradient with a numerical one computed via finite differences.
    !> @param[inout] opt_params The optimization parameters.
    !> @param[in] molecule The molecular system definition.
    !> @param[in] qm_database The database of reference QM calculations.
    !> @param[inout] opt_settings The settings for the optimization run (weights are modified).
    !> @param[inout] FQMolecules The FQ molecules data structure.
    subroutine check_eint_gradient(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        type(optimization_params_t), intent(inout):: opt_params
        type(molecule_t), intent(in):: molecule
        type(qm_datapoint_t), allocatable, intent(in):: qm_database(:)
        type(optimization_settings_t), intent(inout):: opt_settings
        type(FQMol_t), dimension(:), intent(inout):: FQMolecules

        real(wp), parameter:: h = 1.0e-6
        real(wp):: cost_plus, cost_minus, num_grad
        real(wp):: cost
        real(wp), dimension(2*molecule%nAtomTypes):: analytical_grad, dummy_grad
        integer:: i

        ! Save original weights
        real(wp):: original_w_eint, original_w_dip, original_w_quad
        original_w_eint = opt_settings%w_eint
        original_w_dip = opt_settings%w_dip
        original_w_quad = opt_settings%w_quad

        ! Set weights for Eint only
        opt_settings%w_eint = 1.0_wp
        opt_settings%w_dip = 0.0_wp
        opt_settings%w_quad = 0.0_wp

        ! First, get the analytical gradient at the current point
        call calculate_cost_and_gradient(opt_params, molecule, qm_database, opt_settings, cost, analytical_grad, &
                                        FQMolecules, qm_data_idx = 1)

        print *, '--------------------------------------------------'
        print *, 'Eint Gradient Check'
        print *, '--------------------------------------------------'
        print *, 'Parameter', '      Analytical', '      Numerical', '      Abs. Diff.', '   Rel. Diff.'

        ! Check gradients for chi parameters
        do i = 1, molecule%nAtomTypes
            ! Calculate cost at p+h
            opt_params%chi(i) = opt_params%chi(i) + h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_plus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Calculate cost at p-h
            opt_params%chi(i) = opt_params%chi(i) - 2.0_wp*h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_minus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Reset parameter
            opt_params%chi(i) = opt_params%chi(i) + h

            ! Calculate numerical gradient
            num_grad = (cost_plus-cost_minus) / (2.0_wp*h)

            print '(a, i2, a, 4(2x, es12.4e2))', 'chi(', i, ')', &
                analytical_grad(i), num_grad, &
                abs(analytical_grad(i) - num_grad), &
                abs(analytical_grad(i) - num_grad) / abs(analytical_grad(i))
        enddo

        ! Check gradients for eta parameters
        do i = 1, molecule%nAtomTypes
            ! Calculate cost at p+h
            opt_params%eta(i) = opt_params%eta(i) + h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_plus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Calculate cost at p-h
            opt_params%eta(i) = opt_params%eta(i) - 2.0_wp*h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_minus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Reset parameter
            opt_params%eta(i) = opt_params%eta(i) + h

            ! Calculate numerical gradient
            num_grad = (cost_plus-cost_minus) / (2.0_wp*h)

            print '(a, i2, a, 4(2x, es12.4e2))', 'eta(', i, ')', &
                analytical_grad(molecule%nAtomTypes+i), num_grad, &
                abs(analytical_grad(molecule%nAtomTypes+i) - num_grad), &
                abs(analytical_grad(molecule%nAtomTypes+i) - num_grad) / abs(analytical_grad(molecule%nAtomTypes+i))
        enddo

        ! Restore original weights
        opt_settings%w_eint = original_w_eint
        opt_settings%w_dip = original_w_dip
        opt_settings%w_quad = original_w_quad

        print *, '--------------------------------------------------'

    end subroutine check_eint_gradient

    !> @brief Numerically checks the analytical gradient of the dipole moment cost.
    !> @details This subroutine isolates the dipole moment component of the cost function
    !>          and compares its analytical gradient with a numerical one computed via finite differences.
    !> @param[inout] opt_params The optimization parameters.
    !> @param[in] molecule The molecular system definition.
    !> @param[in] qm_database The database of reference QM calculations.
    !> @param[inout] opt_settings The settings for the optimization run (weights are modified).
    !> @param[inout] FQMolecules The FQ molecules data structure.
    subroutine check_dip_gradient(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        type(optimization_params_t), intent(inout):: opt_params
        type(molecule_t), intent(in):: molecule
        type(qm_datapoint_t), allocatable, intent(in):: qm_database(:)
        type(optimization_settings_t), intent(inout):: opt_settings
        type(FQMol_t), dimension(:), intent(inout):: FQMolecules

        real(wp), parameter:: h = 1.0e-6
        real(wp):: cost_plus, cost_minus, num_grad
        real(wp):: cost
        real(wp), dimension(2*molecule%nAtomTypes):: analytical_grad, dummy_grad
        integer:: i

        ! Save original weights
        real(wp):: original_w_eint, original_w_dip, original_w_quad
        original_w_eint = opt_settings%w_eint
        original_w_dip = opt_settings%w_dip
        original_w_quad = opt_settings%w_quad

        ! Set weights for Dipole only
        opt_settings%w_eint = 0.0_wp
        opt_settings%w_dip = 1.0_wp
        opt_settings%w_quad = 0.0_wp

        ! First, get the analytical gradient at the current point
        call calculate_cost_and_gradient(opt_params, molecule, qm_database, opt_settings, cost, analytical_grad, &
                                      FQMolecules, qm_data_idx = 1)

        print *, '--------------------------------------------------'
        print *, 'Dipole Gradient Check'
        print *, '--------------------------------------------------'
        print *, 'Parameter', '      Analytical', '      Numerical', '      Abs. Diff.', '   Rel. Diff.'

        ! Check gradients for chi parameters
        do i = 1, molecule%nAtomTypes
            ! Calculate cost at p+h
            opt_params%chi(i) = opt_params%chi(i) + h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_plus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Calculate cost at p-h
            opt_params%chi(i) = opt_params%chi(i) - 2.0_wp*h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_minus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Reset parameter
            opt_params%chi(i) = opt_params%chi(i) + h

            ! Calculate numerical gradient
            num_grad = (cost_plus-cost_minus) / (2.0_wp*h)

            print '(a, i2, a, 4(2x, es12.4e2))', 'chi(', i, ')', &
                analytical_grad(i), num_grad, &
                abs(analytical_grad(i) - num_grad), &
                abs(analytical_grad(i) - num_grad) / abs(analytical_grad(i))
        enddo

        ! Check gradients for eta parameters
        do i = 1, molecule%nAtomTypes
            ! Calculate cost at p+h
            opt_params%eta(i) = opt_params%eta(i) + h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_plus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Calculate cost at p-h
            opt_params%eta(i) = opt_params%eta(i) - 2.0_wp*h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_minus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Reset parameter
            opt_params%eta(i) = opt_params%eta(i) + h

            ! Calculate numerical gradient
            num_grad = (cost_plus-cost_minus) / (2.0_wp*h)

            print '(a, i2, a, 4(2x, es12.4e2))', 'eta(', i, ')', &
                analytical_grad(molecule%nAtomTypes+i), num_grad, &
                abs(analytical_grad(molecule%nAtomTypes+i) - num_grad), &
                abs(analytical_grad(molecule%nAtomTypes+i) - num_grad) / abs(analytical_grad(molecule%nAtomTypes+i))
        enddo

        ! Restore original weights
        opt_settings%w_eint = original_w_eint
        opt_settings%w_dip = original_w_dip
        opt_settings%w_quad = original_w_quad

        print *, '--------------------------------------------------'

    end subroutine check_dip_gradient

    !> @brief Numerically checks the analytical gradient of the quadrupole moment cost.
    !> @details This subroutine isolates the quadrupole moment component of the cost function
    !>          and compares its analytical gradient with a numerical one computed via finite differences.
    !> @param[inout] opt_params The optimization parameters.
    !> @param[in] molecule The molecular system definition.
    !> @param[in] qm_database The database of reference QM calculations.
    !> @param[inout] opt_settings The settings for the optimization run (weights are modified).
    !> @param[inout] FQMolecules The FQ molecules data structure.
    subroutine check_quad_gradient(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        type(optimization_params_t), intent(inout):: opt_params
        type(molecule_t), intent(in):: molecule
        type(qm_datapoint_t), allocatable, intent(in):: qm_database(:)
        type(optimization_settings_t), intent(inout):: opt_settings
        type(FQMol_t), dimension(:), intent(inout):: FQMolecules

        real(wp), parameter:: h = 1.0e-6
        real(wp):: cost_plus, cost_minus, num_grad
        real(wp):: cost
        real(wp), dimension(2*molecule%nAtomTypes):: analytical_grad, dummy_grad
        integer:: i

        ! Save original weights
        real(wp):: original_w_eint, original_w_dip, original_w_quad
        original_w_eint = opt_settings%w_eint
        original_w_dip = opt_settings%w_dip
        original_w_quad = opt_settings%w_quad

        ! Set weights for Quadrupole only
        opt_settings%w_eint = 0.0_wp
        opt_settings%w_dip = 0.0_wp
        opt_settings%w_quad = 1.0_wp

        ! First, get the analytical gradient at the current point
        call calculate_cost_and_gradient(opt_params, molecule, qm_database, opt_settings, cost, analytical_grad, &
                          FQMolecules, qm_data_idx = 1)

        print *, '--------------------------------------------------'
        print *, 'Quadrupole Gradient Check'
        print *, '--------------------------------------------------'
        print *, 'Parameter', '      Analytical', '      Numerical', '      Abs. Diff.', '   Rel. Diff.'

        ! Check gradients for chi parameters
        do i = 1, molecule%nAtomTypes
            ! Calculate cost at p+h
            opt_params%chi(i) = opt_params%chi(i) + h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_plus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Calculate cost at p-h
            opt_params%chi(i) = opt_params%chi(i) - 2.0_wp*h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_minus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Reset parameter
            opt_params%chi(i) = opt_params%chi(i) + h

            ! Calculate numerical gradient
            num_grad = (cost_plus-cost_minus) / (2.0_wp*h)

            print '(a, i2, a, 4(2x, es12.4e2))', 'chi(', i, ')', &
                analytical_grad(i), num_grad, &
                abs(analytical_grad(i) - num_grad), &
                abs(analytical_grad(i) - num_grad) / abs(analytical_grad(i))
        enddo

        ! Check gradients for eta parameters
        do i = 1, molecule%nAtomTypes
            ! Calculate cost at p+h
            opt_params%eta(i) = opt_params%eta(i) + h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_plus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Calculate cost at p-h
            opt_params%eta(i) = opt_params%eta(i) - 2.0_wp*h
            call calculate_cost_and_gradient(opt_params, molecule, qm_database, &
                opt_settings, cost_minus, dummy_grad, FQMolecules, qm_data_idx = 1)

            ! Reset parameter
            opt_params%eta(i) = opt_params%eta(i) + h

            ! Calculate numerical gradient
            num_grad = (cost_plus-cost_minus) / (2.0_wp*h)

            print '(a, i2, a, 4(2x, es12.4e2))', 'eta(', i, ')', &
                analytical_grad(molecule%nAtomTypes+i), num_grad, &
                abs(analytical_grad(molecule%nAtomTypes+i) - num_grad), &
                abs(analytical_grad(molecule%nAtomTypes+i) - num_grad) / abs(analytical_grad(molecule%nAtomTypes+i))
        enddo

        ! Restore original weights
        opt_settings%w_eint = original_w_eint
        opt_settings%w_dip = original_w_dip
        opt_settings%w_quad = original_w_quad

        print *, '--------------------------------------------------'

    end subroutine check_quad_gradient

    !> @brief Main driver subroutine for running all gradient checks.
    !> @details This routine orchestrates the gradient verification process. It creates the
    !>          FQMolecules data structure and then calls the individual check routines for
    !>          charge derivatives, interaction energy, dipole, and quadrupole gradients.
    !>          Finally, it cleans up the memory.
    !> @param[inout] opt_params The optimization parameters.
    !> @param[in] molecule The molecular system definition.
    !> @param[in] qm_database The database of reference QM calculations.
    !> @param[inout] opt_settings The settings for the optimization run.
    subroutine check_gradients(opt_params, molecule, qm_database, opt_settings)
        type(optimization_params_t), intent(inout):: opt_params
        type(molecule_t), intent(in):: molecule
        type(qm_datapoint_t), allocatable, intent(in):: qm_database(:)
        type(optimization_settings_t), intent(inout):: opt_settings
        type(FQMol_t), dimension(:), allocatable:: FQMolecules

        ! Allocate and initialize the FQMolecules structure
        call CreateFQMolecules(molecule%nMol, molecule%nAtoms, molecule%nAtomsInMol, &
                             molecule%MolID, molecule%MolCharges, molecule%coordinates, &
                             molecule%symbols, FQMolecules)

        ! Call individual gradient checks
        call check_dq_dp_derivatives(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        call check_eint_gradient(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        call check_dip_gradient(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        call check_quad_gradient(opt_params, molecule, qm_database, opt_settings, FQMolecules)

        ! Clean up memory
        call DestroyFQMolecules(FQMolecules)

    end subroutine check_gradients

end module debug_module
