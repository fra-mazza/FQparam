!> @brief Module for calculating the cost function and its gradient.
!> @details This module contains the primary subroutine that evaluates the difference
!>          between the FQ model's predictions and the reference QM data. It computes
!>          the total cost (a sum of squared errors for energy, dipole, and quadrupole)
!>          and the analytical gradient of this cost with respect to the optimizable
!>          parameters (chi and eta).
module cost_function_module
    use types_module
    use fq_util
    implicit none

contains

    !> @brief Calculates the total cost and its gradient for a given set of parameters.
    !> @details This is the core of the optimization. It iterates through the QM database,
    !>          and for each data point, it calculates the FQ properties (interaction energy,
    !>          dipole, quadrupole), compares them to the QM reference values, and accumulates
    !>          the cost and its gradient. The calculation is optimized by computing the expensive
    !>          Jacobian inverse only once per call, as it only depends on the parameters, not
    !>          the individual QM data point.
    !> @param[in] opt_params The current set of optimization parameters (chi, eta).
    !> @param[in] molecule The molecular system definition.
    !> @param[in] qm_database The database of reference QM calculations.
    !> @param[in] opt_settings The settings for the optimization run.
    !> @param[out] cost The calculated total cost.
    !> @param[out] gradient The calculated gradient of the cost function.
    !> @param[inout] FQMolecules The FQ molecules data structure.
    !> @param[in] qm_data_idx Optional index to calculate the cost for only a single QM data point.
    subroutine calculate_cost_and_gradient(opt_params, molecule, qm_database, opt_settings, &
                                           cost, gradient, FQMolecules, qm_data_idx)
        type(optimization_params_t), intent(in):: opt_params
        type(molecule_t), intent(in):: molecule
        type(qm_datapoint_t), allocatable, intent(in):: qm_database(:)
        type(optimization_settings_t), intent(in):: opt_settings
        real(wp), intent(out):: cost
        real(wp), dimension(:), intent(out):: gradient
        type(FQMol_t), dimension(:), intent(inout):: FQMolecules
        integer, optional, intent(in):: qm_data_idx

        real(wp):: EintFQ
        real(wp), dimension(3):: dipFQ
        real(wp), dimension(6):: QuadFQ
        real(wp), dimension(molecule%nAtoms):: V
        real(wp), dimension(:,:), allocatable:: FQJinv
        integer:: i, start_idx, end_idx
        real(wp):: cost_eint, cost_dip, cost_quad

        cost = 0.0_wp
        gradient = 0.0_wp

        if (present(qm_data_idx)) then
            start_idx = qm_data_idx
            end_idx = qm_data_idx
        else
            start_idx = 1
            end_idx = size(qm_database)
        endif

        ! Update FQ parameters (chi, eta) and reset charges
        call UpdateFQParameters(molecule%nAtoms, molecule%atom_types, opt_params, FQMolecules)

        ! Calculate the Jacobian and its inverse once for the current set of parameters
        call MakeFQJ(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv)

        ! Loop over all QM data points to calculate the total cost and gradient
        do i = start_idx, end_idx
            ! Get the external potential from the current point charge
            call GetPotentialPointCharge(molecule%nMol, molecule%nAtoms, FQMolecules, &
                V, qm_database(i)%q_int, qm_database(i)%q_coords, qm_database(i)%ext_dipole)

            ! Solve for the fluctuating charges
            call UpdateFQs(molecule%nMol, molecule%nAtoms, FQMolecules, FQJinv, V)

            ! Calculate FQ properties
            EintFQ = intenergy(molecule%nMol, molecule%nAtoms, FQMolecules, V)
            dipFQ = Dipole(molecule%nMol, molecule%nAtoms, FQMolecules, opt_settings%dipole_origin)
            quadFQ = Quadrupole(molecule%nMol, molecule%nAtoms, FQMolecules, opt_settings%quad_origin)

            ! Accumulate the cost function value
            cost_eint = opt_settings%w_eint*(EintFQ-qm_database(i)%eint_qm)**2
            cost_dip = opt_settings%w_dip*sum((dipFQ-qm_database(i)%dip_qm)**2)
            cost_quad = opt_settings%w_quad*sum((quadFQ-qm_database(i)%quad_qm)**2)

            cost = cost+cost_eint+cost_dip+cost_quad

            ! Accumulate the gradients
            call UpdateGradients(molecule%nMol, molecule%nAtoms, &
                molecule%nAtomTypes, molecule%atom_types, FQMolecules, V, FQJinv, &
                opt_settings%w_eint, EintFQ, qm_database(i)%eint_qm, &
                opt_settings%w_dip, dipFQ, qm_database(i)%dip_qm, &
                opt_settings%w_quad, quadFQ, qm_database(i)%quad_qm, &
                opt_settings%dipole_origin, opt_settings%quad_origin, gradient)
        enddo

        if (allocated(FQJinv)) deallocate(FQJinv)

    end subroutine calculate_cost_and_gradient

end module cost_function_module