module optimizer_module
    use types_module
    use cost_function_module
    use ieee_arithmetic, only: ieee_is_finite
    implicit none

contains

    subroutine optimize_parameters(opt_params, molecule, qm_database, opt_settings, FQMolecules)
        type(optimization_params_t), intent(inout):: opt_params
        type(molecule_t), intent(in):: molecule
        type(qm_datapoint_t), allocatable, intent(in):: qm_database(:)
        type(optimization_settings_t), intent(in):: opt_settings
        type(FQMol_t), dimension(:), intent(inout):: FQMolecules

        integer:: it, ls_iter, n_params
        real(wp):: cost, cost_old
        real(wp), dimension(:), allocatable:: g, g_old, d, p, p_trial
        real(wp):: alpha
        real(wp), parameter:: c1 = 1.0e-4_wp  ! for Armijo condition

        n_params = 2*molecule%nAtomTypes
        allocate(g(n_params), g_old(n_params), d(n_params), p(n_params), p_trial(n_params))

        p(1:molecule%nAtomTypes) = opt_params%chi
        p(molecule%nAtomTypes+1:n_params) = opt_params%eta

        call calculate_cost_and_gradient(opt_params, molecule, qm_database, opt_settings, cost, g, FQMolecules)
        print *, 'Initial Cost (Steepest Descent): ', cost

        do it = 1, opt_settings%max_iter
            if (any(.not. ieee_is_finite(g))) then
                print *, 'Terminated: Gradient is NaN/Inf at iteration ', it
                return
            endif

            g_old = g
            d = -g_old

            cost_old = cost

            ! --- Backtracking line search to find step size alpha ---
            alpha = 1.0_wp
            do ls_iter = 1, 20
                p_trial = p + alpha*d
                opt_params%chi = p_trial(1:molecule%nAtomTypes)
                opt_params%eta = p_trial(molecule%nAtomTypes+1:n_params)

                ! Calculate cost and gradient at the trial point
                call calculate_cost_and_gradient(opt_params, molecule, qm_database, opt_settings, cost, g, FQMolecules)

                ! Check Armijo condition for sufficient decrease
                if (cost <= cost_old + c1*alpha*dot_product(d, g_old)) then
                    exit  ! Condition satisfied, accept step
                endif
                alpha = alpha * 0.5 ! Reduce step size
            enddo

            ! Accept the last tried step, even if it didn't satisfy the Armijo condition.
            ! The cost and gradient are already updated from the last call inside the loop.
            p = p_trial

            print *, 'Iteration: ', it, 'Cost: ', cost, 'alpha: ', alpha

            if (sqrt(dot_product(g, g)) < opt_settings%tolerance) then
                print *, 'Convergence reached after ', it, ' iterations.'
                return
            endif
        enddo
        print *, 'Maximum iterations reached.'
    end subroutine optimize_parameters

end module optimizer_module