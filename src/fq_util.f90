!> @brief Utility module for the Fluctuating Charge (FQ) model.
!> @details This module provides the core subroutines for FQ calculations, including
!>          the construction of the FQ interaction matrix (Jacobian), updating the fluctuating
!>          charges, calculating electrostatic properties (energy, dipole, quadrupole), and
!>          computing the derivatives of charges with respect to model parameters (chi, eta).
!>          It also handles the memory management for the FQMol_t data structure.
Module fq_util

use iso_fortran_env, only: real64, int8
use misc
use types_module, only: FQMol_t, wp

implicit none
save

contains

  !> @brief Allocates and initializes the static parts of the FQMolecules array.
  !> @details This subroutine should be called only ONCE at the beginning of a calculation.
  !>          It allocates the main FQMolecules array and the memory for components that
  !>          do not change during the optimization, such as coordinates and atom symbols.
  !>          The parameter-dependent arrays (chi, eta) and charge arrays (FQs) are allocated
  !>          but not initialized here.
  !> @param[in] nFQMol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] nAtomsInMol Array containing the number of atoms in each molecule.
  !> @param[in] MolID Array mapping each atom to its molecule ID.
  !> @param[in] Charges Array of total charges for each molecule.
  !> @param[in] Coords Cartesian coordinates for all atoms.
  !> @param[in] Symbols Element symbols for all atoms.
  !> @param[out] FQMolecules The freshly allocated and partially initialized array of FQ molecules.
  subroutine CreateFQMolecules(nFQMol, nFQatoms,  nAtomsInMol, MolID, Charges, Coords, Symbols, FQMolecules)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms
    integer, dimension(nFQMol),  intent(in):: nAtomsInMol
    integer, dimension(nFQAtoms),  intent(in):: MolID
    real(kind = wp), dimension(nFQMol), intent(in):: Charges
    real(kind = wp), dimension(3, nFQAtoms), intent(in):: Coords
    character(len = 3), dimension(nFQAtoms), intent(in):: Symbols
    type(FQMol_t), dimension(:), intent(out), allocatable:: FQMolecules
    integer:: i, FirstAtom, LastAtom, Offset

    allocate(FQMolecules(nFQMol))
    Offset = 0
    do i = 1, nFQMol
      FirstAtom = Offset+1
      LastAtom = Offset+nAtomsInMol(i)
      FQMolecules(i)%nAtomsInMol = nAtomsInMol(i)
      FQMolecules(i)%MolID = MolID(FirstAtom)
      FQMolecules(i)%MolCharge = Charges(i)

      allocate(FQMolecules(i)%Coords(3, nAtomsInMol(i)))
      FQMolecules(i)%Coords = Coords(:,FirstAtom:LastAtom)

      allocate(FQMolecules(i)%Symbols(nAtomsInMol(i)))
      FQMolecules(i)%Symbols = Symbols(FirstAtom:LastAtom)

      allocate(FQMolecules(i)%chi(nAtomsInMol(i)))
      allocate(FQMolecules(i)%eta(nAtomsInMol(i)))
      allocate(FQMolecules(i)%FQs(nAtomsInMol(i)))

      Offset = LastAtom
    end do
  end subroutine CreateFQMolecules

  !> @brief Deallocates the FQMolecules array and all its allocatable members.
  !> @details This subroutine should be called at the end of the program to ensure
  !>          proper memory cleanup and prevent memory leaks.
  !> @param[inout] FQMolecules The FQ molecules array to be deallocated.
  subroutine DestroyFQMolecules(FQMolecules)
    implicit none
    type(FQMol_t), dimension(:), allocatable, intent(inout):: FQMolecules
    integer :: i

    if (.not. allocated(FQMolecules)) return

    do i = 1, size(FQMolecules)
        if (allocated(FQMolecules(i)%Symbols)) deallocate(FQMolecules(i)%Symbols)
        if (allocated(FQMolecules(i)%Coords))  deallocate(FQMolecules(i)%Coords)
        if (allocated(FQMolecules(i)%chi))     deallocate(FQMolecules(i)%chi)
        if (allocated(FQMolecules(i)%eta))     deallocate(FQMolecules(i)%eta)
        if (allocated(FQMolecules(i)%FQs))     deallocate(FQMolecules(i)%FQs)
    end do
    deallocate(FQMolecules)
  end subroutine DestroyFQMolecules

  !> @brief Updates the FQ parameters (chi, eta) and resets charges in the FQMolecules structure.
  !> @details This routine is called within the optimization loop whenever the optimizable
  !>          parameters change. It updates the electronegativity and hardness values for each atom.
  !>          It also resets the fluctuating charges to zero before a new calculation.
  !> @param[in] nAtoms Total number of atoms.
  !> @param[in] atom_types Array mapping each atom to its type index.
  !> @param[in] opt_params The current set of optimization parameters (chi and eta values for each type).
  !> @param[inout] FQMolecules The FQ molecules array to be updated.
  subroutine UpdateFQParameters(nAtoms, atom_types, opt_params, FQMolecules)
      use types_module, only: optimization_params_t
      implicit none
      integer, intent(in) :: nAtoms
      integer, dimension(nAtoms), intent(in) :: atom_types
      type(optimization_params_t), intent(in) :: opt_params
      type(FQMol_t), dimension(:), intent(inout) :: FQMolecules

      integer :: iMol, iAtom, atom_idx

      atom_idx = 0
      do iMol = 1, size(FQMolecules)
          do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
              atom_idx = atom_idx + 1
              FQMolecules(iMol)%chi(iAtom) = opt_params%chi(atom_types(atom_idx))
              FQMolecules(iMol)%eta(iAtom) = opt_params%eta(atom_types(atom_idx))
          end do
          FQMolecules(iMol)%FQs = 0.0_wp
      end do
  end subroutine UpdateFQParameters

  !> @brief Constructs the FQ matrix J and its inverse J^-1.
  !> @details This is a critical and computationally expensive step. The J matrix represents the
  !>          linear response of the energy to changes in charge. It includes the atomic hardness
  !>          terms on the diagonal and screened Coulomb interactions on the off-diagonal.
  !>          The matrix is augmented with Lagrange multiplier constraints to enforce charge
  !>          conservation on each molecule.
  !> @param[in] nFQMol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] FQMolecules The array of FQ molecules, containing current eta values and coordinates.
  !> @param[out] FQJinv The calculated inverse of the FQ Jacobian matrix.
  subroutine MakeFQJ(nFQMol, nFQAtoms, FQMolecules, FQJinv)
    implicit none
    type(FQMol_t), dimension(:), intent(in):: FQMolecules
    integer, intent(in):: nFQMol, nFQAtoms
    integer:: i, j, iAtom, iMol, jAtom, jMol
    real(kind = wp), dimension(3):: R_i, R_j
    real(kind = wp), dimension(:,:), intent(out), allocatable:: FQJinv
    real(kind = wp), dimension(:,:), allocatable:: FQJ
    real(kind = wp):: eta_ij, R_ij

    allocate(FQJ(nFQAtoms+nFQMol, nFQAtoms+nFQMol))
    allocate(FQJinv(nFQAtoms+nFQMol, nFQAtoms+nFQMol))
    FQJ = 0.0D0

    ! J_ii =  eta_i
    ! J_ij = eta_ij / (1+eta_ij^2 r_ij^2)^(1/2) with eta_ij = (eta_i+eta_j) / 2
    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        j = 0
        R_i = FQMolecules(iMol)%Coords(:,iAtom)
        do jMol = 1, nFQMol
          do jAtom = 1, FQMolecules(jMol)%nAtomsInMol
            j = j+1
            if (i == j) then
              FQJ(i, j) = FQMolecules(iMol)%eta(iAtom)
            else
              eta_ij = (FQMolecules(iMol)%eta(iAtom) + FQMolecules(jMol)%eta(jAtom))/2.0D0
              R_j =  FQMolecules(jMol)%Coords(:,jAtom)
              R_ij = Sqrt((R_i(1)-R_j(1))**2.0D0 + (R_i(2)-R_j(2))**2.0D0 + (R_i(3)-R_j(3))**2.0D0)
              FQJ(i, j) = eta_ij/Sqrt((1.0D0 + (eta_ij**2 * R_ij**2)))
            end if
          end do
        end do
      end do
    end do

    ! Add charge constraints
    i = 0
    do iMol = 1, nFQMol
      do iAtom =  1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        FQJ(i, nFQAtoms+iMol) = 1.0D0
        FQJ(nFQAtoms+iMol, i) = 1.0D0
      end do
    end do

    FQJinv = inv(FQJ)

  end subroutine MakeFQJ


  !> @brief Solves the FQ linear system to find the fluctuating charges.
  !> @details This subroutine solves the equation `J * q = -chi - V`, where J is the FQ matrix,
  !>          q are the charges, chi are the electronegativities, and V is the external potential.
  !>          It uses the pre-computed inverse of J (FQJinv) for efficiency.
  !> @param[in] nFQMol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[inout] FQMolecules The FQ molecules array. FQs will be updated.
  !> @param[in] FQJinv The inverse of the FQ Jacobian matrix.
  !> @param[in] V The external electrostatic potential at each atomic site.
  Subroutine  UpdateFQs(nFQMol, nFQAtoms, FQMolecules, FQJinv, V)
    implicit none
    type(FQMol_t), dimension(:), intent(inout):: FQMolecules
    integer, intent(in):: nFQMol, nFQAtoms
    real(kind = wp), dimension(:,:), intent(in):: FQJinv
    real(kind = wp), dimension(:), intent(in):: V

    real(kind = wp), dimension(nFQMol+nFQatoms):: RHS
    real(kind = wp), dimension(nFQMol+nFQatoms):: q
    integer:: i, j, iAtom, iMol, jAtom, jMol


    ! Compute the right hand side
    RHS = 0.0_wp

    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        RHS(i) = -FQMolecules(iMol)%Chi(iAtom) - V(i)
      end do
    end do

    do iMol = 1, nFQMol
        i = i+1
        RHS(i) = -FQMolecules(iMol)%MolCharge
    enddo

    ! q = FQJinv @ RHS
    q = 0.0_wp
    call DGEMV('N', nFQMol+nFQAtoms, nFQMol+nFQAtoms, 1.0_wp, &
                 FQJinv, nFQMol+nFQAtoms, RHS, 1, 0.0_wp, q, 1 )

    ! Update FQs in FQMolecules
    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        FQMolecules(iMol)%FQs(iAtom) = q(i)
      end do
    end do

  End Subroutine UpdateFQs

  !> @brief Calculates the FQ interaction energy.
  !> @details The interaction energy is calculated as the sum over all atoms of q_i * V_i,
  !>          where q_i is the fluctuating charge and V_i is the external potential.
  !> @param[in] nFQMol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] FQMolecules The FQ molecules array, containing the calculated FQs.
  !> @param[in] V The external electrostatic potential.
  !> @return The total FQ interaction energy.
  real(wp) function  IntEnergy(nFQMol, nFQatoms, FQMolecules, V)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol_t), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:), intent(in):: V
    integer:: i, iAtom, iMol

    IntEnergy = 0.0_wp
    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        IntEnergy = IntEnErgy+FQMolecules(iMol)%FQs(iAtom)*V(i)
      end do
    end do

  End Function IntEnergy

  !> @brief Calculates the external potential from a point charge and a point dipole.
  !> @details This subroutine computes the electrostatic potential at each atomic site arising
  !>          from the external QM perturbation, which consists of a point charge and a point dipole.
  !> @param[in] nFQMol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] FQMolecules The FQ molecules array (for atomic coordinates).
  !> @param[inout] V The potential array to be filled.
  !> @param[in] qint The magnitude of the external point charge.
  !> @param[in] qcoords The coordinates of the external point charge.
  !> @param[in] ext_dipole The vector of the external point dipole.
  Subroutine  GetPotentialPointCharge(nFQMol, nFQAtoms, FQMolecules, V, qint, qcoords, ext_dipole)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol_t), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:), intent(inout):: V
    real(kind = wp), dimension(3), intent(in):: qcoords
    real(kind = wp),  intent(in):: qint
    real(kind = wp), dimension(3), intent(in):: ext_dipole
    integer:: i, iAtom, iMol
    real(wp):: R, R_cubed
    real(wp), dimension(3):: R_vec

    V = 0.0_wp

    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        R_vec = FQMolecules(iMol)%Coords(:,iAtom) - qcoords
        R = sqrt(sum(R_vec**2))
        R_cubed = R**3

        ! Contribution from point charge
        V(i) = V(i) + qint/R

        ! Contribution from point dipole
        V(i) = V(i) + dot_product(ext_dipole, R_vec) / R_cubed

      end do
    end do
  End Subroutine GetPotentialPointCharge

  !> @brief Calculates the derivatives of fluctuating charges with respect to FQ parameters.
  !> @details This is a key routine for gradient-based optimization. It computes d(q_i)/d(chi_k)
  !>          and d(q_i)/d(eta_k) for all atoms i and atom types k. These derivatives are essential
  !>          for calculating the gradient of the cost function.
  !> @param[in] nFQmol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] nAtomtypes Number of unique atom types.
  !> @param[in] atomtypes Array mapping each atom to its type index.
  !> @param[in] FQMolecules The FQ molecules array.
  !> @param[in] FQJinv The inverse of the FQ Jacobian matrix.
  !> @param[in] V The external potential.
  !> @param[out] dqdchi The calculated derivatives w.r.t. electronegativity.
  !> @param[out] dqdeta The calculated derivatives w.r.t. self-hardness.
  Subroutine  GetFQsDerivatives(nFQmol, nFQAtoms, nAtomtypes,  atomtypes, FQMolecules, FQJinv, V, dqdchi, dqdeta)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms, nAtomtypes
    integer, dimension(:), intent(in):: atomtypes
    type(FQMol_t), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:,:), intent(in):: FQJinv
    real(kind = wp), dimension(:), intent(in):: V
    real(kind = wp), dimension(:,:), intent(out):: dqdchi, dqdeta
    real(kind = wp), dimension(nFQMol+nFQatoms):: RHS, q, dq_dp_temp
    integer:: i, j, iAtom, iMol, jAtom, jMol, jat, type_A
    real(wp):: temp_val
    real(wp), dimension(nFQatoms+nFQMol, nFQatoms+nFQMol):: dFQJ_deta_typeA_matrix
    real(wp), dimension(nFQatoms+nFQMol):: temp_vec_full
    real(wp):: d_eta_ij_d_eta_typeA, eta_ij, R_ij
    real(kind = wp), dimension(3):: R_i, R_j


    ! Compute the right hand side
    RHS = 0.0_wp

    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        RHS(i) = -FQMolecules(iMol)%Chi(iAtom) - V(i)
      end do
    end do

    do iMol = 1, nFQMol
        i = i+1
        RHS(i) = -FQMolecules(iMol)%MolCharge
    enddo


    ! q = D^{-1}_{ij} * RHS{j}
    q = 0.0_wp
    call DGEMV('N', nFQMol+nFQAtoms, nFQMol+nFQAtoms, 1.0_wp, &
                 FQJinv, nFQMol+nFQAtoms, RHS, 1, 0.0_wp, q, 1 )

    ! dFQ/dchi

    dqdchi = 0.0_wp
    do i = 1, nFQatoms  ! q
      do j = 1, nFQatoms  ! chi
        jat  = atomtypes(j)
        dqdchi(i, jat) = dqdchi(i, jat) - FQJinv(i, j)
      enddo
    enddo

    ! dFQ/deta
    ! NB: this is valid only for the ohno kernel

    dqdeta = 0.0_wp

    do type_A = 1, nAtomtypes
      dFQJ_deta_typeA_matrix = 0.0_wp

      ! Populate dFQJ_deta_typeA_matrix
      i = 0
      do iMol = 1, nFQMol
        do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
          i = i+1
          j = 0
          do jMol = 1, nFQMol
            do jAtom = 1, FQMolecules(jMol)%nAtomsInMol
              j = j+1
              if (i == j) then
                ! Diagonal elements: d(eta_i)/d(eta_type_A)
                if (atomtypes(i) == type_A) then
                  dFQJ_deta_typeA_matrix(i, i) = 1.0_wp
                endif
              else
                ! Off-diagonal elements: d(eta_ij / (1+eta_ij^2 r_ij^2)^(1/2)) / d(eta_type_A)
                eta_ij = (FQMolecules(iMol)%eta(iAtom) + FQMolecules(jMol)%eta(jAtom))/2.0D0
                R_i = FQMolecules(iMol)%Coords(:,iAtom)
                R_j = FQMolecules(jMol)%Coords(:,jAtom)
                R_ij = Sqrt((R_i(1)-R_j(1))**2.0D0 + (R_i(2)-R_j(2))**2.0D0 + (R_i(3)-R_j(3))**2.0D0)
                temp_val = (1.0D0 + (eta_ij**2 * R_ij**2))**(-1.5)

                d_eta_ij_d_eta_typeA = 0.0_wp
                if (atomtypes(i) == type_A) then
                  d_eta_ij_d_eta_typeA = d_eta_ij_d_eta_typeA+0.5_wp
                endif
                if (atomtypes(j) == type_A) then
                  d_eta_ij_d_eta_typeA = d_eta_ij_d_eta_typeA+0.5_wp
                endif
                dFQJ_deta_typeA_matrix(i, j) = d_eta_ij_d_eta_typeA*temp_val
              endif
            end do
          end do
        end do
      end do

      ! Calculate temp_vec_full = (dFQJ_deta_typeA_matrix*q)
      call DGEMV('N', nFQatoms+nFQMol, nFQatoms+nFQMol, 1.0_wp, &
                   dFQJ_deta_typeA_matrix, nFQatoms+nFQMol, q, 1, 0.0_wp, temp_vec_full, 1)

      ! Calculate dq/d(eta_typeA) = - J^-1 * (dJ/d(eta_typeA) * q) into a temporary vector
      call DGEMV('N', nFQatoms+nFQMol, nFQatoms+nFQMol, -1.0_wp, &
                   FQJinv, nFQatoms+nFQMol, temp_vec_full, 1, 0.0_wp, dq_dp_temp, 1)

      ! Copy the first nFQatoms elements (the charge derivatives) to the output array
      dqdeta(:, type_A) = dq_dp_temp(1:nFQatoms)

    enddo

  End Subroutine GetFQsDerivatives


  !> @brief Calculates the total dipole moment of the system.
  !> @param[in] nFQmol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] FQMolecules The FQ molecules array, containing charges and coordinates.
  !> @param[in] DipoleOrigin The origin point for the dipole calculation.
  !> @return A 3-element array containing the (x, y, z) components of the dipole moment.
  function Dipole(nFQmol, nFQAtoms, FQMolecules, DipoleOrigin) 
    implicit none
    real(wp), dimension(3):: Dipole
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol_t), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(3), intent(in):: DipoleOrigin
    integer:: i, iAtom, iMol

    Dipole = 0.0_wp

    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1

        Dipole = Dipole+FQMolecules(iMol)%FQs(iAtom)*(FQMolecules(iMol)%Coords(:,iAtom)-DipoleOrigin)

      end do
    end do

  end function Dipole

  !> @brief Calculates the total quadrupole moment of the system.
  !> @param[in] nFQmol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] FQMolecules The FQ molecules array, containing charges and coordinates.
  !> @param[in] QuadOrigin The origin point for the quadrupole calculation.
  !> @return A 6-element array containing the (xx, xy, xz, yy, yz, zz) components of the quadrupole moment.
  function Quadrupole(nFQmol, nFQAtoms, FQMolecules, QuadOrigin) 
    implicit none
    real(wp), dimension(6):: Quadrupole
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol_t), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(3), intent(in):: QuadOrigin
    integer:: i, iAtom, iMol
    real(wp), dimension(3):: R

    Quadrupole = 0.0_wp

    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        R = FQMolecules(iMol)%Coords(:,iAtom) - QuadOrigin

        ! Qxx
        Quadrupole(1) = Quadrupole(1) + FQMolecules(iMol)%FQs(iAtom)*R(1)*R(1)
        ! Qxy
        Quadrupole(2) = Quadrupole(2) + FQMolecules(iMol)%FQs(iAtom)*R(1)*R(2)
        ! Qxz
        Quadrupole(3) = Quadrupole(3) + FQMolecules(iMol)%FQs(iAtom)*R(1)*R(3)
        ! Qyy
        Quadrupole(4) = Quadrupole(4) + FQMolecules(iMol)%FQs(iAtom)*R(2)*R(2)
        ! Qyz
        Quadrupole(5) = Quadrupole(5) + FQMolecules(iMol)%FQs(iAtom)*R(2)*R(3)
        ! Qzz
        Quadrupole(6) = Quadrupole(6) + FQMolecules(iMol)%FQs(iAtom)*R(3)*R(3)

      end do
    end do

  end function Quadrupole

  !> @brief Updates the gradient vector with contributions from all error terms.
  !> @details This subroutine calculates the derivative of the total cost function with respect to
  !>          each of the optimizable parameters (chi and eta for each atom type). It does this by
  !>          applying the chain rule, combining the derivatives of the error terms (e.g., d(cost_Eint)/d(EintFQ))
  !>          with the derivatives of the FQ properties w.r.t. the charges (e.g., d(EintFQ)/d(q))
  !>          and the derivatives of the charges w.r.t. the parameters (d(q)/d(p), from GetFQsDerivatives).
  !> @param[in] nFQMol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[in] nAtomtypes Number of unique atom types.
  !> @param[in] atomtypes Array mapping each atom to its type index.
  !> @param[in] FQMolecules The FQ molecules array.
  !> @param[in] V The external potential.
  !> @param[in] FQJinv The inverse of the FQ Jacobian matrix.
  !> @param[in] wEint Weight for the interaction energy error.
  !> @param[in] EintFQ The calculated FQ interaction energy.
  !> @param[in] EintQM The reference QM interaction energy.
  !> @param[in] wdip Weight for the dipole error.
  !> @param[in] dipFQ The calculated FQ dipole.
  !> @param[in] QMdip The reference QM dipole.
  !> @param[in] wquad Weight for the quadrupole error.
  !> @param[in] quadFQ The calculated FQ quadrupole.
  !> @param[in] QMquad The reference QM quadrupole.
  !> @param[in] DipoleOrigin The origin for the dipole calculation.
  !> @param[in] QuadOrigin The origin for the quadrupole calculation.
  !> @param[inout] gradients The gradient vector to be updated.
  Subroutine  UpdateGradients(nFQMol, nFQAtoms, nAtomTypes, atomtypes,  FQMolecules, V, FQJinv, wEint, EintFQ, EintQM, &
      wdip, dipFQ, QMdip, wquad, quadFQ, QMquad, DipoleOrigin, QuadOrigin, gradients)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms, nAtomtypes
    integer, dimension(:), intent(in):: atomtypes
    type(FQMol_t), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:,:), intent(in):: FQJinv
    real(kind = wp), dimension(:), intent(in):: V
    real(kind = wp),  intent(in):: wEint, EintFQ, EintQM, wdip,  wquad
    real(kind = wp), dimension(3),  intent(in):: dipFQ, QMdip, DipoleOrigin, QuadOrigin
    real(kind = wp), dimension(6),  intent(in):: quadFQ, QMquad
    real(kind = wp), dimension(:), intent(inout):: gradients

    real(kind = wp), dimension(nFQAtoms, nAtomtypes)  :: dqdchi, dqdeta
    integer:: i, j, iAtom, iMol, jAtom, jMol, k
    real(wp), dimension(3):: R

    dqdchi = 0.0_wp
    dqdeta = 0.0_wp

    ! compute derivatives of charges
    call  GetFQsDerivatives(nFQmol, nFQAtoms, nAtomtypes,  atomtypes, FQMolecules, FQJinv, V, dqdchi, dqdeta)

    ! compute dF/dchi and dF/deta
    do i = 1, nAtomtypes
      ! Eint contribution
      do j = 1, nFQAtoms
        gradients(i) = gradients(i) + 2*wEint*(EintFQ-EintQM)*dqdchi(j, i)*V(j)
        gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wEint*(EintFQ-EintQM)*dqdeta(j, i)*V(j)
      enddo

      ! dipole contribution
      j = 0
      do jMol = 1, nFQMol
        do jAtom = 1, FQMolecules(jMol)%nAtomsInMol
          j = j+1
          R = FQMolecules(jMol)%Coords(:,jAtom) - DipoleOrigin
          do k = 1, 3
            gradients(i) = gradients(i) + 2*wdip*(dipFQ(k) - QMdip(k))*dqdchi(j, i)*R(k)
            gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wdip*(dipFQ(k) - QMdip(k))*dqdeta(j, i)*R(k)
          enddo
        enddo
      enddo

      ! quadrupole contribution
      j = 0
      do jMol = 1, nFQMol
        do jAtom = 1, FQMolecules(jMol)%nAtomsInMol
          j = j+1
          R = FQMolecules(jMol)%Coords(:,jAtom) - QuadOrigin
          ! Qxx
          gradients(i) = gradients(i) + 2*wquad*(quadFQ(1) - QMquad(1))*dqdchi(j, i)*R(1)*R(1)
          gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wquad*(quadFQ(1) - QMquad(1))*dqdeta(j, i)*R(1)*R(1)
          ! Qxy
          gradients(i) = gradients(i) + 2*wquad*(quadFQ(2) - QMquad(2))*dqdchi(j, i)*R(1)*R(2)
          gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wquad*(quadFQ(2) - QMquad(2))*dqdeta(j, i)*R(1)*R(2)
          ! Qxz
          gradients(i) = gradients(i) + 2*wquad*(quadFQ(3) - QMquad(3))*dqdchi(j, i)*R(1)*R(3)
          gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wquad*(quadFQ(3) - QMquad(3))*dqdeta(j, i)*R(1)*R(3)
          ! Qyy
          gradients(i) = gradients(i) + 2*wquad*(quadFQ(4) - QMquad(4))*dqdchi(j, i)*R(2)*R(2)
          gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wquad*(quadFQ(4) - QMquad(4))*dqdeta(j, i)*R(2)*R(2)
          ! Qyz
          gradients(i) = gradients(i) + 2*wquad*(quadFQ(5) - QMquad(5))*dqdchi(j, i)*R(2)*R(3)
          gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wquad*(quadFQ(5) - QMquad(5))*dqdeta(j, i)*R(2)*R(3)
          ! Qzz
          gradients(i) = gradients(i) + 2*wquad*(quadFQ(6) - QMquad(6))*dqdchi(j, i)*R(3)*R(3)
          gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wquad*(quadFQ(6) - QMquad(6))*dqdeta(j, i)*R(3)*R(3)
        enddo
      enddo
    enddo


  End Subroutine UpdateGradients

  !> @brief Calculates and prints the fluctuating charges for the unperturbed system.
  !> @details This subroutine is used to see the intrinsic charge distribution of the molecule
  !>          in the absence of any external electric field. It sets the external potential V to zero
  !>          and solves for the charges.
  !> @param[in] nFQMol Total number of molecules.
  !> @param[in] nFQAtoms Total number of atoms.
  !> @param[inout] FQMolecules The FQ molecules array. The FQs will be updated and printed.
  !> @param[in] FQJinv The inverse of the FQ Jacobian matrix.
  subroutine calculate_and_print_unperturbed_charges(nFQMol, nFQAtoms, FQMolecules, FQJinv)
    implicit none
    integer, intent(in) :: nFQMol, nFQAtoms
    type(FQMol_t), dimension(:), intent(inout) :: FQMolecules
    real(kind=wp), dimension(:,:), intent(in) :: FQJinv
    real(kind=wp), dimension(nFQAtoms) :: V
    integer :: i, iMol, iAtom

    ! Set external potential to zero
    V = 0.0_wp

    ! Solve for the fluctuating charges
    call UpdateFQs(nFQMol, nFQAtoms, FQMolecules, FQJinv, V)

    ! Print the charges
    print *, '--------------------------------------------------'
    print *, ' Unperturbed Fluctuating Charges'
    print *, '--------------------------------------------------'
    print *, ' Atom  Symbol   Charge (a.u.)'
    i = 0
    do iMol = 1, nFQMol
        do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
            i = i + 1
            print '(i5, 4x, a3, 4x, f12.6)', i, FQMolecules(iMol)%Symbols(iAtom), FQMolecules(iMol)%FQs(iAtom)
        end do
    end do
    print *, '--------------------------------------------------'

  end subroutine calculate_and_print_unperturbed_charges

End Module  fq_util