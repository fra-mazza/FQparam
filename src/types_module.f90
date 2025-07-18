module types_module
  implicit none
  integer, parameter :: wp = selected_real_kind(15, 307)

  !> Holds atom-specific information for the FQ calculation
  type :: FQMol_t
    integer :: nAtomsInMol, MolID
    real(kind = wp) :: MolCharge
    character(len = 3), dimension(:), allocatable :: Symbols
    real(kind = wp), dimension(:,:), allocatable :: Coords
    real(kind = wp), dimension(:), allocatable :: chi, eta
    real(kind = wp), dimension(:), allocatable :: FQs
  end type FQMol_t

  !> Holds the data for a single QM reference calculation
  type :: qm_datapoint_t
      real(wp) :: q_int
      real(wp), dimension(3) :: q_coords
      real(wp), dimension(3) :: ext_dipole
      real(wp) :: eint_qm
      real(wp), dimension(3) :: dip_qm
      real(wp), dimension(6) :: quad_qm
  end type qm_datapoint_t

  !> Holds the definition of the molecule being parametrized
  type :: molecule_t
      integer :: nAtoms, nAtomTypes
      integer, dimension(:), allocatable :: atom_types ! map each atom to an atom type index
      character(len=3), dimension(:), allocatable :: symbols
      real(wp), dimension(:,:), allocatable :: coordinates
      integer :: nMol
      integer, dimension(:), allocatable :: nAtomsInMol
      integer, dimension(:), allocatable :: MolID
      real(wp), dimension(:), allocatable :: MolCharges
  end type molecule_t

  !> Holds the parameters for the optimization
  type :: optimization_params_t
      real(wp), dimension(:), allocatable :: chi
      real(wp), dimension(:), allocatable :: eta
  end type optimization_params_t

  !> Holds settings for the optimization run
  type :: optimization_settings_t
      integer :: max_iter
      real(wp) :: tolerance
      real(wp) :: w_eint, w_dip, w_quad
      real(wp), dimension(3) :: dipole_origin
      real(wp), dimension(3) :: quad_origin
  end type optimization_settings_t

end module types_module