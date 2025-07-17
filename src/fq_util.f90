Module  fq_util  !

use iso_fortran_env, only: real64, int8
use misc


implicit none
save

integer, parameter:: wp = real64, iwp = int8

!
type:: FQMol
  integer:: nAtomsInMol, MolID
  real(kind = wp):: MolCharge
  character(len = 3), dimension(:), allocatable:: Symbols
  real(kind = wp), dimension(:,:), allocatable:: Coords
  real(kind = wp), dimension(:), allocatable:: chi, eta
  real(kind = wp), dimension(:), allocatable:: FQs
end type FQMol


!
contains
!
  subroutine InitializeMolecules(nFQMol, nFQatoms,  nAtomsInMol, MolID, Charges, chi, eta, Coords, Symbols, FQs, FQMolecules)
    implicit none
    integer, parameter:: wp = real64, iwp = int8
    integer, intent(in):: nFQMol, nFQAtoms
    integer, dimension(nFQMol),  intent(in):: nAtomsInMol
    integer, dimension(nFQAtoms),  intent(in):: MolID
    real(kind = wp), dimension(nFQAtoms), intent(in):: eta, chi, Charges, FQs
    real(kind = wp), dimension(3, nFQAtoms), intent(in):: Coords
    character(len = 3), dimension(nFQAtoms), intent(in):: Symbols
    type(FQMol), dimension(:), intent(out), allocatable:: FQMolecules
    integer:: i, FirstAtom, LastAtom, Offset
    
    allocate(FQMolecules(nFQMol))
    Offset = 0
    do i = 1, nFQMol
      FirstAtom = Offset+1
      LastAtom = Offset+nAtomsInMol(i)
      FQMolecules(i) = FQMol(nAtomsInMol = nAtomsInMol(i),            &
                             MolID = MolID(i),                        &
                             MolCharge = Charges(i),                  &
                             Coords = Coords(:,FirstAtom:LastAtom),   &
                             chi = chi(FirstAtom:LastAtom),           &
                             eta = eta(FirstAtom:LastAtom),           &
                             Symbols = Symbols(FirstAtom:LastAtom),   &
                             FQs = FQs(FirstAtom:LastAtom))
      Offset = LastAtom
    end do
  end subroutine InitializeMolecules
!
  subroutine MakeFQJ(nFQMol, nFQAtoms, FQMolecules, FQJinv)
!
! Construct the FQ matrix J and invert it to J-1 (used throughout the calculation)
!
    implicit none
    integer, parameter:: wp = real64, iwp = int8
    type(FQMol), dimension(:), intent(in):: FQMolecules
    integer, intent(in):: nFQMol, nFQAtoms
    integer:: i, j, iAtom, iMol, jAtom, jMol
    real(kind = wp), dimension(3):: R_i, R_j
    real(kind = wp), dimension(:,:), intent(out), allocatable:: FQJinv
    real(kind = wp), dimension(:,:), allocatable:: FQJ
    real(kind = wp):: eta_ij, R_ij
    

    !
    allocate(FQJ(nFQAtoms+nFQMol, nFQAtoms+nFQMol))
    allocate(FQJinv(nFQAtoms+nFQMol, nFQAtoms+nFQMol))
    FQJ = 0.0D0
!
!   J_ii =  eta_i
!   J_ij = eta_ij / (1+eta_ij^2 r_ij^2)^(1/2) with eta_ij = (eta_i+eta_j) / 2
!
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
!
!   Add charge constraints
!
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


  Subroutine  UpdateFQs(nFQMol, nFQAtoms, FQMolecules, FQJinv, V)
    implicit none 
    type(FQMol), dimension(:), intent(inout):: FQMolecules
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
    call DGEMV('N', nFQMol+nFQAtoms, nFQMol+nFQAtoms, 1.0_wp, FQJinv, nFQMol+nFQAtoms, RHS, 1, 0.0_wp, q, 1 )

    ! Update FQs in FQMolecules
    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        FQMolecules(iMol)%FQs(iAtom) = q(i)
      end do
    end do




  End Subroutine UpdateFQs 

  real(wp) function  IntEnergy(nFQMol, nFQatoms, FQMolecules, V)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:), intent(in):: V
    integer:: i, j, iAtom, iMol, jAtom, jMol

    IntEnergy = 0.0_wp
    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        IntEnergy = IntEnErgy+FQMolecules(iMol)%FQs(iAtom)*V(i)
      end do
    end do

  End Function IntEnergy 

  Subroutine  GetPotentialPointCharge(nFQMol, nFQAtoms, FQMolecules, V, qint, qcoords)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:), intent(inout):: V
    real(kind = wp), dimension(3), intent(in):: qcoords
    real(kind = wp),  intent(in):: qint
    integer:: i, j, iAtom, iMol, jAtom, jMol, iq, nq
    real(wp):: R

    V = 0.0_wp

    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1
        R = sqrt((qcoords(1) - FQMolecules(iMol)%Coords(1, iAtom))**2 + &
          (qcoords(2) - FQMolecules(iMol)%Coords(2, iAtom))**2 &
        +(qcoords(3) - FQMolecules(iMol)%Coords(3, iAtom))**2) 
        
        V(i) = V(i) + qint/R
      end do
    end do
  End Subroutine GetPotentialPointCharge 

  Subroutine  GetFQsDerivatives(nFQmol, nFQAtoms, nAtomtypes,  atomtypes, FQMolecules, FQJinv, V, dqdchi, dqdmu)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms, nAtomtypes
    integer, dimension(:), intent(in):: atomtypes
    type(FQMol), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:,:), intent(in):: FQJinv
    real(kind = wp), dimension(:), intent(in):: V
    real(kind = wp), dimension(:,:), intent(out):: dqdchi, dqdmu
    real(kind = wp), dimension(nFQMol+nFQatoms):: RHS, q
    integer:: i, j, iAtom, iMol, jAtom, jMol, jat, iat, k
    real(wp):: R, temp
    real(kind = wp), dimension(3):: R_i, R_j
    real(kind = wp):: eta_ij, R_ij


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
    call DGEMV('N', nFQMol+nFQAtoms, nFQMol+nFQAtoms, 1.0_wp, FQJinv, nFQMol+nFQAtoms, RHS, 1, 0.0_wp, q, 1 )

    ! dFQ/dchi

    dqdchi = 0.0_wp
    do i = 1, nFQatoms  ! q
      do j = 1, nFQatoms  ! chi
        jat  = atomtypes(j)
        dqdchi(i, jat) = dqdchi(i, jat) - FQJinv(i, j)
      enddo
    enddo

    ! dFQ/dmu
    ! NB: this is valid only for the ohno kernel

    dqdmu = 0.0_wp
    do k = 1, nFQatoms
      i = 0
      do iMol = 1, nFQMol
        do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
          i = i+1
          iat  = atomtypes(i)
          ! - D^{-1}_{ki} D^{-1}_{ij} RHS_{j} = -D^{-1}_{ki}q{i}
          dqdmu(k, iat) = dqdmu(k, iat) - FQJinv(k, i)*q(i)
          
          j = 0
          R_i = FQMolecules(iMol)%Coords(:,iAtom)
          do jMol = 1, nFQMol
            do jAtom = 1, FQMolecules(jMol)%nAtomsInMol
              j = j+1
              if (i .ne. j) then
                ! -1/2 D^{-1}_{ki} * (1/(1+\eta_{ij}^2 r_{ij}^2)^3/2) D^{-1}_{jm}RHS{m} =
                ! = -1/2 D^{-1}_{ki} * (1/(1+\eta_{ij}^2 r_{ij}^2)^3/2) * q_{j}        


                eta_ij = (FQMolecules(iMol)%eta(iAtom) + FQMolecules(jMol)%eta(jAtom))/2.0D0
                R_j =  FQMolecules(jMol)%Coords(:,jAtom)
                R_ij = Sqrt((R_i(1)-R_j(1))**2.0D0 + (R_i(2)-R_j(2))**2.0D0 + (R_i(3)-R_j(3))**2.0D0)
                temp = (1.0D0 + (eta_ij**2 * R_ij**2))**(-1.5)


                dqdmu(k, iat) = dqdmu(k, iat) - 0.5_wp*FQJinv(k, i) * temp*q(j)
                
                ! -1/2 D^{-1}_{kj} * (1/(1+\eta_{ji}^2 r_{ji}^2)^3/2) D^{-1}_{im}RHS{m} =
                ! = -1/2 D^{-1}_{kj} * (1/(1+\eta_{ji}^2 r_{ji}^2)^3/2) * q_{i}        
                
                dqdmu(k, iat) = dqdmu(k, iat) - 0.5_wp*FQJinv(k, j) * temp*q(i)

              end if
            end do
          end do
        end do
      end do
    enddo


    !print*, dqdchi
    !print*, dqdmu

  
  End Subroutine GetFQsDerivatives 


  function Dipole(nFQmol, nFQAtoms, FQMolecules, DipoleOrigin) 
    implicit none
    real(wp), dimension(3):: Dipole
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(3), intent(in):: DipoleOrigin
    integer:: i, j, iAtom, iMol, jAtom, jMol, iq, nq

    Dipole = 0.0_wp

    i = 0
    do iMol = 1, nFQMol
      do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
        i = i+1

        Dipole = Dipole+FQMolecules(iMol)%FQs(iAtom)*(FQMolecules(iMol)%Coords(:,iAtom)-DipoleOrigin)

      end do
    end do
    
  end function Dipole
  
  function Quadrupole(nFQmol, nFQAtoms, FQMolecules, QuadOrigin) 
    implicit none
    real(wp), dimension(6):: Quadrupole
    integer, intent(in):: nFQMol, nFQAtoms
    type(FQMol), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(3), intent(in):: QuadOrigin
    integer:: i, j, iAtom, iMol, jAtom, jMol, iq, nq
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

  Subroutine  UpdateGradients(nFQMol, nFQAtoms, nAtomTypes, atomtypes,  FQMolecules, V, FQJinv, wEint, EintFQ, EintQM, &
      wdip, dipFQ, QMdip, wquad, quadFQ, QMquad, DipoleOrigin, QuadOrigin, gradients)
    implicit none
    integer, intent(in):: nFQMol, nFQAtoms, nAtomtypes
    integer, dimension(:), intent(in):: atomtypes
    type(FQMol), dimension(:), intent(in):: FQMolecules
    real(kind = wp), dimension(:,:), intent(in):: FQJinv
    real(kind = wp), dimension(:), intent(in):: V
    real(kind = wp),  intent(in):: wEint, EintFQ, EintQM, wdip,  wquad
    real(kind = wp), dimension(3),  intent(in):: dipFQ, QMdip, DipoleOrigin, QuadOrigin
    real(kind = wp), dimension(6),  intent(in):: quadFQ, QMquad
    real(kind = wp), dimension(:), intent(inout):: gradients

    real(kind = wp), dimension(nFQAtoms, nAtomtypes)  :: dqdchi, dqdmu
    integer:: i, j, iAtom, iMol, jAtom, jMol, iq, nq, k, l
    real(wp), dimension(3):: R

    dqdchi = 0.0_wp
    dqdmu = 0.0_wp

    ! compute derivatives of charges
    call  GetFQsDerivatives(nFQmol, nFQAtoms, nAtomtypes,  atomtypes, FQMolecules, FQJinv, V, dqdchi, dqdmu)

    ! compute dF/dchi and dF/dmu
    do i = 1, nAtomtypes
      ! Eint contribution
      do j = 1, nFQAtoms
        gradients(i) = gradients(i) + 2*wEint*(EintFQ-EintQM)*dqdchi(j, i)*V(j)
        gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wEint*(EintFQ-EintQM)*dqdmu(j, i)*V(j)
      enddo

      ! dipole contribution
      j = 0
      do jMol = 1, nFQMol
        do jAtom = 1, FQMolecules(jMol)%nAtomsInMol
          j = j+1
          R = FQMolecules(jMol)%Coords(:,jAtom) - DipoleOrigin
          do k = 1, 3
            gradients(i) = gradients(i) + 2*wdip*(dipFQ(k) - QMdip(k))*dqdchi(j, i)*R(k)
            gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wdip*(dipFQ(k) - QMdip(k))*dqdmu(j, i)*R(k)
          enddo
        enddo
      enddo

      ! dipole contribution
      j = 0
      do jMol = 1, nFQMol
        do jAtom = 1, FQMolecules(jMol)%nAtomsInMol
          j = j+1
          R = FQMolecules(jMol)%Coords(:,jAtom) - QuadOrigin
          do k = 1, 3
            do l = k, 3
              gradients(i) = gradients(i) + 2*wquad*(quadFQ(k+l) - QMquad(k+l))*dqdchi(j, i)*R(k)*R(l)
              gradients(nAtomTypes+i) = gradients(nAtomTypes+i) + 2*wquad*(quadFQ(k+l) - QMquad(k+l))*dqdmu(j, i)*R(k)*R(l)
            enddo  
          enddo
        enddo
      enddo
    enddo

  
  End Subroutine UpdateGradients 

End Module  fq_util 
