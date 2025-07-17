!This is file : optimizer
! Author = francescomazza
! Started at: 03.07.2025
! Last Modified: Fri Jul  4 18:26:29 2025
!
Program  optimizer

use fq_util
! use parsing

Implicit None

! provvisorio
real(wp), parameter:: Angstrom = 0.52917721090299996
integer:: nMol, nAtoms, nAtomTypes, nq, nQMcalc
integer, dimension(:), allocatable:: atomtypes, nAtomsInMol, MolID
real(wp), dimension(:), allocatable:: ChiParam, EtaParam, chi, eta, molCharges, V, qint, EintQM, gradients
real(wp), dimension(:,:), allocatable:: coordinates, qcoords, dqdchi, dqdmu, QMdip, QMQuad, dipint
character(len = 3), dimension(:), allocatable:: Symbols
real(wp):: EintFQ,  totalw, wEint, wdip, wquad, Functional, singleQ
real(wp), dimension(3)::  DipoleOrigin, QuadOrigin, dipFQ, Qcoord
real(wp), dimension(6):: QuadFQ

real(kind = wp), dimension(:), allocatable:: FQs 
real(kind = wp), dimension(:,:), allocatable:: FQJinv
type(FQMol), dimension(:), allocatable:: FQMolecules


integer:: i, iostat, it, maxit, iflag, iatom, imol, l

! Read coordinates, molecules and atom types from input file
nMol = 1
nAtoms = 3
nAtomTypes = 2
nQMcalc = 1

wEint = 10.0_wp
wdip = 2.0_wp
wquad = 1.0_wp 

maxit = 20000
iflag = 1


open(unit = 16, file='FQpar_data.dat' , iostat = iostat, action='read')
read(16, *) nQMcalc
read(16, '(a)') 

allocate(atomtypes(nAtomTypes), nAtomsInMol(nMol), MolID(nAtoms), ChiParam(nAtomTypes), EtaPAram(nAtomTypes), chi(nAtoms), &
  eta(natoms), coordinates(3, nAtoms), symbols(nAtoms), molcharges(nMol), FQs(nAtoms), V(natoms), qint(nQMcalc), & 
  qcoords(3, nQMcalc), dqdmu(nAtoms, natomtypes), dqdchi(nAtoms, natomtypes), EintQM(nQMcalc), QMdip(3, nQMcalc), &
  QMquad(6, nQMcalc), gradients(2*nAtomTypes), dipint(3, nQMcalc))

atomtypes = (/ 1, 1, 2 /)
nAtomsInMol = (/ 3 /)
MolCharges = (/ 0.0_wp /)
MolID =(/ 1, 1, 1 /)
ChiParam = (/ 0.0000001_wp, 0.11_wp/)
EtaParam = (/ 0.63_wp, 0.58_wp/)
coordinates(:, 1) = (/    0.75308062_wp,     0.60025412_wp,     0.00000000_wp /)
coordinates(:, 2) = (/   -0.75308062_wp,     0.60025412_wp,     0.00000000_wp /) 
coordinates(:, 3) = (/    0.00000000_wp,    -0.00025412_wp,     0.00000000_wp /)
symbols = (/ 'HW ', 'HW ', 'OW ' /)


DipoleOrigin = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
QuadOrigin = (/ 0.0_wp, 0.0_wp, 0.0_wp /)


! Read form file QM data
do i = 1, nQMcalc
  read(16, *) qcoords(:,i), qint(i), dipint(:,i), EintQM(i), QMdip(:,i), QMquad(:,i)
enddo

close(16)

QuadOrigin = (/ 0.0_wp, 0.12652035_wp, 0.0_wp /)


! convert coordinates to atomic units
coordinates = coordinates/Angstrom
qcoords = qcoords/Angstrom



! Build chi and eta list

do i = 1, nAtoms
  chi(i) = ChiParam(atomtypes(i))
  eta(i) = EtaParam(atomtypes(i))
enddo

FQs = 0.0_wp

call InitializeMolecules(nMol, nAtoms,  nAtomsInMol, MolID, MolCharges, chi, eta, coordinates, symbols, FQs, FQMolecules)

do i = 1, nmol
  print*, 'we are looking at ', i
  print*, 'natomsinmol', fqmolecules(i)%natomsinmol
  print*, 'symbols', fqmolecules(i)%symbols(:)
  print*, 'coords', fqmolecules(i)%coords(:,:)
  print*, 'molid', fqmolecules(i)%molid
  print*, 'chi', fqmolecules(i)%chi(:)
  print*, 'eta', fqmolecules(i)%eta(:)
  print*, 'charge', fqmolecules(i)%molcharge
end do

call MakeFQJ(nMol, nAtoms, FQMolecules, FQJinv)

call GetPotentialPointCharge(nMol, nAtoms, FQMolecules, V, qint(1), qcoords(:,1))

call UpdateFQs(nMol, nAtoms, FQMolecules, FQJinv, V)

do i = 1, nMol
  print*, 'We are looking at ', i
  print*, 'FQs', FQMolecules(i)%FQs(:)
enddo

EintFQ = intenergy(nMol, nAtoms, FQMolecules, V)
print*, 'E int: ', EintFQ

call  GetFQsDerivatives(nmol, nAtoms, nAtomtypes,  atomtypes, FQMolecules, FQJinv, V, dqdchi, dqdmu)

! update the Functional and the gradient by summing the contribution from all the charges
totalw = (wEint+wdip+wquad)
wEint = wEint/(totalw*nQMcalc)
wdip = wdip/(totalw*nQMcalc*3)
wquad = wquad/(totalw*nQMcalc*6)


do it = 1, maxit
  
  
  

  Functional = 0.0_wp
  gradients = 0.0_wp
  
  do i = 1, nQMcalc
    singleQ = qint(i) 
    Qcoord = Qcoords(:,i)
    call GetPotentialPointCharge(nMol, nAtoms, FQMolecules, V, SingleQ, Qcoord)
    call UpdateFQs(nMol, nAtoms, FQMolecules, FQJinv, V)
    EintFQ = intenergy(nMol, nAtoms, FQMolecules, V)
    dipFQ = Dipole(nMol, nAtoms, FQMolecules, DipoleOrigin)
    quadFQ = Quadrupole(nMol, nAtoms, FQMolecules, QuadOrigin)
  
    Functional = Functional+wEint*(EintFQ-EintQM(i))**2 
    Functional = Functional+wdip*norm2(DipFQ-QMDip(:,i))
    Functional = Functional+wquad*norm2(quadFQ-QMQuad(:,i))
  
    call  UpdateGradients(nMol, nAtoms, nAtomTypes, atomtypes,  FQMolecules, V, FQJinv, wEint, EintFQ, EintQM(i), &
      wdip, dipFQ, QMdip(:,i), wquad, quadFQ, QMquad(:,i), DipoleOrigin, QuadOrigin, gradients)
  
  
  enddo
  
  print*, Functional
  print*, gradients

  call GD(natomtypes, ChiParam, EtaParam, gradients, 1.0D-5, iflag)
  if (iflag .eq. 0) then
    print*, 'Convergence reached', it
    exit
  endif

  ! update parameters
  do i = 1, nAtoms
    chi(i) = ChiParam(atomtypes(i))
    eta(i) = EtaParam(atomtypes(i))
  enddo
  print*, 'chieta', chi, eta

  l = 0
  do iMol = 1, nMol
    do iAtom = 1, FQMolecules(iMol)%nAtomsInMol
      l = l+1
      FQMolecules(IMol)%chi(iAtom) = chi(l)
      FQMolecules(IMol)%eta(iAtom) = eta(l)
    enddo
  enddo




enddo
  
do i = 1, nmol
  print*, 'we are looking at ', i
  print*, 'natomsinmol', fqmolecules(i)%natomsinmol
  print*, 'symbols', fqmolecules(i)%symbols(:)
  print*, 'coords', fqmolecules(i)%coords(:,:)
  print*, 'molid', fqmolecules(i)%molid
  print*, 'chi', fqmolecules(i)%chi(:)
  print*, 'eta', fqmolecules(i)%eta(:)
  print*, 'charge', fqmolecules(i)%molcharge
end do


End Program  optimizer 
