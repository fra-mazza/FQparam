! This module contains useful subroutine for reading command line arguments 
! and the input file.

module parsing
use omp_lib
use iso_fortran_env, only : real64, output_unit
Implicit None

contains
  ! This subroutine reads command line arguments:
  ! -h--> display help message
  ! -i--> name of the input file. If not present, input.inp is chosen
  ! -o--> name of the output file. if not present, the name is built taking the radix of the input file and appending '.log' suffix
  ! -t--> name of the trajectory file. if not present, the name is built taking the radix of the input file and appending '_traj.log' suffix
  Subroutine  read_command_line_argument(onlyinput)
    use defaults, only : input_file, output_file, traj_file, FULLOUTPUT, output
    use misc, only: print_help
    implicit none
    logical, intent(in)  :: onlyinput
    integer    :: i
    character(256)    ::  commands, radix
    logical:: exist

    call get_command(commands)


    ! Search '-h' option and print help message
    i = -1
    i = index(commands, '-h') 
    if (i .gt. 0) then
      call print_help()
      stop
    endif

    ! Set only input and output file, trajectory file is optional and is set after having read the input
    if (onlyinput) then
      ! Search '-i' option and set input file name
      i = -1
      i = index(commands, '-i') 
      if (i .eq. 0) then
        write(unit = output_unit, fmt='(a)') 'No input file given. The software will search for input.inp.'
        write(unit = output_unit, fmt='(a, /)') "Please consider using the '-h' option to display the help message..."
      endif
        
      ! Read input file
      if (i .gt. 0) read(commands(i+2:), *) input_file



      write(unit = output_unit, fmt='(a, a)') 'Input file: ', trim(input_file)

      inquire(file = trim(input_file), exist = exist)
      if (.not. exist) then
        write(unit = output_unit, fmt='(a)') 'Input file does not exist.' 
        error stop
      endif


      ! Construct output and trajectory file name, to be used if-o and-t not present
      i = -1
      i = index(input_file, '.', back=.true.)
      radix = input_file(:i-1)
      output_file = trim(radix)//'.log'
      traj_file = trim(radix)//'_traj.xyz'

      ! Search '-o' option and set output file name
      i = -1
      i = index(commands, '-o') 
      if (i .eq. 0) write(unit = output_unit, fmt='(a, a, a)') 'Error while reading output file from command line: ',&
                                                             trim(output_file), '  will be used.'
      if (i .gt. 0 ) read(commands(i+2:), *) output_file
      write(unit = output_unit, fmt='(a, a)') 'Output file: ', trim(output_file)
    endif

      
    ! If required, search '-t'  option and set trajectory file name
    if (.not. onlyinput) then    
      if (output .eq. FULLOUTPUT) then
        i = -1
        i = index(commands, '-t') 
        if (i .eq. 0) write(unit = output_unit, fmt='(a, a, a)') 'Error while reading trajectory file from command line: ' &
                                                              , trim(traj_file), ' will be used.'
        if (i .gt. 0) read(commands(i+2:), *) traj_file
        write(unit = output_unit, fmt='(a, a)') 'Trajectory file: ', trim(traj_file)
      endif
    endif
  
  End Subroutine read_command_line_argument 

  ! This subroutine read the input file and recover the parameters of the calculation
  Subroutine  read_input()
    use defaults
    use misc 
    implicit none
    character(256)    :: value, line
    character(5)      :: key
    integer           :: i, iostat, j

    ! Already stopped at command line options reading
    !inquire(file = input_file, exist = exists)
    !if (.not. exists) then
    !  write(unit = iout, fmt='(a)')  'Input file does not exists'
    !  flush(iout)
    !  error stop
    !endif
    
    ! Open input file
    open(unit = iinp, file = input_file, iostat = iostat, action='read')
    if (iostat .ne. 0) write(unit = iout, fmt='(a)')  'Error while opening input file.'
    
    ! Read each line of the input, ignore left blank characters, skip empty and comment lines, and convert to uppercase to be case
    ! insensitive
    do 
      read(iinp, '(A)', End = 101) line
      line = adjustl(line)
      if (line(1:1) .eq. ' ') cycle
      if (line(1:1) .eq. '#') cycle
      ! capitalize line and consider only the first 5 non blank character
      key = to_upper(line(1:5))

      ! Select the corresponding keyword and initialize the parameters. Then write a copy of the input parameters in the output file
      ! If errors while reading, stop the program.
      select case (key)
    
        case ('CUTOF')
          lcutoff = .true.
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  "Error in CUTOFF input. Missing '=' character."
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) cutoff
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading CUTOFF input.'
            flush(iout)
            error stop
          endif
          cutoff = cutoff*1.0E-10_real64
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, f0.3)')  'CUTOFf', cutoff*1.0E10_real64
    
        case ('NUMIT')
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  "Error in NUMITERATIONS input. Missing '=' character."
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) num_it
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading NUMITERATIONS input. Missing ''='' character.'
            flush(iout)
            error stop
          endif
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, i0)')  'NUMITerations', num_it
    
        case ('COORD')
          ! Read the number of atoms
          read(iinp, *, iostat = iostat) n
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading COORDinates and velocities input.'
            flush(iout)
            error stop
          endif
          ! Read the empty line after n
          read(iinp, *) 
          allocate(atoms(n), x(3, n), v(3, n), f1(3, n), f2(3, n))  ! altre allocazioni da fare!!
          write(unit = iout, fmt='(2x, a)') 'COORDinates and velocities'
          write(unit = iout, fmt='(2x, i0, /)') n

          ! Read symbol, coordinates and velocities for each atom
          do i = 1, n
            read(iinp, '(A)') line
            read(line, *, iostat = iostat) atoms(i)%symbol, x(1, i), x(2, i), x(3, i), v(1, i), v(2, i), v(3, i)
            if (iostat .ne. 0) then
              write(unit = iout, fmt='(a)')  'Error while reading COORDinates and velocities input.'
              flush(iout)
              error stop
            endif
            do j = 1, 118
              if (atomic_symbols(j) .eq. atoms(i)%symbol) exit
            enddo
            atoms(i)%atnum = j
            atoms(i)%mass = atomic_masses(j)*1.66054E-27_real64  ! from u.m.a. to Kg
            atoms(i)%eps = atomic_eps(j)*4.359748199E-18_real64  ! fram Hartree to J
            atoms(i)%sig = atomic_sig(j)*1.0E-10_real64  ! from Angstrom to m

            !convert angstrom to m
            x(:,i) = x(:,i)*1.0E-10_real64

            write(unit = iout, fmt = fmt_xyz)  atoms(i)%symbol, x(:,i)*1.0E10_real64, v(:,i)
          end do

        case ('NOTRA') 
          output = ONLYLOG
          write(unit = iout, fmt='(2x, a)')  'NOTRAjectory'

        case ('LENBO')
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  'Error in LENBOX input. Missing ''='' character.'
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) l_box
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading LENBOX input.'
            flush(iout)
            error stop
          endif
          ! l_box is the half of the size of the box
          l_box = l_box/2.0_real64
          l_box = l_box*1.0E-10_real64
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, f0.3)')  'LENBOx', l_box*2.0E10_real64

        case ('TIMES') 
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  'Error in TIMESTEP input. Missing ''='' character.'
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) dt
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading TIMESTEP input.'
            flush(iout)
            error stop
          endif
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, es0.2)')  'TIMEStep', dt

        case ('NPROC') 
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  'Error in NPROC input. Missing ''='' character.'
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) nproc 
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading NPROCESSORS input.'
            flush(iout)
            error stop
          endif
          if (nproc .gt. omp_get_max_threads()) then
            write(unit = iout, fmt='(a, i0, a)')  'Error: NPROC larger than maximum thread available (',omp_get_max_threads(), ')'
            flush(iout)
            error stop
          endif
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, i0)')  'NPROCessors', nproc
          call omp_set_num_threads(nproc)

        case ('NPRIN')
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  'Error in NPRINT input. Missing ''='' character.'
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) nprint
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading NPRINT input.'
            flush(iout)
            error stop
          endif
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, i0)')  'NPRINT', nprint



        case ('INTEG')
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  'Error in integrator input. Missing ''='' character.'
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          if (to_upper(value(1:4)) .eq. 'VVER') then
            integrator = VVERLET
            write(unit = iout, fmt='(2x, a, t30, ''='', 2x, a)')  'INTEGRATOR', 'VVerlet'
          else if (to_upper(value(1:4)) .eq. 'LEAP') then
            integrator = LEAPFROG
            write(unit = iout, fmt='(2x, a, t30, ''='', 2x, a)')  'INTEGRATOR', 'Leap Frog'
          else
            write(unit = iout, fmt='(a)')  'Error while reading INTEGRATOR input.'
            flush(iout)
            error stop
          endif


        case ('THERM')
          ltermostat = .true.
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  'Error in THERMOSTAT input. Missing ''='' character.'
            flush(iout)
            error stop
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) temp_target
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading THERMOSTAT input.'
            flush(iout)
            error stop
          endif
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, f0.2)')  'THERMostat temperature', temp_target

        case ('AVERA')
          laverage = .true.
          i = -1
          i = index(line, '=')
          if (i .eq. 0) then
            write(unit = iout, fmt='(a)')  'No starting point for average. Using default: 80%.'
            cycle
          endif
          value = adjustl(line(i+1:))
          read(value, *, iostat = iostat) start_average
          if (iostat .ne. 0) then
            write(unit = iout, fmt='(a)')  'Error while reading AVERAGE START input.'
            flush(iout)
            error stop
          endif
          if (start_average .gt. 100.0_real64 .or. start_average .lt. 0.0_real64) then
            write(unit = iout, fmt='(a)')  'The average starting point should be a real number between 0 and 100'
            flush(iout)
            error stop
          endif
          write(unit = iout, fmt='(2x, a, t30, ''='', 2x, e0.2)')  'AVERAge start', start_average


        case ('NOAVE')
          laverage = .false.
          write(unit = iout, fmt='(2x, a)')  'NOAVErage'


     
        case DEFAULT
            write(unit = iout, fmt='(a, a)')  'Error in INPUT line: ', line
            flush(iout)
            error stop 
    
      end select
    
    
    
    end do 
    
    
    101 write(unit = iout, fmt='(a)')  'Input file read.'

    close(iinp)
  
  End Subroutine read_input 

  ! This subroutine check if some parameters are incompatible
  ! First it checks if the COORD keyword is present, then check that the box contains all the atoms.
  ! Finally it suggests not to use parallelization if n is small and suggests to choose an appropriate timestep
  Subroutine  check_warning_input()
    use defaults
    real (kind = real64)    :: x_max, x_min, abs_max

    ! Check that COORD is present
    if (.not. allocated(x) .or. n .eq. 0) then
      write(unit = iout, fmt='(/,a)')  'WARNING: Coordinates and velocities not provided!' 
      write(unit = iout, fmt='(a, /)') 'Please provide them in the input file using the keyword COORD'
      flush(iout)
      error stop
    endif


    ! Check that the box is large enought
    x_max = abs(maxval(x))
    x_min = abs(minval(x))

    abs_max = max(x_max, x_min)

    if (abs_max .gt. l_box) then
      write(unit = iout, fmt='(/,a)') 'WARNING: some atoms are outside the simulation box!'
      l_box = abs_max*1.2_real64
      write(unit = iout, fmt='(a, f0.4, a, /)')   'Changing the box lenght to LENBOX=', l_box*2.0E10_real64, ' Å'
    endif

    ! Check if the timestep provided is large
    if (dt .gt. 1.0E-14_real64) write(unit = iout, fmt='(/,a, /)') 'WARNING: Timestep is large. Consider choosing a &
                                                                 &smaller timestep to obtain accurate results.'

    ! Suggest to use serial code if n small
    if (nproc .gt. 1 .and. n .lt. 100) then
      write(unit = iout, fmt='(/,a, /)') 'WARNING: with a small number of atoms, it could be possible that a PARALLEL &
                                 &simulation is more expensive than a SERIAL simulation. Consider setting NPROC = 1'
    endif 
  
  End Subroutine check_warning_input 

End module parsing


