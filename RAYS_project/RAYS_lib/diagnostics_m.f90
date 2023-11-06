 module diagnostics_m
!   Controls output to files, stdio, etc
!
!   Contains and initializes fortran output unit numbers.
!   Contains diagnostic switches, error flags, and other miscellaneous data
!   Provides two generic routines (message and text_message) to produce formatted output
!   requiring only one line of code.
!   A switch -> 'verbosity' gives central control of how much output is produced.
!   messages are suppressed unless the optional argument 'threshold' is >= 'verbosity'.
!   Routines are also provided to toggle 'message_unit' the originally specified
!   unit number over to a diffferent unit and back so messages could be sent to other files
!   'set_message_unit_stdio' and 'reset_message_unit_file'

!   External procedures: cpu_time (intrinsic)

! Specifications

    implicit none

! generic procedure: message(mess/character, value/generic, threshold/integer)
! Prints "caller: mess= value" when threshold > verbosity
! value can be integer,real, complex, or matrix
    interface message
    module procedure blank_message, &
    &   i_message, logical_message, r_message, rdbl_message, c_message, &
    &   ivn_message,  rvn_message, rvndbl_message, cvn_message, &
    &   cvndbl_message, rmatrix_message, rmatrixdbl_message, cmatrix_message, &
    &   cmatrixdbl_message
    end interface

! generic procedure: text_message(mess/character, threshold/integer)
!                 or text_message(mess1/character,mess2/character, threshold/integer)

    interface text_message
    module procedure text_message, two_texts_message
    end interface

! Time and date vector - loaded in subroutine initialize()
    integer :: date_v(8)

! Unit where output from message() goes
    integer :: message_unit

! Place to temporarily stash message_unit to send messages to alternate file
    integer :: save_message_unit

! Unit for alternate message file
    integer :: message_unit_alternate ! Superseded, stdout -> write(*...

!  Default file name for message output.  At end copied to "log.RAYS.<run label>"
    character(len=80) :: message_file = 'messages'

!  verbosity = a switch to set the level of output from message() and text_message()
!  verbosity = 0 gives minimum output
    integer :: verbosity

!   messages_to_stdout = If true write messages to stdout as well as to messages file
    logical :: messages_to_stdout = .false.

!   write_formatted_ray_files = If true write ray data to formatted ASCII file
    logical :: write_formatted_ray_files = .false.

!   unit numbers for persistent files (i.e. ones that stay open all through the run)
    integer :: output_unit   ! Unit for formatted ray data output.  Set in initialize()
    integer :: ray_list_unit ! Unit for formatted ray list output.  Set in initialize()

!  Run label (N.B. should be legal in file name, e.g.no blanks allowed)
    character(len=60) :: run_label = ''

!  Run description
    character(len=80) :: run_description = ''

!   Switch for diagnostics.
!   integrate_eq_gradients = false: default
!   integrate_eq_gradients = true: integrate gradients of B, Te, and ne along the ray.
    logical :: integrate_eq_gradients = .false.

!   Counter for general use, mostly for use in debugging
    integer :: diag_count = 0
    integer :: max_diag_count = 5

!   Timing variables
    real :: t_start_rays, t_finish_rays, t_start_tracing, t_finish_tracing

    namelist /diagnostics_list/ verbosity, messages_to_stdout, write_formatted_ray_files, &
           & run_description, run_label, integrate_eq_gradients

!******************************

contains

!******************************
  subroutine initialize_diagnostics(read_input)

    implicit none
    logical, intent(in) :: read_input

    integer :: n_args ! number of command line arguments
    character(len=80) :: namelist_file = ''
	integer :: input_unit, get_unit_number ! External, free unit finder

    call cpu_time(t_start_rays)

! Default filename is 'rays.in'.  Optionally get input file name from command line. then
! copy that file to 'ray.in'
    n_args = command_argument_count()
    if(n_args > 1) then
        write(*,*) 'RAYS takes zero or one command line argument -> namelist filename'
        stop 'incorrect command line arguments'
    else if (n_args == 1) then
        call get_command_argument(1,namelist_file)
        call system('cp '//trim(namelist_file)//' rays.in')
    end if

	! Get unit number and open file for output messages
		message_unit = get_unit_number()
		open(unit=message_unit, file=trim(message_file),                    &
		   & action='write', status='replace', form='formatted')

    if (read_input .eqv. .true.) then
    ! Read input namelist
    	input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, diagnostics_list)
        close(unit=input_unit)
    end if

    call text_message('Initializing diagnostics', 1)

! Write input namelist
    if (verbosity > 0) then
		write(message_unit, diagnostics_list)
		if (messages_to_stdout) write(*, diagnostics_list)
		call message(1)
    end if

    return

  end subroutine initialize_diagnostics


! Print blank line  i.e. skip a line
    subroutine blank_message(threshold)
		implicit none
		integer, optional, intent (in) :: threshold

		if (present(threshold)) then
		   if (verbosity < threshold) then
			 return
		   end if
		end if

		write(message_unit, '(" ")')
		if (messages_to_stdout) write(*, '(" ")')

		return
    end subroutine blank_message



! Print text line if threshold greater than verbosity
    subroutine text_message(text, threshold)

        implicit none
        character (len=*), intent (in) :: text
        integer, optional, intent (in) :: threshold

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if

        write(message_unit, '(a)') trim(text)
        if (messages_to_stdout) write(*, '(a)') trim(text)

        return
    end subroutine text_message


! Print text line with 2 character inputs e.g. a string and a character variable
    subroutine two_texts_message(text1,text2, threshold)

        implicit none
        character (len=*), intent (in) :: text1, text2
        integer, optional, intent (in) :: threshold

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if

        write(message_unit, '(a," ", a)') trim(text1), trim(text2)
        if (messages_to_stdout) write(*, '(a," ", a)') trim(text1), trim(text2)

        return
    end subroutine two_texts_message


! Print message with integer value
    subroutine i_message (mess, value, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    integer, intent (in) :: value
    integer, optional, intent (in) :: threshold

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if

        if (abs(value) < 1.e4) then
            write(message_unit, '(a, " = ", I6)') &
            &    mess, value
            if (messages_to_stdout) then
				write(*, '(a, " = ", I6)') &
				&    mess, value
            end if

        else
            write(message_unit, '(a, " = ", 1pe11.4)') &
            &    mess, real(value)
            if (messages_to_stdout) then
				write(*, '(a, " = ", 1pe11.4)') &
				&    mess, real(value)
            end if

        end if

        return
    end subroutine i_message


! Print message with logical value
    subroutine logical_message (mess, value, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    logical, intent (in) :: value
    integer, optional, intent (in) :: threshold

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if

        if ( value ) then
            write(message_unit, '(a, " = .true. ")') &
            & mess
            if (messages_to_stdout) then
				write(*, '(a, " = .true. ")') &
				& mess
            end if
        else
            write(message_unit, '(a, " = .false.")') &
            & mess
            if (messages_to_stdout) then
				write(*, '(a, " = .false.")') &
				& mess
            end if
        end if
    return
    end subroutine logical_message


! Print message with default real value
    subroutine r_message (mess, value, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    real, intent (in) :: value
    integer, optional, intent (in) :: threshold

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if

        if ( (abs(value) <  1.e4) .and. (abs(value) > 1.e-3) ) then
            write(message_unit, '(a, " = ", f11.4)') &
            & mess, value
            if (messages_to_stdout) then
				write(*, '(a, " = ", f11.4)') &
				& mess, value
            end if
        else
            write(message_unit, '(a, " = ", 1pe11.4)') &
            & mess, value
            if (messages_to_stdout) then
				write(*, '(a, " = ", 1pe11.4)') &
				& mess, value
            end if
        end if
    return
    end subroutine r_message


! Print message with default complex value
    subroutine c_message (mess, value, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    integer, optional, intent (in) :: threshold
    complex, intent (in) :: value

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
        if ( (abs(value) <  1.e4) .and. (abs(value) > 1.e-3) ) then
            write(message_unit, '(a, " = (", f11.4,", ",f11.4,")")')&
            & mess, value
            if (messages_to_stdout) then
            write(*, '(a, " = (", f11.4,", ",f11.4,")")')&
            & mess, value
            end if
        else
            write(message_unit, '(a, " = (", 1pe11.4,", ",1pe11.4,")")')&
            & mess, value
            if (messages_to_stdout) then
				write(*, '(a, " = (", 1pe11.4,", ",1pe11.4,")")')&
				& mess, value
            end if
        end if
    return
    end subroutine c_message



!   print messsage with integer n-vector value
    subroutine ivn_message (mess, value, length_n, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    integer, dimension(*), intent (in) :: value
    integer, intent (in) :: length_n
    integer, optional, intent (in) :: threshold


    real :: v_max
    integer :: i

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
        v_max = 0

        do i = 1, length_n
            if (abs(value(i)) > v_max ) v_max = abs(value(i))
        end do

        if ((v_max < 1.e5)) then
            write(message_unit, '(a, " = ")') mess
            write(message_unit, '(10i8)') (value(i), i=1, length_n)
            if (messages_to_stdout) then
				write(*, '(a, " = ")') mess
				write(*, '(10i8)') (value(i), i=1, length_n)
            end if
        else
            write(message_unit, '(a, " = ")') mess
            write(message_unit, '(10(1pe12.4))') (value(i), i=1, length_n)
            if (messages_to_stdout) then
				write(*, '(a, " = ")') mess
				write(*, '(10(1pe12.4))') (value(i), i=1, length_n)
            end if
        end if
    return
    end subroutine ivn_message




!   print messsage with default real n-vector value
    subroutine rvn_message (mess, value, length_n, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    real, dimension(*), intent (in) :: value
    integer, intent (in) :: length_n
    integer, optional, intent (in) :: threshold

    real :: v_min, v_max
    integer :: i

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if

        v_min = huge(v_min)
        v_max = tiny(v_max)
        do i = 1, length_n
            if ((abs(value(i)) > 0.).and.(abs(value(i)) < v_min)) v_min = abs(value(i))
            if (abs(value(i)) > v_max ) v_max = abs(value(i))
        end do

        if (length_n <= 10 ) then

            if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                write(message_unit, '(a, " = ", 10f12.4,")")') &
                & mess, (value(i), i=1, length_n)
            if (messages_to_stdout) then
                write(*, '(a, " = ", 10f12.4,")")') &
                & mess, (value(i), i=1, length_n)
            end if
            else
                write(message_unit, '(a, " = ", 10(1pe12.4),")")') &
                & mess, (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ", 10(1pe12.4),")")') &
					& mess, (value(i), i=1, length_n)
				end if
            end if
        else

            if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                write(message_unit, '(10f12.4)') (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(10f12.4)') (value(i), i=1, length_n)
				end if
            else
                write(message_unit, '(a, " = ")') mess
                write(message_unit, '(10(1pe12.4))') (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ")') mess
					write(*, '(10(1pe12.4))') (value(i), i=1, length_n)
				end if
            end if

        end if
     return
     end subroutine rvn_message


!   print messsage with default complex n-vector value
    subroutine cvn_message (mess, value, length_n, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    complex, dimension(*), intent (in) :: value
    integer, intent (in) :: length_n
    integer, optional, intent (in) :: threshold

    real :: Re_min, Re_max, Im_min, Im_max, absRe, absIm, v_min, v_max
    integer :: i

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        Re_min = huge(v_min)
        Im_min = huge(v_min)
        Re_max = tiny(v_min)
        Im_max = tiny(v_min)

        do i = 1, length_n
            absRe = abs(real(value(i)))
            absIm = abs(aimag(value(i)))

            if ((absRe > 0.).and.(absRe < Re_min)) Re_min = absRe
            if ((absIm > 0.).and.(absIm < Im_min)) Im_min = absIm
            if (absRe > Re_max ) Re_max = absRe
            if (absIm > Im_max ) Im_max = absIm
        end do

        v_min = min(Re_min, Im_min)
        v_max = max(Re_max, Im_max)

        if (length_n <= 5 ) then

            if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                write(message_unit, '(a, " = ", 5(2f12.5,5x),")")') &
                & mess, (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ", 5(2f12.5,5x),")")') &
					& mess, (value(i), i=1, length_n)
				end if
				else
                write(message_unit, '(a, " = ", 5(1pe12.4,1x,1pe12.4,4x),")")') &
                & mess, (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ", 5(1pe12.4,1x,1pe12.4,4x),")")') &
					& mess, (value(i), i=1, length_n)
				end if
            end if

        else

            if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                write(message_unit, '(a, " = ")') mess
                write(message_unit, '(5(2f10.4,5x))') (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ")') mess
					write(*, '(5(2f10.4,5x))') (value(i), i=1, length_n)
				end if
           else
                write(message_unit, '(a, " = ")') mess
                write(message_unit, '(5(2(1pe11.4),5x))') (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ")') mess
					write(*, '(5(2(1pe11.4),5x))') (value(i), i=1, length_n)
				end if
            end if

        end if
!******************************
    return
    end subroutine cvn_message



!   print messsage with complex double precision n-vector value
    subroutine cvndbl_message (mess, value, length_n, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    complex(kind=kind(1.d0)), dimension(*), intent (in) :: value
    integer, intent (in) :: length_n
    integer, optional, intent (in) :: threshold

    real(kind=kind(1.d0)) :: Re_min, Re_max, Im_min, Im_max, absRe, absIm, &
        & v_min, v_max
    integer :: i

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        Re_min = huge(v_min)
        Im_min = huge(v_min)
        Re_max = tiny(v_min)
        Im_max = tiny(v_min)

        do i = 1, length_n
            absRe = abs(real(value(i)))
            absIm = abs(aimag(value(i)))

            if ((absRe > 0.).and.(absRe < Re_min)) Re_min = absRe
            if ((absIm > 0.).and.(absIm < Im_min)) Im_min = absIm
            if (absRe > Re_max ) Re_max = absRe
            if (absIm > Im_max ) Im_max = absIm
        end do

        v_min = min(Re_min, Im_min)
        v_max = max(Re_max, Im_max)

        if (length_n <= 5 ) then

            if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                write(message_unit, '(a, " = ", 5(2f12.5,5x),")")') &
                & mess, (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ", 5(2f12.5,5x),")")') &
					& mess, (value(i), i=1, length_n)
				end if
            else
                write(message_unit, '(a, " = ", 5(1pe12.4,1x,1pe12.4,4x),")")') &
                & mess, (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ", 5(1pe12.4,1x,1pe12.4,4x),")")') &
					& mess, (value(i), i=1, length_n)
				end if
            end if

        else

            if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                write(message_unit, '(a, " = ")') mess
                write(message_unit, '(5(2f10.4,5x))') (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ")') mess
					write(*, '(5(2f10.4,5x))') (value(i), i=1, length_n)
				end if
            else
                write(message_unit, '(a, " = ")') mess
                write(message_unit, '(5(2(1pe11.4),5x))') (value(i), i=1, length_n)
				if (messages_to_stdout) then
					write(*, '(a, " = ")') mess
					write(*, '(5(2(1pe11.4),5x))') (value(i), i=1, length_n)
				end if
            end if

        end if
!******************************
    return
    end subroutine cvndbl_message


!   print messsage with real mxn matrix value
    subroutine rmatrix_message (mess, value, m_dim, n_dim, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    real, intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: threshold

    real :: v_min, v_max
    integer :: i, j

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        v_min = huge(v_min)
        v_max = tiny(v_max)

        do i = 1, m_dim
            do j = 1, n_dim

                if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
                    & v_min = abs(value(i,j))
                if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
            end do
        end do


        if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then

            write(message_unit, '(a, " = ")' ) mess
			if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
			end if
            do i = 1, m_dim
                write(message_unit, '(10f12.5)' )  (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(10f12.5)' )  (value(i, j), j=1, n_dim)
				end if
            end do

        else

            write(message_unit, '(a, " = ")' ) mess
			if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
			end if
            do i = 1, m_dim
                write(message_unit, '(10(1pe12.4))' )  (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(10(1pe12.4))' )  (value(i, j), j=1, n_dim)
				end if
            end do

        end if
!******************************
    return
    end subroutine rmatrix_message



! Print message with real double precision value

    subroutine rdbl_message (mess, value, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    real(kind=kind(1.D0)), intent (in) :: value
    integer, optional, intent (in) :: threshold

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        if ( (abs(value) <  1.e4) .and. (abs(value) > 1.e-3) ) then
            write(message_unit, '(a, " = ", f11.4)') &
            & mess, value
            if (messages_to_stdout) then
				write(*, '(a, " = ", f11.4)') &
				& mess, value
            end if
        else
            write(message_unit, '(a, " = ", 1pe11.4)') &
            & mess, value
            if (messages_to_stdout) then
				write(*, '(a, " = ", 1pe11.4)') &
				& mess, value
            end if
        end if
!******************************
    return
    end subroutine rdbl_message


!   print messsage with real double precision n-vector value

    subroutine rvndbl_message (mess, value, length_n, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    real(kind=kind(1.D0)), dimension(*), intent (in) :: value
    integer, intent (in) :: length_n
    integer, optional, intent (in) :: threshold


    real(kind=kind(1.D0)) :: v_min, v_max
    integer :: i

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        v_min = huge(v_min)
        v_max = tiny(v_max)

        do i = 1, length_n
            if ((abs(value(i)) > 0.).and.(abs(value(i)) < v_min)) v_min = abs(value(i))
            if (abs(value(i)) > v_max ) v_max = abs(value(i))
        end do

        if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
            write(message_unit, '(a, " = ", 10f12.4,")")') &
            & mess, (value(i), i=1, length_n)
            if (messages_to_stdout) then
				write(*, '(a, " = ", 10f12.4,")")') &
				& mess, (value(i), i=1, length_n)
            end if
       else
            write(message_unit, '(a, " = ", 10(1pe12.4),")")') &
            & mess, (value(i), i=1, length_n)
            if (messages_to_stdout) then
				write(*, '(a, " = ", 10(1pe12.4),")")') &
				& mess, (value(i), i=1, length_n)
            end if
        end if
!******************************
    return
    end subroutine rvndbl_message



!   print messsage with real double precision mxn matrix value

    subroutine rmatrixdbl_message (mess, value, m_dim, n_dim, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    real(kind=kind(1.D0)), intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: threshold

    real(kind=kind(1.D0)) :: v_min, v_max
    integer :: i, j

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        v_min = huge(v_min)
        v_max = tiny(v_max)

        do i = 1, m_dim
            do j = 1, n_dim

                if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
                    & v_min = abs(value(i,j))
                if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
            end do
        end do


        if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then

            write(message_unit, '(a, " = ")' ) mess
            if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
            end if
            do i = 1, m_dim
                write(message_unit, '(10f12.5)' )  (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(10f12.5)' )  (value(i, j), j=1, n_dim)
				end if
            end do

        else

            write(message_unit, '(a, " = ")' ) mess
			if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
			end if
           do i = 1, m_dim
                write(message_unit, '(10(1pe12.4))' )  (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(10(1pe12.4))' )  (value(i, j), j=1, n_dim)
				end if
            end do

        end if
!******************************
    return
    end subroutine rmatrixdbl_message




!   print messsage with complex mxn matrix value

    subroutine cmatrix_message (mess, value, m_dim, n_dim, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    complex(kind=kind(1.E0)), intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: threshold


    real :: v_min, v_max
    integer :: i, j

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        v_min = huge(v_min)
        v_max = tiny(v_max)

        do i = 1, m_dim
            do j = 1, n_dim
                if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
                    & v_min = abs(value(i,j))
                if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
            end do
        end do


        if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then

            write(message_unit, '(a, " = ")' ) mess
			if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
			end if
            do i = 1, m_dim
                write(message_unit, '(i3,8(3x,2f10.4))' )  i, (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(i3,8(3x,2f10.4))' )  i, (value(i, j), j=1, n_dim)
				end if
            end do

        else

            write(message_unit, '(a, " = ")' ) mess
            if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
			end if
            do i = 1, m_dim
                write(message_unit, '(i3,8(2x,2(1pe12.4)))' )  i, (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(i3,8(2x,2(1pe12.4)))' )  i, (value(i, j), j=1, n_dim)
				end if
            end do

        end if
!******************************
    return
    end subroutine cmatrix_message


!   print messsage with complex double precision mxn matrix value

    subroutine cmatrixdbl_message (mess, value, m_dim, n_dim, threshold)

    implicit none
    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    complex(kind=kind(1.D0)), intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: threshold

    real(kind=kind(1.D0)) :: v_min, v_max
    integer :: i, j

        if (present(threshold)) then
           if (verbosity < threshold) then
             return
           end if
        end if
!******************************
        v_min = huge(v_min)
        v_max = tiny(v_max)

        do i = 1, m_dim
            do j = 1, n_dim

                if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
                    & v_min = abs(value(i,j))
                if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
            end do
        end do


        if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then

            write(message_unit, '(a, " = ")' ) mess
            if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
            end if
            do i = 1, m_dim
                write(message_unit, '(i3,2x,8(3x,2f10.4))' )  i, (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(i3,2x,8(3x,2f10.4))' )  i, (value(i, j), j=1, n_dim)
				end if
            end do

        else

            write(message_unit, '(a, " = ")' ) mess
            if (messages_to_stdout) then
				write(*, '(a, " = ")' ) mess
            end if
            do i = 1, m_dim
                write(message_unit, '(i3,2x,8(3x,2(1pe11.4)))' )  i, (value(i, j), j=1, n_dim)
				if (messages_to_stdout) then
					write(*, '(i3,2x,8(3x,2(1pe11.4)))' )  i, (value(i, j), j=1, n_dim)
				end if
            end do

        end if
!******************************
    return
    end subroutine cmatrixdbl_message


    subroutine set_message_unit_alternate
    implicit none
        save_message_unit = message_unit
        message_unit = message_unit_alternate
    return
    end subroutine set_message_unit_alternate



    subroutine set_message_unit_file(i, filename)
    implicit none

! With both 'i' and 'filename' present it opens file for message output on unit 'i'
! With only 'i' present it sets 'message_unit' to 'i'
! With no arguments it resets 'message_unit' to previously stored value

    integer, intent(in), optional :: i
    character*(*), intent(in), optional :: filename

    if ( present(i) .and. present(filename) ) then
        message_unit = i
    !   Open file for output messages.
        open(unit=message_unit, file=filename,                  &
       & action='write', status='replace', form='formatted')
    end if

    if  ( present(i) .and. (.not.present(filename)) ) message_unit = i

    if  ( (.not.present(i)) .and. (.not.present(filename)) )  &
                & message_unit = save_message_unit

    return
    end subroutine set_message_unit_file

! Counter for debugging purposes.  Count from zero to max_diag_count then stop.
    subroutine diagnostic_counter()
        diag_count = diag_count + 1
        if (diag_count >= max_diag_count) then
            write(message_unit, *) 'diagnostic_counter issued stop: count = ',&
               &diag_count, '  = max_count'
            stop
        end if
    end subroutine diagnostic_counter

!********************************************************************

    subroutine deallocate_diagnostics_m
        return ! Nothing to deallocate
    end subroutine deallocate_diagnostics_m

 end module diagnostics_m
