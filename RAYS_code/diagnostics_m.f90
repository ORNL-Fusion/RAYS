 module diagnostics_m
!   contains diagnostic switches, error flags, and a routine, message(), to conveniently
!   output messages and numbers
!   The fortran I/O unit number -> 'message_unit' must be set in initialization.
!   A switch -> 'verbosity' gives central control of how much output is produced.
!   since messages are skipped unless the argument 'threshold' is >= 'verbosity'.
!   Routines are also provided to toggle 'message_unit' the originally specified
!   unit number over to std I/O (i.e. unit=6) and back.
!   'set_message_unit_stdio' and 'reset_message_unit_file'

! Specifications

    implicit none



! generic procedure: message(caller/character, mess/character, value/generic,
!  threshold/integer)
!    
! Prints "caller: mess= value" when threshold > verbosity
! value can be integer,real or complex

    interface message
    module procedure blank_message, &
    &   i_message, logical_message, r_message, rdbl_message, c_message, &
    &   ivn_message,  rvn_message, rvndbl_message, cvn_message, &
    &   cvndbl_message, rmatrix_message, rmatrixdbl_message, cmatrix_message, &
    &   cmatrixdbl_message
    end interface

    interface text_message
    module procedure text_message, two_texts_message
    end interface
    
! Time and date vector - loaded in subroutine initialize()
    integer :: date_v(8)

! message_unit = device number where output from message() goes
    integer :: message_unit = 11
    integer :: save_message_unit

!  File name for message output
    character(len=80) :: message_file = 'messages'

! message_unit = for standard I/O for writting to screen
    integer, parameter :: message_unit_stdio = 6

!  verbosity = a switch to set the level of output from message()
!  verbosity = 0 gives minimum output 
    integer :: verbosity

!  Run label (N.B. should be legal in file name, e.g.no blanks allowed)
    character(len=15) :: run_label = ''

!  Run description
    character(len=80) :: run_description = ''

!   Switch for diagnostics.
!   integrate_eq_gradients = false: default
!   integrate_eq_gradients = true: integrate gradients of B, Te, and ne along the ray.
    logical :: integrate_eq_gradients = .false.

!   Switch for numerically computing the derivatives of D.
!   num_deriv = 0: default
!   num_deriv = 1: do numerical computation of the derivatives of D.
    integer :: num_deriv

!   Error returns
    character(len=20) :: equib_err = ''
    character(len=20) :: ray_stop_flag = ''

!   Abort ode integration on error condition
    logical :: stop_ode = .false.
    
!   Counter for general use, mostly for use in debugging
    integer :: diag_count = 0
    integer :: max_diag_count = 5
    
!   Timing variables
    real :: t_start_rays, t_finish_rays, t_start_tracing, t_finish_tracing    

    namelist /diagnostics_list/ message_unit, verbosity, run_description, &
           & run_label, integrate_eq_gradients

!******************************    
    
contains

!******************************
  subroutine initialize_diagnostics

    use constants_m, only : input_unit
    
    implicit none

    call cpu_time(t_start_rays)
    
! Read input namelist
    open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
    read(input_unit, diagnostics_list)
    close(unit=input_unit)

!   Open file for output messages 
    open(unit=message_unit, file=trim(message_file),                    &
       & action='write', status='replace', form='formatted')     

! Write input namelist
    write(message_unit, diagnostics_list)

    return
    
  end subroutine initialize_diagnostics


! Print blank line  i.e. skip a line 
    subroutine blank_message()

    implicit none

    write (message_unit, '(" ")')

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
           else             
             write (message_unit, '(a)') text
           end if
        end if

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
           else             
             write (message_unit, '(a," ", a)') text1, text2
           end if
        end if

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
           else             
                if (abs(value) < 1.e4) then
                    write (message_unit, '(a, " = ", I6)') &
                &    mess, value
                else
                    write (message_unit, '(a, " = ", 1pe11.4)') &
                &    mess, real(value)
                end if
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
           else             
                if ( value ) then
                    write (message_unit, '(a, " = .true. ")') &
                    & mess
                else
                    write (message_unit, '(a, " = .false.")') &
                    & mess
            end if
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
       else             
            if ( (abs(value) <  1.e4) .and. (abs(value) > 1.e-3) ) then
                write (message_unit, '(a, " = ", f11.4)') &
                & mess, value
            else
                write (message_unit, '(a, " = ", 1pe11.4)') &
                & mess, value 
            end if
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
       else             
            if ( (abs(value) <  1.e4) .and. (abs(value) > 1.e-3) ) then
                write (message_unit, '(a, " = (", f11.4,", ",f11.4,")")')&
                & mess, value
            else
                write (message_unit, '(a, " = (", 1pe11.4,", ",1pe11.4,")")')&
                & mess, value   
            end if
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
       else             
            v_max = 0
    
            do i = 1, length_n
                if (abs(value(i)) > v_max ) v_max = abs(value(i))
            end do

            if ((v_max < 1.e5)) then
                write (message_unit, '(a, " = ")') mess
                write (message_unit, '(10i8)') (value(i), i=1, length_n)

            else
                write (message_unit, '(a, " = ")') mess
                write (message_unit, '(10(1pe12.4))') (value(i), i=1, length_n)
            end if
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
       else             
            v_min = huge(v_min)
            v_max = tiny(v_max)
            do i = 1, length_n
                if ((abs(value(i)) > 0.).and.(abs(value(i)) < v_min)) v_min = abs(value(i))
                if (abs(value(i)) > v_max ) v_max = abs(value(i))
            end do

            if (length_n <= 10 ) then
    
                if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                    write (message_unit, '(a, " = ", 10f12.4,")")') &
                    & mess, (value(i), i=1, length_n)
                else
                    write (message_unit, '(a, " = ", 10(1pe12.4),")")') &
                    & mess, (value(i), i=1, length_n) 
                end if
            else
    
                if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                    write (message_unit, '(10f12.4)') (value(i), i=1, length_n)
                else
                    write (message_unit, '(a, " = ")') mess
                    write (message_unit, '(10(1pe12.4))') (value(i), i=1, length_n) 
                end if
        
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
       else             
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
                    write (message_unit, '(a, " = ", 5(2f12.5,5x),")")') &
                    & mess, (value(i), i=1, length_n)
                else
                    write (message_unit, '(a, " = ", 5(1pe12.4,1x,1pe12.4,4x),")")') &
                    & mess, (value(i), i=1, length_n) 
                end if
    
            else
    
                if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                    write (message_unit, '(a, " = ")') mess
                    write (message_unit, '(5(2f10.4,5x))') (value(i), i=1, length_n)
                else
                    write (message_unit, '(a, " = ")') mess
                    write (message_unit, '(5(2(1pe11.4),5x))') (value(i), i=1, length_n) 
                end if
    
            end if
!******************************
       end if
     end if
    
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
       else             
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
                    write (message_unit, '(a, " = ", 5(2f12.5,5x),")")') &
                    & mess, (value(i), i=1, length_n)
                else
                    write (message_unit, '(a, " = ", 5(1pe12.4,1x,1pe12.4,4x),")")') &
                    & mess, (value(i), i=1, length_n) 
                end if
    
            else
    
                if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                    write (message_unit, '(a, " = ")') mess
                    write (message_unit, '(5(2f10.4,5x))') (value(i), i=1, length_n)
                else
                    write (message_unit, '(a, " = ")') mess
                    write (message_unit, '(5(2(1pe11.4),5x))') (value(i), i=1, length_n) 
                end if
    
            end if
!******************************
       end if
    end if
    
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
       else             
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
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(10f12.5)' )  (value(i, j), j=1, n_dim)
                end do
            
            else
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(10(1pe12.4))' )  (value(i, j), j=1, n_dim)
                end do
            
            end if
!******************************
       end if
    end if
    

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
       else             
!******************************
            if ( (abs(value) <  1.e4) .and. (abs(value) > 1.e-3) ) then
                write (message_unit, '(a, " = ", f11.4)') &
                & mess, value
            else
                write (message_unit, '(a, " = ", 1pe11.4)') &
                & mess, value 
            end if
!******************************
       end if
    end if

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
       else             
!******************************
            v_min = huge(v_min)
            v_max = tiny(v_max)
    
            do i = 1, length_n
                if ((abs(value(i)) > 0.).and.(abs(value(i)) < v_min)) v_min = abs(value(i))
                if (abs(value(i)) > v_max ) v_max = abs(value(i))
            end do

            if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
                write (message_unit, '(a, " = ", 10f12.4,")")') &
                & mess, (value(i), i=1, length_n)
            else
                write (message_unit, '(a, " = ", 10(1pe12.4),")")') &
                & mess, (value(i), i=1, length_n) 
            end if
!******************************
       end if
    end if
    

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
       else             
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
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(10f12.5)' )  (value(i, j), j=1, n_dim)
                end do
            
            else
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(10(1pe12.4))' )  (value(i, j), j=1, n_dim)
                end do
            
            end if
!******************************
       end if
    end if
    

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
       else             
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
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(i3,8(3x,2f10.4))' )  i, (value(i, j), j=1, n_dim)
                end do
            
            else
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(i3,8(2x,2(1pe12.4)))' )  i, (value(i, j), j=1, n_dim)
                end do
            
            end if
!******************************
       end if
    end if
    
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
       else             
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
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(i3,2x,8(3x,2f10.4))' )  i, (value(i, j), j=1, n_dim)
                end do
            
            else
    
                write (message_unit, '(a, " = ")' ) mess
                do i = 1, m_dim
                    write (message_unit, '(i3,2x,8(3x,2(1pe11.4)))' )  i, (value(i, j), j=1, n_dim)
                end do
            
            end if
!******************************
       end if
    end if
    
    return
    end subroutine cmatrixdbl_message 
    
    
    subroutine set_message_unit_stdio
    implicit none
        save_message_unit = message_unit
        message_unit = message_unit_stdio
    return
    end subroutine set_message_unit_stdio



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
            write (message_unit, *) 'diagnostic_counter issued stop: count = ',&
               &diag_count, '  = max_count'
            stop
        end if
    end subroutine diagnostic_counter

 end module diagnostics_m
