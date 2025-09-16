module template_m
! This is a template for formating the header of a module file.  The header format is
! intended to allow a parser script to automatically generate documentation files.

! The header consists of:
! 1) the top-level comments
! 2) Working notes comments
! 3) declaration of module data (not including variables that are input from namelists)
!    definitions of derived types and interfaces
! 4) declaration of variables that come in from namelists.
!
! So the header is the material above the "contains" line.  Each of these blocks is
! considered to be a section of the header. Each section is delineated by marker lines
! of the form "!_____ ...".  These marker lines are signals to the parser. Except for the
! top-level comments, which follow immediately after the "module <module name>" line,
! the beginning, of a section is indicated by a 3 line boundary marker of the form:

!_________________________________________________________________________________________
! <section name>
!_________________________________________________________________________________________

N.B. If you want to add other separators or decorations please use something else, such as
"!************...".

As an example, the working notes section is demarked by:
!_________________________________________________________________________________________
! Working notes
!_________________________________________________________________________________________

There sould be at least one marker line below the material in the header. (e.g. a contains
line)

Top-level comments: should describe the module, tell what the module does, and give other
useful info.

Working notes: should serve as a change log describing major additions and changes to
implementation. Or could be notes to self on issues or plans for more stuff.  If the
changes add functionality or affect the interfaces they should also appear in the
top-level comments section at top. Working notes are not included in the automatically
generated documentation, so you can be as informal as you like. Example:

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

! (9/4/2025 DBB) Added some signal lines for parser "!_____..."

Module data: Here put module-level "use" statements,  declaration of module variables,
and other entities such as derived type definitions and interface blocks. Be sure to
comment the definition of each entity.  Example:

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use ABC_m , only :
    use XYZ_m , only :

    implicit none

! number of widgets used
	integer : n_widgets

! zork level as defined in Plan 9 from outer space
	real(KIND=rkind) :: x_zork

! These may not need a comment.  Be real.
	integer :: i, j, n

! Derived type containing data galaxy definition
    type galaxy

    !   Number of stars
        integer :: n_stars

	!   Standard name of star
        character(len = 256), allocatable :: star_ID(:)

	!   Average density of dark matter
		real(KIND=rkind) :: dens_dark_matter

    !   Error returns
        character(len=60) :: err = ''

    end type galaxy

Namelist data: Put declaration of namelist variables at end of declarations, i.e. just
before the namelist itself. Include the name of the namelist in the marker line.  There
can be zero or more than one Namelist data section.  Example:

!_________________________________________________________________________________________
! Namelist data for /template_list/
!_________________________________________________________________________________________

! Number of galaxies
    integer :: n_galaxies

! Name of fusion model
	character(len = 80) :: fusion_model

 namelist /template_list/&
     & n_galaxies, &
     & fusion_model

The "!_________ ..." line below, i.e. above "contains" indicates end of the header

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

Below is typical shape of a RAYS module.  As of now 9/2025 the parser doesn''t look
below here

  subroutine initialize_template_m

	  Whatever it takes

    return
  end subroutine initialize_template_m

Below is a separator line that won''t  comfuse the parser
!********************************************************************

  subroutine calculate_impossible(a,b,c)

	Clever calculation with lots of comments

    return
 end subroutine calculate_impossible

!********************************************************************

    subroutine deallocate_template_m

	Clean up mess

    end subroutine deallocate_template_m

!********************************************************************

end module template_m
