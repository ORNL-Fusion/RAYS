subroutine ode ( f, neqn, y, t, tout, relerr, abserr, iflag, work, iwork, ray_stop)
! A fortran 90 implementation of the vemnerable Shampine and Gordon ODE package, 
! subsequently modified to support ray (integration) termination conditions from the
! derivative subroutine f()

! Working notes:
! (DBB 2/4/2022) To allow thread safety removed 'stop_ode' from module ode_m and made
! a derived type 'ode_stop', made 'ray_stop' and argument of all subroutines.
! (DBB 12/2021) Converted to use KIND = rkind for all reals.
! (DBB 10/28/2021)  Added error return 'stop_ode' to trap error conditions in 
! derivative subroutine (referred to as f() here, is deriv() in RAYS code)
!
!*****************************************************************************80
!
!! ODE is the user interface to an ordinary differential equation solver.
!
!  Discussion:
!
!    ODE integrates a system of NEQN first order ordinary differential
!    equations of the form:
!      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
!      Y(i) given at T.
!    The subroutine integrates from T to TOUT.  On return, the
!    parameters in the call list are set for continuing the integration.
!    The user has only to define a new value TOUT and call ODE again.
!
!    The differential equations are actually solved by a suite of codes
!    DE, STEP, and INTRP.  ODE allocates virtual storage in the
!    arrays WORK and IWORK and calls DE.  DE is a supervisor which
!    directs the solution.  It calls the routines STEP and INTRP
!    to advance the integration and to interpolate at output points.
!
!    STEP uses a modified divided difference form of the Adams PECE
!    formulas and local extrapolation.  It adjusts the order and step
!    size to control the local error per unit step in a generalized
!    sense.  Normally each call to STEP advances the solution one step
!    in the direction of TOUT.  For reasons of efficiency, DE integrates
!    beyond TOUT internally, though never beyond T+10*(TOUT-T), and
!    calls INTRP to interpolate the solution at TOUT.  An option is
!    provided to stop the integration at TOUT but it should be used
!    only if it is impossible to continue the integration beyond TOUT.
!
!    On the first call to ODE, the user must provide storage in the calling
!    program for the arrays in the call list,
!      Y(NEQN), WORK(100+21*NEQN), IWORK(5),
!    declare F in an external statement, supply the double precision
!      SUBROUTINE F ( T, Y, YP )
!    to evaluate dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
!    and initialize the parameters:
!    * NEQN, the number of equations to be integrated;
!    * Y(1:NEQN), the vector of initial conditions;
!    * T, the starting point of integration;
!    * TOUT, the point at which a solution is desired;
!    * RELERR, ABSERR, the relative and absolute local error tolerances;
!    * IFLAG, an indicator to initialize the code.  Normal input
!      is +1.  The user should set IFLAG = -1 only if it is
!      impossible to continue the integration beyond TOUT.
!    All parameters except F, NEQN and TOUT may be altered by the
!    code on output, and so must be variables in the calling program.
!
!    On normal return from ODE, IFLAG is 2, indicating that T has been
!    set to TOUT, and Y has been set to the approximate solution at TOUT.
!
!    If IFLAG is 3, then the program noticed that RELERR or ABSERR was
!    too small; the output values of these variables are more appropriate,
!    and integration can be resumed by setting IFLAG to 1.
!
!    IFLAG is -2 if the user input IFLAG = -1, and the code was able to
!    reach TOUT exactly.  In that case, the output value of T is TOUT,
!    and the output value of Y is the solution at TOUT, which was computed
!    directly, and not by interpolation.
!
!    Other values of IFLAG generally indicate an error.
!
!    Normally, it is desirable to compute many points along the solution
!    curve.  After the first successful step, more steps may be taken
!    simply by updating the value of TOUT and calling ODE again.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2009
!
!  Author:
!
!    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Lawrence Shampine, Marilyn Gordon,
!    Computer Solution of Ordinary Differential Equations:
!    The Initial Value Problem,
!    Freeman, 1975,
!    ISBN: 0716704617,
!    LC: QA372.S416.
!
!  Parameters:
!
!    Input, external F, the name of a user-supplied routine of the form
!      subroutine f ( t, y, yp )
!      real(KIND=rkind) t
!      real(KIND=rkind) y(neqn)
!      real(KIND=rkind) yp(neqn)
!    which accepts input values T and Y(1:NEQN), evaluates the right hand
!    sides of the ODE, and stores the result in YP(1:NEQN).
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input/output, real(KIND=rkind) Y(NEQN), the current solution.
!
!    Input/output, real(KIND=rkind) T, the current value of the independent
!    variable.
!
!    Input, real(KIND=rkind) TOUT, the desired value of T on output.
!
!    Input, real(KIND=rkind) RELERR, ABSERR, the relative and absolute error
!    tolerances.  At each step, the code requires
!      abs ( local error ) <= abs ( y ) * relerr + abserr
!    for each component of the local error and solution vectors.
!
!    Input/output, integer ( kind = 4 ) IFLAG, indicates the status of 
!    integration.  On input, IFLAG is normally 1 (or -1 in the special case 
!    where TOUT is not to be exceeded.)  On normal output, IFLAG is 2.  Other 
!    output values are:
!    * 3, integration did not reach TOUT because the error tolerances
!      were too small.  But RELERR and ABSERR were increased appropriately
!      for continuing;
!    * 4, integration did not reach TOUT because more than 500 steps were taken;
!    * 5, integration did not reach TOUT because the equations appear to
!      be stiff;
!    * 6, invalid input parameters (fatal error).
!    The value of IFLAG is returned negative when the input value is
!    negative and the integration does not reach TOUT.
!
!    Input/output, real(KIND=rkind) WORK(100+21*NEQN), workspace.
!
!    Input/output, integer ( kind = 4 ) IWORK(5), workspace.
!

  

  use constants_m, only : rkind
  use ode_m, only : ode_stop
    
  implicit none

  integer ( kind = 4 ) neqn

  real(KIND=rkind) abserr
  external f
  
  type(ode_stop)  :: ray_stop

  integer ( kind = 4 ), parameter :: ialpha = 1
  integer ( kind = 4 ), parameter :: ibeta = 13
  integer ( kind = 4 ), parameter :: idelsn = 93
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), parameter :: ig = 62
  integer ( kind = 4 ), parameter :: ih = 89
  integer ( kind = 4 ), parameter :: ihold = 90
  integer ( kind = 4 ) ip
  integer ( kind = 4 ), parameter :: iphase = 75
  integer ( kind = 4 ) iphi
  integer ( kind = 4 ), parameter :: ipsi = 76
  integer ( kind = 4 ), parameter :: isig = 25
  integer ( kind = 4 ), parameter :: istart = 91
  integer ( kind = 4 ), parameter :: itold = 92
  integer ( kind = 4 ), parameter :: iv = 38
  integer ( kind = 4 ), parameter :: iw = 50
  integer ( kind = 4 ) iwt
  integer ( kind = 4 ), parameter :: ix = 88
  integer ( kind = 4 ) iyp
  integer ( kind = 4 ) iypout
  integer ( kind = 4 ), parameter :: iyy = 100
  integer ( kind = 4 ) iwork(5)
  logical nornd
  logical phase1
  real(KIND=rkind) relerr
  logical start
  real(KIND=rkind) t
  real(KIND=rkind) tout
  real(KIND=rkind) work(100+21*neqn)
  real(KIND=rkind) y(neqn)

  iwt = iyy + neqn
  ip = iwt + neqn
  iyp = ip + neqn
  iypout = iyp + neqn
  iphi = iypout + neqn

  if ( abs ( iflag ) /= 1 ) then
    start = ( 0.0e+00 < work(istart) )
    phase1 = ( 0.0e+00 < work(iphase) )
    nornd = ( iwork(2) /= -1 )
  end if

  call de ( f, neqn, y, t, tout, relerr, abserr, iflag, work(iyy), &
    work(iwt), work(ip), work(iyp), work(iypout), work(iphi), &
    work(ialpha), work(ibeta), work(isig), work(iv), work(iw), work(ig), &
    phase1, work(ipsi), work(ix), work(ih), work(ihold), start, &
    work(itold), work(idelsn), iwork(1), nornd, iwork(3), iwork(4), &
    iwork(5), ray_stop)

! Error return (DBB)
  if (ray_stop%stop_ode .eqv. .true.) return
!
  if ( start ) then
    work(istart) = 1.0e+00
  else
    work(istart) = -1.0e+00
  end if

  if ( phase1 ) then
    work(iphase) = 1.0e+00
  else
    work(iphase) = -1.0e+00
  end if

  if ( nornd ) then
    iwork(2) = 1
  else
    iwork(2) = -1
  end if

  return
end

subroutine de ( f, neqn, y, t, tout, relerr, abserr, iflag, yy, wt, p, yp, &
  ypout, phi, alpha, beta, sig, v, w, g, phase1, psi, x, h, hold, start, &
  told, delsgn, ns, nornd, k, kold, isnold, ray_stop)

!*****************************************************************************80
!
!! DE carries out the ODE solution algorithm.
!
!  Discussion:
!
!    ODE merely allocates storage for DE, to relieve the user of the
!    inconvenience of a long call list.  Consequently, DE is used as
!    described in the comments for ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2009
!
!  Author:
!
!    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Lawrence Shampine, Marilyn Gordon,
!    Computer Solution of Ordinary Differential Equations:
!    The Initial Value Problem,
!    Freeman, 1975,
!    ISBN: 0716704617,
!    LC: QA372.S416.
!
!  Parameters:
!
!    Input, external F, the name of a user-supplied routine of the form
!      subroutine f ( t, y, yp )
!      real(KIND=rkind) t
!      real(KIND=rkind) y(neqn)
!      real(KIND=rkind) yp(neqn)
!    which accepts input values T and Y(1:NEQN), evaluates the right hand
!    sides of the ODE, and stores the result in YP(1:NEQN).
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input/output, real(KIND=rkind) Y(NEQN), the current solution.
!
!    Input/output, real(KIND=rkind) T, the current value of the independent
!    variable.
!
!    Input, real(KIND=rkind) TOUT, the desired value of T on output.
!
!    Input, real(KIND=rkind) RELERR, ABSERR, the relative and absolute error
!    tolerances.  At each step, the code requires
!      abs ( local error ) <= abs ( Y ) * RELERR + ABSERR
!    for each component of the local error and solution vectors.
!
!    Input/output, integer ( kind = 4 ) IFLAG, indicates the status of 
!    integration.  On input, IFLAG is normally 1 (or -1 in the special case 
!    where TOUT is not to be exceeded.)  On normal output, IFLAG is 2.  Other 
!    output values are:
!    * 3, integration did not reach TOUT because the error tolerances were
!         too small.
!         But RELERR and ABSERR were increased appropriately for continuing;
!    * 4, integration did not reach TOUT because more than 500 steps were taken;
!    * 5, integration did not reach TOUT because the equations appear to be
!         stiff;
!    * 6, invalid input parameters (fatal error).
!    The value of IFLAG is returned negative when the input value is negative
!    and the integration does not reach TOUT.
!
!    Workspace, real(KIND=rkind) YY(NEQN), used to hold old solution data.
!
!    Input, real(KIND=rkind) WT(NEQN), the error weight vector.
!
!    Workspace, real(KIND=rkind) P(NEQN).
!
!    Workspace, real(KIND=rkind) YP(NEQN), used to hold values of the
!    solution derivative.
!
!    Workspace, real(KIND=rkind) YPOUT(NEQN), used to hold values of the
!    solution derivative.
!
!    Workspace, real(KIND=rkind) PHI(NEQN,16), contains divided difference
!    information about the polynomial interpolant to the solution.
!
!    Workspace, real(KIND=rkind) ALPHA(12), BETA(12), SIG(13).
!
!    Workspace, real(KIND=rkind) V(12), W(12), G(13).
!
!    Input/output, logical PHASE1, indicates whether the program is in the
!    first phase, when it always wants to increase the ODE method order.
!
!    Workspace, real(KIND=rkind) PSI(12), contains information about
!    the polynomial interpolant to the solution.
!
!    Input/output, real(KIND=rkind) X, a "working copy" of T, the current value
!    of the independent variable, which is adjusted as the code attempts
!    to take a step.
!
!    Input/output, real(KIND=rkind) H, the current stepsize.
!
!    Input/output, real(KIND=rkind) HOLD, the last successful stepsize.
!
!    Input/output, logical START, is TRUE on input for the first step.
!    The program initializes data, and sets START to FALSE.
!
!    Input/output, real(KIND=rkind) TOLD, the previous value of T.
!
!    Input/output, real(KIND=rkind) DELSGN, the sign (+1 or -1) of
!    TOUT - T.
!
!    Input/output, integer ( kind = 4 ) NS, the number of steps taken with 
!    stepsize H.
!
!    Input/output, logical NORND, ?
!
!    Input, integer ( kind = 4 ) K, the order of the current ODE method.
!
!    Input, integer ( kind = 4 ) KOLD, the order of the ODE method on the 
!    previous step.
!
!    Input/output, integer ( kind = 4 ) ISNOLD, the previous value of ISN, the 
!    sign of IFLAG.
!
!  Local parameters:
!
!    Local, integer MAXNUM, the maximum number of steps allowed in one
!    call to DE.
!

  use constants_m, only : rkind
  use ode_m, only : ode_stop
  
  implicit none

  integer ( kind = 4 ) neqn
  
  type(ode_stop), intent(out)  :: ray_stop

  real(KIND=rkind) absdel
  real(KIND=rkind) abseps
  real(KIND=rkind) abserr
  real(KIND=rkind) alpha(12)
  real(KIND=rkind) beta(12)
  logical crash
  real(KIND=rkind) del
  real(KIND=rkind) delsgn
  real(KIND=rkind) eps
  external f
  real(KIND=rkind) fouru
  real(KIND=rkind) g(13)
  real(KIND=rkind) h
  real(KIND=rkind) hold
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) isn
  integer ( kind = 4 ) isnold
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kle4
  integer ( kind = 4 ) kold
  integer ( kind = 4 ), parameter :: maxnum = 500
  logical nornd
  integer ( kind = 4 ) nostep
  integer ( kind = 4 ) ns
  real(KIND=rkind) p(neqn)
  real(KIND=rkind) phi(neqn,16)
  logical phase1
  real(KIND=rkind) psi(12)
  real(KIND=rkind) releps
  real(KIND=rkind) relerr
  real(KIND=rkind) sig(13)
  logical start
  logical stiff
  real(KIND=rkind) t
  real(KIND=rkind) tend
  real(KIND=rkind) told
  real(KIND=rkind) tout
  real(KIND=rkind) v(12)
  real(KIND=rkind) w(12)
  real(KIND=rkind) wt(neqn)
  real(KIND=rkind) x
  real(KIND=rkind) y(neqn)
  real(KIND=rkind) yp(neqn)
  real(KIND=rkind) ypout(neqn)
  real(KIND=rkind) yy(neqn)
!
!  Test for improper parameters.
!
  fouru = 4.0e+00 * epsilon ( fouru )

  if ( neqn < 1 ) then
    iflag = 6
    return
  end if

  if ( t == tout ) then
    iflag = 6
  write(*,*) 't == tout'
    return
  end if

  if ( relerr < 0.0e+00 .or. abserr < 0.0e+00 ) then
    iflag = 6
    return
  end if

  eps = max ( relerr, abserr )

  if ( eps <= 0.0e+00 ) then
    iflag = 6
    return
  end if

  if ( iflag == 0 ) then
    iflag = 6
    return
  end if

  isn = sign ( 1, iflag )
  iflag = abs ( iflag )

  if ( iflag /= 1 ) then

    if ( t /= told ) then
      iflag = 6
      return
    end if

    if ( iflag < 2 .or. 5 < iflag ) then
      iflag = 6
      return
    end if
  end if
!
!  On each call set interval of integration and counter for number of
!  steps.  Adjust input error tolerances to define weight vector for
!  subroutine STEP.
!
  del = tout - t
  absdel = abs ( del )

  if ( isn < 0 ) then
    tend = tout
  else
    tend = t + 10.0e+00 * del
  end if

  nostep = 0
  kle4 = 0
  stiff = .false.
  releps = relerr / eps
  abseps = abserr / eps
!
!  On start and restart, also set work variables X and YY(*), store the
!  direction of integration, and initialize the step size.
!
  if ( iflag == 1 .or. isnold < 0 .or. delsgn * del <= 0.0e+00 ) then

    start = .true.
    x = t
    yy(1:neqn) = y(1:neqn)
    delsgn = sign ( real(1.,KIND=rkind), del )
    h = sign ( max ( abs ( tout - x ), fouru * abs ( x ) ), tout - x )

  end if
!
!  If already past the output point, then interpolate and return.
!
  do

    if ( absdel <= abs ( x - t ) ) then
      call intrp ( x, yy, tout, y, ypout, neqn, kold, phi, psi )
      iflag = 2
      t = tout
      told = t
      isnold = isn
      exit
    end if
!
!  If we cannot go past the output point, and we are sufficiently
!  close to it, then extrapolate and return.
!
    if ( isn <= 0 .and. abs ( tout - x ) < fouru * abs ( x ) ) then
      h = tout - x
      call f ( x, yy, yp, ray_stop )
      y(1:neqn) = yy(1:neqn) + h * yp(1:neqn)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      exit
    end if
!
!  Test for too many steps.
!
    if ( maxnum <= nostep ) then
      iflag = isn * 4
      if ( stiff ) then
        iflag = isn * 5
      end if
      y(1:neqn) = yy(1:neqn)
      t = x
      told = t
      isnold = 1
      exit
    end if
!
!  Limit the step size, set the weight vector and take a step.
!
    h = sign ( min ( abs ( h ), abs ( tend - x ) ), h )
    wt(1:neqn) = releps * abs ( yy(1:neqn) ) + abseps

     call step ( x, yy, f, neqn, h, eps, wt, start, &
      hold, k, kold, crash, phi, p, yp, psi, &
      alpha, beta, sig, v, w, g, phase1, ns, nornd, ray_stop)

! Error return (DBB)
    if (ray_stop%stop_ode .eqv. .true.) return
!

!
!  Test for tolerances too small.
!
    if ( crash ) then
      iflag = isn * 3
      relerr = eps * releps
      abserr = eps * abseps
      y(1:neqn) = yy(1:neqn)
      t = x
      told = t
      isnold = 1
      exit
    end if
!
!  Augment the step counter and test for stiffness.
!
    nostep = nostep + 1
    kle4 = kle4 + 1

    if ( 4 < kold ) then
      kle4 = 0
    end if

    if ( 50 <= kle4 ) then
      stiff = .true.
    end if

  end do

  return
end

subroutine step ( x, y, f, neqn, h, eps, wt, start, hold, k, kold, crash, &
  phi, p, yp, psi, alpha, beta, sig, v, w, g, phase1, ns, nornd, ray_stop)

!*****************************************************************************80
!
!! STEP integrates the system of ODE's one step, from X to X+H.
!
!  Discussion:
!
!    This routine integrates a system of first order ordinary differential
!    equations one step, normally from x to x+h, using a modified divided
!    difference form of the Adams PECE formulas.  Local extrapolation is
!    used to improve absolute stability and accuracy.  The code adjusts its
!    order and step size to control the local error per unit step in a
!    generalized sense.  Special devices are included to control roundoff
!    error and to detect when the user is requesting too much accuracy.
!
!    STEP is normally not called directly by the user.  However, it is
!    possible to do so.
!
!    On the first call to STEP, the user must pass in values for:
!    * X, the initial value of the independent variable;
!    * Y, the vector of initial values of dependent variables;
!    * NEQN, the number of equations to be integrated;
!    * H, the nominal step size indicating direction of integration
!      and maximum size of step.  H must be a variable, not a constant;
!    * EPS, the local error tolerance per step.  EPS must be variable;
!    * WT, the vector of non-zero weights for error criterion;
!    * START, set to TRUE.
!
!    STEP requires the L2 norm of the vector with components
!      local error(1:NEQN) / WT(1:NEQN)
!    to be less than EPS for a successful step.  The array WT allows the user
!    to specify an error test appropriate for the problem.  For example,
!    if WT(L):
!    = 1.0, specifies absolute error,
!    = abs(Y(L)), specifies error relative to the most recent value of
!      the L-th component of the solution,
!    = abs(YP(L)), specifies error relative to the most recent value of
!      the L-th component of the derivative,
!    = max (WT(L),abs(Y(L))), specifies error relative to the largest
!      magnitude of L-th component obtained so far,
!    = abs(Y(L))*RELERR/EPS + ABSERR/EPS, specifies a mixed
!      relative-absolute test where EPS = max ( RELERR, ABSERR ).
!
!    On subsequent calls to STEP, the routine is designed so that all
!    information needed to continue the integration, including the next step
!    size H and the next order K, is returned with each step.  With the
!    exception of the step size, the error tolerance, and the weights, none
!    of the parameters should be altered.  The array WT must be updated after
!    each step to maintain relative error tests like those above.
!
!    Normally the integration is continued just beyond the desired endpoint
!    and the solution interpolated there with subroutine INTRP.  If it is
!    impossible to integrate beyond the endpoint, the step size may be
!    reduced to hit the endpoint since the code will not take a step
!    larger than the H input.
!
!    Changing the direction of integration, that is, the sign of H, requires
!    the user to set START = TRUE before calling STEP again.  This is the
!    only situation in which START should be altered.
!
!    A successful step: the subroutine returns after each successful step with
!    START and CRASH both set to FALSE.  X represents the independent variable
!    advanced by one step of length HOLD from its input value; Y has been
!    updated to the solution vector at the new value of X.  All other parameters
!    represent information corresponding to the new X needed to continue the
!    integration.
!
!    Unsuccessful steps: when the error tolerance is too small, the subroutine
!    returns without taking a step and sets CRASH to TRUE. An appropriate step
!    size and error tolerance for continuing are estimated and all other
!    information is restored as upon input before returning.  To continue
!    with the larger tolerance, the user just calls the code again.  A
!    restart is neither required nor desirable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2009
!
!  Author:
!
!    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Lawrence Shampine, Marilyn Gordon,
!    Computer Solution of Ordinary Differential Equations:
!    The Initial Value Problem,
!    Freeman, 1975,
!    ISBN: 0716704617,
!    LC: QA372.S416.
!
!  Parameters:
!
!    Input/output, real(KIND=rkind) X, the value of the independent variable.
!
!    Input/output, real(KIND=rkind) Y(NEQN), the approximate solution at the
!    current value of the independent variable.
!
!    Input, external F, the name of a user-supplied routine of the form
!      subroutine f ( t, y, yp )
!      real(KIND=rkind) t
!      real(KIND=rkind) y(neqn)
!      real(KIND=rkind) yp(neqn)
!    which accepts input values T and Y(1:NEQN), evaluates the right hand
!    sides of the ODE, and stores the result in YP(1:NEQN).
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input/output, real(KIND=rkind) H, the suggested stepsize.
!
!    Input/output, real(KIND=rkind) EPS, the local error tolerance.
!
!    Input, real(KIND=rkind) WT(NEQN), the vector of error weights.
!
!    Input/output, logical START, is set to TRUE before the first step.
!    The program initializes data, and resets START to FALSE.
!
!    Input/output, real(KIND=rkind) HOLD, the step size used on the last
!    successful step.
!
!    Input/output, integer ( kind = 4 ) K, the appropriate order for the 
!    next step.
!
!    Input/output, integer ( kind = 4 ) KOLD, the order used on the last
!    successful step.
!
!    Output, logical CRASH, is set to TRUE if no step can be taken.
!
!    Workspace, real(KIND=rkind) PHI(NEQN,16), contains divided difference
!    information about the polynomial interpolant to the solution.
!
!    Workspace, real(KIND=rkind) P(NEQN).
!
!    Workspace, real(KIND=rkind) YP(NEQN), used to hold values of the
!    solution derivative.
!
!    Workspace, real(KIND=rkind) PSI(12), contains information about
!    the polynomial interpolant to the solution.
!
!    Workspace, real(KIND=rkind) ALPHA(12), BETA(12), SIG(13).
!
!    Workspace, real(KIND=rkind) V(12), W(12), G(13).
!
!    Input/output, logical PHASE1, indicates whether the program is in the
!    first phase, when it always wants to increase the ODE method order.
!
!    Input/output, integer ( kind = 4 ) NS, the number of steps taken with
!    stepsize H.
!
!    Input/output, logical NORND, ?
!

  use constants_m, only : rkind
  use ode_m, only : ode_stop
  
  implicit none

  integer ( kind = 4 ) neqn
  
  type(ode_stop)  :: ray_stop

  real(KIND=rkind) absh
  real(KIND=rkind) alpha(12)
  real(KIND=rkind) beta(12)
  logical crash
  real(KIND=rkind) eps
  real(KIND=rkind) erk
  real(KIND=rkind) erkm1
  real(KIND=rkind) erkm2
  real(KIND=rkind) erkp1
  real(KIND=rkind) err
  external f
  real(KIND=rkind) fouru
  real(KIND=rkind) g(13)
  real(KIND=rkind), dimension ( 13 ) :: gstr = (/ &
    0.50e+00,    0.0833e+00,  0.0417e+00,  0.0264e+00,  0.0188e+00, &
    0.0143e+00,  0.0114e+00,  0.00936e+00, 0.00789e+00, 0.00679e+00, &
    0.00592e+00, 0.00524e+00, 0.00468e+00 /)
  real(KIND=rkind) h
  real(KIND=rkind) hnew
  real(KIND=rkind) hold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) km2
  integer ( kind = 4 ) knew
  integer ( kind = 4 ) kold
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) l
  logical nornd
  integer ( kind = 4 ) ns
  integer ( kind = 4 ) nsp1
  real(KIND=rkind) p(neqn)
  real(KIND=rkind) p5eps
  logical phase1
  real(KIND=rkind) phi(neqn,16)
  real(KIND=rkind) psi(12)
  real(KIND=rkind) r
  real(KIND=rkind) rho
  real(KIND=rkind) round
  real(KIND=rkind) sig(13)
  logical start
  real(KIND=rkind) total
  real(KIND=rkind) tau
  real(KIND=rkind) temp1
  real(KIND=rkind) temp2
  real(KIND=rkind), dimension ( 13 ) :: two = (/ &
       2.0e+00,    4.0e+00,    8.0e+00,  16.0e+00,   32.0e+00, &
      64.0e+00,  128.0e+00,  256.0e+00, 512.0e+00, 1024.0e+00, &
    2048.0e+00, 4096.0e+00, 8192.0e+00/)
  real(KIND=rkind) twou
  real(KIND=rkind) v(12)
  real(KIND=rkind) w(12)
  real(KIND=rkind) wt(neqn)
  real(KIND=rkind) x
  real(KIND=rkind) xold
  real(KIND=rkind) y(neqn)
  real(KIND=rkind) yp(neqn)

  twou = 2.0e+00 * epsilon ( twou )
  fouru = 2.0e+00 * twou
!
!  Check if the step size or error tolerance is too small.  If this is the
!  first step, initialize the PHI array and estimate a starting step size.
!
!  If the step size is too small, determine an acceptable one.
!
  crash = .true.

  if ( abs ( h ) < fouru * abs ( x ) ) then
    h = sign ( fouru * abs ( x ), h )
    return
  end if

  p5eps = 0.5e+00 * eps
!
!  If the error tolerance is too small, increase it to an acceptable value.
!
  round = twou * sqrt ( sum ( ( y(1:neqn) / wt(1:neqn) )**2 ) )

  if ( p5eps < round ) then
    eps = 2.0e+00 * round * ( 1.0e+00 + fouru )
    return
  end if

  crash = .false.
  g(1) = 1.0e+00
  g(2) = 0.5e+00
  sig(1) = 1.0e+00
!
!  Initialize.  Compute an appropriate step size for the first step.
!
  if ( start ) then
  
     call f ( x, y, yp, ray_stop)

! Error return (DBB)
     if (ray_stop%stop_ode .eqv. .true.) return
!  
    phi(1:neqn,1) = yp(1:neqn)
    phi(1:neqn,2) = 0.0e+00
    total = sqrt ( sum ( ( yp(1:neqn) / wt(1:neqn) )**2 ) )
    absh = abs ( h )
    if ( eps < 16.0e+00 * total * h * h ) then
      absh = 0.25D+00 * sqrt ( eps / total )
    end if
    h = sign ( max ( absh, fouru * abs ( x ) ), h )
    hold = 0.0e+00
    k = 1
    kold = 0
    start = .false.
    phase1 = .true.
    nornd = .true.

    if ( p5eps <= 100.0e+00 * round ) then
      nornd = .false.
      phi(1:neqn,15) = 0.0e+00
    end if

  end if

  ifail = 0
!
!  Compute coefficients of formulas for this step.  Avoid computing
!  those quantities not changed when step size is not changed.
!
  do

    kp1 = k + 1
    kp2 = k + 2
    km1 = k - 1
    km2 = k - 2
!
!  NS is the number of steps taken with size H, including the current
!  one.  When K < NS, no coefficients change.
!
    if ( h /= hold ) then
      ns = 0
    end if

    if ( ns <= kold ) then
      ns = ns + 1
    end if

    nsp1 = ns + 1
!
!  Compute those components of ALPHA, BETA, PSI and SIG which change.
!
    if ( ns <= k ) then

      beta(ns) = 1.0e+00
      alpha(ns) = 1.0e+00 / real ( ns, kind = rkind )
      temp1 = h * real ( ns, kind = rkind )
      sig(nsp1) = 1.0e+00

      do i = nsp1, k
        temp2 = psi(i-1)
        psi(i-1) = temp1
        beta(i) = beta(i-1) * psi(i-1) / temp2
        temp1 = temp2 + h
        alpha(i) = h / temp1
        sig(i+1) = real ( i, kind = rkind ) * alpha(i) * sig(i)
      end do

      psi(k) = temp1
!
!  Compute coefficients G.
!
!  Initialize V and set W.
!
      if ( ns <= 1 ) then

        do iq = 1, k
          v(iq) = 1.0e+00 / real ( iq * ( iq + 1 ), kind = rkind )
          w(iq) = v(iq)
        end do
!
!  If order was raised, update the diagonal part of V.
!
      else

        if ( kold < k ) then

          v(k) = 1.0e+00 / real ( k * kp1, kind = rkind )

          do j = 1, ns - 2
            i = k - j
            v(i) = v(i) - alpha(j+1) * v(i+1)
          end do

        end if
!
!  Update V and set W.
!
        do iq = 1, kp1 - ns
          v(iq) = v(iq) - alpha(ns) * v(iq+1)
          w(iq) = v(iq)
        end do
        g(nsp1) = w(1)

      end if
!
!  Compute the G in the work vector W.
!
      do i = ns + 2, kp1
        do iq = 1, kp2 - i
          w(iq) = w(iq) - alpha(i-1) * w(iq+1)
        end do
        g(i) = w(1)
      end do

    end if
!
!  Predict a solution P, evaluate derivatives using predicted
!  solution, estimate local error at order K and errors at orders K,
!  K-1, K-2 as if a constant step size were used.
!
!  Change PHI to PHI star.
!
    do i = nsp1, k
      phi(1:neqn,i) = beta(i) * phi(1:neqn,i)
    end do
!
!  Predict solution and differences.
!
    phi(1:neqn,kp2) = phi(1:neqn,kp1)
    phi(1:neqn,kp1) = 0.0e+00
    p(1:neqn) = 0.0e+00

    do j = 1, k
      i = kp1 - j
      do l = 1, neqn
        p(l) = p(l) + phi(l,i) * g(i)
        phi(l,i) = phi(l,i) + phi(l,i+1)
      end do
    end do

    if ( .not. nornd ) then
      do l = 1, neqn
        tau = h * p(l) - phi(l,15)
        p(l) = y(l) + tau
        phi(l,16) = ( p(l) - y(l) ) - tau
      end do
    else
      p(1:neqn) = y(1:neqn) + h * p(1:neqn)
    end if

    xold = x
    x = x + h
    absh = abs ( h )

    call f ( x, p, yp, ray_stop )

! Error return (DBB)
    if (ray_stop%stop_ode .eqv. .true.) return
!

!
!  Estimate the errors at orders K, K-1 and K-2.
!
    erkm2 = 0.0e+00
    erkm1 = 0.0e+00
    erk = 0.0e+00

    do l = 1, neqn

      if ( 0 < km2 ) then
        erkm2 = erkm2 + ( ( phi(l,km1) + yp(l) - phi(l,1) ) / wt(l) )**2
      end if

      if ( 0 <= km2 ) then
        erkm1 = erkm1 + ( ( phi(l,k) + yp(l) - phi(l,1) ) / wt(l) )**2
      end if

      erk = erk + ( ( yp(l) - phi(l,1) ) / wt(l) )**2

    end do

    if ( 0 < km2 ) then
      erkm2 = absh * sig(km1) * gstr(km2) * sqrt ( erkm2 )
    end if

    if ( 0 <= km2 ) then
      erkm1 = absh * sig(k) * gstr(km1) * sqrt ( erkm1 )
    end if

    err = absh * sqrt ( erk ) * ( g(k) - g(kp1) )
    erk = absh * sqrt ( erk ) * sig(kp1) * gstr(k)
    knew = k
!
!  Test if the order should be lowered.
!
    if ( 0 < km2 ) then

      if ( max ( erkm1, erkm2 ) <= erk ) then
        knew = km1
      end if

    else if ( 0 == km2 ) then

      if ( erkm1 <= 0.5D+00 * erk ) then
        knew = km1
      end if

    end if
!
!  Test if the step was successful.
!
    if ( err <= eps ) then
      exit
    end if
!
!  The step is unsuccessful.  Restore X, PHI and PSI.
!  If third consecutive failure, set order to one.  If the step fails more
!  than three times, consider an optimal step size.  Double the error
!  tolerance and return if the estimated step size is too small for machine
!  precision.
!
!  Restore X, PHI and PSI.
!
    phase1 = .false.
    x = xold
    do i = 1, k
      phi(1:neqn,i) = ( phi(1:neqn,i) - phi(1:neqn,i+1) ) / beta(i)
    end do

    do i = 2, k
      psi(i-1) = psi(i) - h
    end do
!
!  On third failure, set the order to one.  Thereafter, use optimal step size.
!
    ifail = ifail + 1
    temp2 = 0.5D+00

    if ( 3 < ifail ) then

      if ( p5eps < 0.25D+00 * erk ) then
        temp2 = sqrt ( p5eps / erk )
      end if
    end if

    if ( 3 <= ifail ) then
      knew = 1
    end if

    h = temp2 * h
    k = knew

    if ( abs ( h ) < fouru * abs ( x ) ) then
      crash = .true.
      h = sign ( fouru * abs ( x ), h )
      eps = eps + eps
      return
    end if

  end do
!
!  The step is successful.  Correct the predicted solution, evaluate
!  the derivatives using the corrected solution and update the
!  differences.  Determine best order and step size for next step.
!
  kold = k
  hold = h
!
!  Correct and evaluate.
!
  if ( .not. nornd ) then
    do l = 1, neqn
      rho = h * g(kp1) * ( yp(l) - phi(l,1) ) - phi(l,16)
      y(l) = p(l) + rho
      phi(l,15) = ( y(l) - p(l) ) - rho
    end do
  else
    y(1:neqn) = p(1:neqn) + h * g(kp1) * ( yp(1:neqn) - phi(1:neqn,1) )
  end if
  call f ( x, y, yp, ray_stop)
 
! Error return (DBB)
    if (ray_stop%stop_ode .eqv. .true.) return
!
  
!
!  Update differences for the next step.
!
  phi(1:neqn,kp1) = yp(1:neqn) - phi(1:neqn,1)
  phi(1:neqn,kp2) = phi(1:neqn,kp1) - phi(1:neqn,kp2)

  do i = 1, k
    phi(1:neqn,i) = phi(1:neqn,i) + phi(1:neqn,kp1)
  end do
!
!  Estimate error at order K+1 unless:
!  * in first phase when always raise order, or,
!  * already decided to lower order, or,
!  * step size not constant so estimate unreliable.
!
  erkp1 = 0.0e+00

  if ( knew == km1 .or. k == 12 ) then
    phase1 = .false.
  end if

  if ( phase1 ) then

    k = kp1
    erk = erkp1

  else if ( knew == km1 ) then

    k = km1
    erk = erkm1

  else if ( kp1 <= ns ) then

    do l = 1, neqn
      erkp1 = erkp1 + ( phi(l,kp2) / wt(l) )**2
    end do
    erkp1 = absh * gstr(kp1) * sqrt ( erkp1 )
!
!  Using estimated error at order K+1, determine appropriate order
!  for next step.
!
    if ( k == 1 ) then

      if ( erkp1 < 0.5D+00 * erk ) then
        k = kp1
        erk = erkp1
      end if

    else if ( erkm1 <= min ( erk, erkp1 ) ) then

      k = km1
      erk = erkm1

    else if ( erkp1 < erk .and. k < 12 ) then

      k = kp1
      erk = erkp1

    end if

  end if
!
!  With the new order, determine appropriate step size for next step.
!
  hnew = h + h

  if ( .not. phase1 ) then

    if ( p5eps < erk * two(k+1) ) then

      hnew = h

      if ( p5eps < erk ) then
        temp2 = real ( k + 1, kind = rkind )
        r = ( p5eps / erk )**( 1.0e+00 / temp2 )
        hnew = absh * max ( real(0.5, KIND=rkind), min ( real(0.9, KIND=rkind), r ) )
        hnew = sign ( max ( hnew, fouru * abs ( x ) ), h )
      end if

    end if

  end if

  h = hnew

  return
end
subroutine intrp ( x, y, xout, yout, ypout, neqn, kold, phi, psi )

!*****************************************************************************80
!
!! INTRP approximates the solution at XOUT by polynomial interpolation.
!
!  Discussion:
!
!    The methods in STEP approximate the solution near X by a polynomial.
!    This routine approximates the solution at XOUT by evaluating the
!    polynomial there.  Information defining this polynomial is passed
!    from STEP, so INTRP cannot be used alone.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2009
!
!  Author:
!
!    Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Lawrence Shampine, Marilyn Gordon,
!    Computer Solution of Ordinary Differential Equations:
!    The Initial Value Problem,
!    Freeman, 1975,
!    ISBN: 0716704617,
!    LC: QA372.S416.
!
!  Parameters:
!
!    Input, real(KIND=rkind) X, the point where the solution has been computed.
!
!    Input, real(KIND=rkind) Y(NEQN), the computed solution at X.
!
!    Input, real(KIND=rkind) XOUT, the point at which the solution is desired.
!
!    Output, real(KIND=rkind) YOUT(NEQN), the solution at XOUT.
!
!    Output, real(KIND=rkind) YPOUT(NEQN), the derivative of the solution
!    at XOUT.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) KOLD, the order used for the last
!    successful step.
!
!    Input, real(KIND=rkind) PHI(NEQN,16), contains information about the
!    interpolating polynomial.
!
!    Input, real(KIND=rkind) PSI(12), contains information about the
!    interpolating polynomial.
!

  use constants_m, only : rkind

  implicit none

  integer ( kind = 4 ) neqn

  real(KIND=rkind) eta
  real(KIND=rkind) g(13)
  real(KIND=rkind) gamma
  real(KIND=rkind) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) kold
  real(KIND=rkind) phi(neqn,16)
  real(KIND=rkind) psi(12)
  real(KIND=rkind) psijm1
  real(KIND=rkind) rho(13)
  real(KIND=rkind) term
  real(KIND=rkind) w(13)
  real(KIND=rkind) x
  real(KIND=rkind) xout
  real(KIND=rkind) y(neqn)
  real(KIND=rkind) yout(neqn)
  real(KIND=rkind) ypout(neqn)

  hi = xout - x
  ki = kold + 1
!
!  Initialize W for computing G.
!
  do i = 1, ki
    w(i) = 1.0e+00 / real ( i, kind = rkind )
  end do
!
!  Compute G.
!
  g(1) = 1.0e+00
  rho(1) = 1.0e+00
  term = 0.0e+00

  do j = 2, ki
    psijm1 = psi(j-1)
    gamma = ( hi + term ) / psijm1
    eta = hi / psijm1
    do i = 1, ki + 1 - j
      w(i) = gamma * w(i) - eta * w(i+1)
    end do
    g(j) = w(1)
    rho(j) = gamma * rho(j-1)
    term = psijm1
  end do
!
!  Interpolate.
!
  ypout(1:neqn) = 0.0e+00
  yout(1:neqn) = 0.0e+00

  do j = 1, ki
    i = ki + 1 - j
    yout(1:neqn) = yout(1:neqn) + g(i) * phi(1:neqn,i)
    ypout(1:neqn) = ypout(1:neqn) + rho(i) * phi(1:neqn,i)
  end do

  yout(1:neqn) = y(1:neqn) + hi * yout(1:neqn)

  return
end


subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
