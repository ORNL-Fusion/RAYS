MODULE rpoly_m

!*****************************************************************************************
!
!  Finds the zeros of a general real polynomial using the Jenkins & Traub algorithm.
!
!### History
!  * M. A. Jenkins, "[Algorithm 493: Zeros of a Real Polynomial](https://dl.acm.org/doi/10.1145/355637.355643)",
!    ACM Transactions on Mathematical SoftwareVolume 1, Issue 2, June 1975, pp 178-189
!  * code converted using to_f90 by alan miller, 2003-06-02
!  * Jacob Williams, 9/13/2022 : modernized this code
! https://github.com/jacobwilliams/polyroots-fortran/blob/master/src/polyroots_module.F90

!______________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

! DBB 12/5/2025
! Slightly modified (definition of numeric kinds) to match RAYS conventions

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use iso_fortran_env
    use ieee_arithmetic

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    integer, parameter :: wp = rkind  !! a shorter name

    real(wp), parameter :: eps = epsilon(1.0_wp) !! machine epsilon
    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp), parameter :: deg2rad = pi/180.0_wp

contains

subroutine rpoly(op, degree, zeror, zeroi, istat)

    implicit none

    integer, intent(in) :: degree !! degree of polynomial
    real(wp), dimension(degree+1), intent(in) :: op !! vector of coefficients in order of decreasing powers
    real(wp), dimension(degree), intent(out) :: zeror !! output vector of real parts of the zeros
    real(wp), dimension(degree), intent(out) :: zeroi !! output vector of imaginary parts of the zeros
    integer, intent(out) :: istat !! status output:
                                  !!
                                  !! * `0` : success
                                  !! * `-1` : leading coefficient is zero
                                  !! * `-2` : no roots found
                                  !! * `>0` : the number of zeros found

    real(wp), dimension(:), allocatable :: p, qp, k, qk, svk, temp, pt
    real(wp) :: sr, si, u, v, a, b, c, d, a1, a3, &
                a7, e, f, g, h, szr, szi, lzr, lzi, &
                t, aa, bb, cc, factor, mx, mn, xx, yy, &
                xxx, x, sc, bnd, xm, ff, df, dx
    integer :: cnt, nz, i, j, jj, l, nm1, n, nn
    logical :: zerok

    real(wp), parameter :: cosr = cos(94.0_wp*deg2rad)
    real(wp), parameter :: sinr = sin(86.0_wp*deg2rad)
    real(wp), parameter :: base = radix(1.0_wp)
    real(wp), parameter :: eta = eps
    real(wp), parameter :: infin = huge(1.0_wp)
    real(wp), parameter :: smalno = tiny(1.0_wp)
    real(wp), parameter :: sqrthalf = sqrt(0.5_wp)
    real(wp), parameter :: are = eta !! unit error in +
    real(wp), parameter :: mre = eta !! unit error in *
    real(wp), parameter :: lo = smalno/eta

    ! initialization of constants for shift rotation
    xx = sqrthalf
    yy = -xx
    istat = 0
    n = degree
    nn = n + 1

    ! algorithm fails if the leading coefficient is zero.
    if (op(1) == 0.0_wp) then
        istat = -1
        return
    end if

    ! remove the zeros at the origin if any
    do
        if (op(nn) /= 0.0_wp) exit
        j = degree - n + 1
        zeror(j) = 0.0_wp
        zeroi(j) = 0.0_wp
        nn = nn - 1
        n = n - 1
    end do

    ! allocate various arrays
    if (allocated(p)) deallocate (p, qp, k, qk, svk)
    allocate (p(nn), qp(nn), k(nn), qk(nn), svk(nn), temp(nn), pt(nn))

    ! make a copy of the coefficients
    p(1:nn) = op(1:nn)

    main: do

        ! start the algorithm for one zero
        if (n <= 2) then
            if (n < 1) return

            ! calculate the final zero or pair of zeros
            if (n /= 2) then
                zeror(degree) = -p(2)/p(1)
                zeroi(degree) = 0.0_wp
                return
            end if
            call quad(p(1), p(2), p(3), zeror(degree - 1), zeroi(degree - 1), &
                      zeror(degree), zeroi(degree))
            return
        end if

        ! find largest and smallest moduli of coefficients.
        mx = 0.0_wp ! max
        mn = infin  ! min
        do i = 1, nn
            x = abs(real(p(i), wp))
            if (x > mx) mx = x
            if (x /= 0.0_wp .and. x < mn) mn = x
        end do

        ! scale if there are large or very small coefficients computes a scale
        ! factor to multiply the coefficients of the polynomial.
        ! the scaling is done to avoid overflow and to avoid undetected underflow
        ! interfering with the convergence criterion.
        ! the factor is a power of the base
        scale: block
            sc = lo/mn
            if (sc <= 1.0_wp) then
                if (mx < 10.0_wp) exit scale
                if (sc == 0.0_wp) sc = smalno
            else
                if (infin/sc < mx) exit scale
            end if
            l = log(sc)/log(base) + 0.5_wp
            factor = (base*1.0_wp)**l
            if (factor /= 1.0_wp) then
                p(1:nn) = factor*p(1:nn)
            end if
        end block scale

        ! compute lower bound on moduli of zeros.
        pt(1:nn) = abs(p(1:nn))
        pt(nn) = -pt(nn)

        ! compute upper estimate of bound
        x = exp((log(-pt(nn)) - log(pt(1)))/n)
        if (pt(n) /= 0.0_wp) then
            ! if newton step at the origin is better, use it.
            xm = -pt(nn)/pt(n)
            if (xm < x) x = xm
        end if

        ! chop the interval (0,x) until ff <= 0
        do
            xm = x*0.1_wp
            ff = pt(1)
            do i = 2, nn
                ff = ff*xm + pt(i)
            end do
            if (ff > 0.0_wp) then
                x = xm
            else
                exit
            end if
        end do
        dx = x

        ! do newton iteration until x converges to two decimal places
        do
            if (abs(dx/x) <= 0.005_wp) exit
            ff = pt(1)
            df = ff
            do i = 2, n
                ff = ff*x + pt(i)
                df = df*x + ff
            end do
            ff = ff*x + pt(nn)
            dx = ff/df
            x = x - dx
        end do
        bnd = x

        ! compute the derivative as the intial k polynomial
        ! and do 5 steps with no shift
        nm1 = n - 1
        do i = 2, n
            k(i) = (nn - i)*p(i)/n
        end do
        k(1) = p(1)
        aa = p(nn)
        bb = p(n)
        zerok = k(n) == 0.0_wp
        do jj = 1, 5
            cc = k(n)
            if (.not. zerok) then
                ! use scaled form of recurrence if value of k at 0 is nonzero
                t = -aa/cc
                do i = 1, nm1
                    j = nn - i
                    k(j) = t*k(j - 1) + p(j)
                end do
                k(1) = p(1)
                zerok = abs(k(n)) <= abs(bb)*eta*10.0_wp
            else
                ! use unscaled form of recurrence
                do i = 1, nm1
                    j = nn - i
                    k(j) = k(j - 1)
                end do
                k(1) = 0.0_wp
                zerok = k(n) == 0.0_wp
            end if
        end do

        ! save k for restarts with new shifts
        temp(1:n) = k(1:n)

        ! loop to select the quadratic  corresponding to each
        ! new shift
        do cnt = 1, 20
            ! quadratic corresponds to a double shift to a non-real point and its complex
            ! conjugate.  the point has modulus bnd and amplitude rotated by 94 degrees
            ! from the previous shift
            xxx = cosr*xx - sinr*yy
            yy = sinr*xx + cosr*yy
            xx = xxx
            sr = bnd*xx
            si = bnd*yy
            u = -2.0_wp*sr
            v = bnd

            ! second stage calculation, fixed quadratic
            call fxshfr(20*cnt, nz)
            if (nz /= 0) then
                ! the second stage jumps directly to one of the third stage iterations and
                ! returns here if successful.
                ! deflate the polynomial, store the zero or zeros and return to the main
                ! algorithm.
                j = degree - n + 1
                zeror(j) = szr
                zeroi(j) = szi
                nn = nn - nz
                n = nn - 1
                p(1:nn) = qp(1:nn)
                if (nz /= 1) then
                    zeror(j + 1) = lzr
                    zeroi(j + 1) = lzi
                end if
                cycle main
            end if

            ! if the iteration is unsuccessful another quadratic
            ! is chosen after restoring k
            k(1:nn) = temp(1:nn)
        end do

        exit
    end do main

    ! return with failure if no convergence with 20 shifts
    istat = degree - n
    if (istat == 0) istat = -2  ! if not roots found

contains

    subroutine fxshfr(l2, nz)

      !! computes up to  l2  fixed shift k-polynomials, testing for convergence in
      !! the linear or quadratic case.  initiates one of the variable shift
      !! iterations and returns with the number of zeros found.

        integer, intent(in) :: l2 !! limit of fixed shift steps
        integer, intent(out) :: nz !! number of zeros found

        real(wp) :: svu, svv, ui, vi, s, betas, betav, oss, ovv, &
                    ss, vv, ts, tv, ots, otv, tvv, tss
        integer :: type, j, iflag
        logical :: vpass, spass, vtry, stry, skip

        nz = 0
        betav = 0.25_wp
        betas = 0.25_wp
        oss = sr
        ovv = v

        ! evaluate polynomial by synthetic division
        call quadsd(nn, u, v, p, qp, a, b)
        call calcsc(type)
        do j = 1, l2
            ! calculate next k polynomial and estimate v
            call nextk(type)
            call calcsc(type)
            call newest(type, ui, vi)
            vv = vi

            ! estimate s
            ss = 0.0_wp
            if (k(n) /= 0.0_wp) ss = -p(nn)/k(n)
            tv = 1.0_wp
            ts = 1.0_wp
            if (j /= 1 .and. type /= 3) then
                ! compute relative measures of convergence of s and v sequences
                if (vv /= 0.0_wp) tv = abs((vv - ovv)/vv)
                if (ss /= 0.0_wp) ts = abs((ss - oss)/ss)

                ! if decreasing, multiply two most recent convergence measures
                tvv = 1.0_wp
                if (tv < otv) tvv = tv*otv
                tss = 1.0_wp
                if (ts < ots) tss = ts*ots

                ! compare with convergence criteria
                vpass = tvv < betav
                spass = tss < betas
                if (spass .or. vpass) then

                    ! at least one sequence has passed the convergence test.
                    ! store variables before iterating
                    svu = u
                    svv = v
                    svk(1:n) = k(1:n)
                    s = ss

                    ! choose iteration according to the fastest converging sequence
                    vtry = .false.
                    stry = .false.
                    skip = (spass .and. ((.not. vpass) .or. tss < tvv))

                    do

                        try: block

                            if (.not. skip) then
                                call quadit(ui, vi, nz)
                                if (nz > 0) return

                                ! quadratic iteration has failed. flag that it has
                                ! been tried and decrease the convergence criterion.
                                vtry = .true.
                                betav = betav*0.25_wp

                                ! try linear iteration if it has not been tried and
                                ! the s sequence is converging
                                if (stry .or. (.not. spass)) exit try
                                k(1:n) = svk(1:n)
                            end if
                            skip = .false.

                            call realit(s, nz, iflag)
                            if (nz > 0) return

                            ! linear iteration has failed.  flag that it has been
                            ! tried and decrease the convergence criterion
                            stry = .true.
                            betas = betas*0.25_wp
                            if (iflag /= 0) then
                                ! if linear iteration signals an almost double real
                                ! zero attempt quadratic interation
                                ui = -(s + s)
                                vi = s*s
                                cycle
                            end if

                        end block try

                        ! restore variables
                        u = svu
                        v = svv
                        k(1:n) = svk(1:n)

                        ! try quadratic iteration if it has not been tried
                        ! and the v sequence is converging
                        if (.not. (vpass .and. (.not. vtry))) exit

                    end do

                    ! recompute qp and scalar values to continue the second stage
                    call quadsd(nn, u, v, p, qp, a, b)
                    call calcsc(type)
                end if
            end if
            ovv = vv
            oss = ss
            otv = tv
            ots = ts
        end do

    end subroutine fxshfr

    subroutine quadit(uu, vv, nz)

         !! variable-shift k-polynomial iteration for a quadratic factor, converges
         !! only if the zeros are equimodular or nearly so.

        real(wp), intent(in) :: uu !! coefficients of starting quadratic
        real(wp), intent(in) :: vv !! coefficients of starting quadratic
        integer, intent(out) :: nz !! number of zero found

        real(wp) :: ui, vi, mp, omp, ee, relstp, t, zm
        integer :: type, i, j
        logical :: tried

        nz = 0
        tried = .false.
        u = uu
        v = vv
        j = 0

        ! main loop
        main: do
            call quad(1.0_wp, u, v, szr, szi, lzr, lzi)

            ! return if roots of the quadratic are real and not
            ! close to multiple or nearly equal and  of opposite sign.
            if (abs(abs(szr) - abs(lzr)) > 0.01_wp*abs(lzr)) return

            ! evaluate polynomial by quadratic synthetic division
            call quadsd(nn, u, v, p, qp, a, b)
            mp = abs(a - szr*b) + abs(szi*b)

            ! compute a rigorous  bound on the rounding error in evaluting p
            zm = sqrt(abs(v))
            ee = 2.0_wp*abs(qp(1))
            t = -szr*b
            do i = 2, n
                ee = ee*zm + abs(qp(i))
            end do
            ee = ee*zm + abs(a + t)
            ee = (5.0_wp*mre + 4.0_wp*are)*ee - &
                 (5.0_wp*mre + 2.0_wp*are)*(abs(a + t) + &
                                            abs(b)*zm) + 2.0_wp*are*abs(t)

            ! iteration has converged sufficiently if the
            ! polynomial value is less than 20 times this bound
            if (mp <= 20.0_wp*ee) then
                nz = 2
                return
            end if
            j = j + 1

            ! stop iteration after 20 steps
            if (j > 20) return
            if (j >= 2) then
                if (.not. (relstp > 0.01_wp .or. mp < omp .or. tried)) then

                    ! a cluster appears to be stalling the convergence.
                    ! five fixed shift steps are taken with a u,v close to the cluster
                    if (relstp < eta) relstp = eta
                    relstp = sqrt(relstp)
                    u = u - u*relstp
                    v = v + v*relstp
                    call quadsd(nn, u, v, p, qp, a, b)
                    do i = 1, 5
                        call calcsc(type)
                        call nextk(type)
                    end do
                    tried = .true.
                    j = 0
                end if
            end if
            omp = mp

            ! calculate next k polynomial and new u and v
            call calcsc(type)
            call nextk(type)
            call calcsc(type)
            call newest(type, ui, vi)

            ! if vi is zero the iteration is not converging
            if (vi == 0.0_wp) exit
            relstp = abs((vi - v)/vi)
            u = ui
            v = vi

        end do main

    end subroutine quadit

    subroutine realit(sss, nz, iflag)

         !! variable-shift h polynomial iteration for a real zero.

        real(wp), intent(inout) :: sss !! starting iterate
        integer, intent(out) :: nz !! number of zero found
        integer, intent(out) :: iflag !! flag to indicate a pair of zeros near real axis.

        real(wp) :: pv, kv, t, s, ms, mp, omp, ee
        integer :: i, j

        nz = 0
        s = sss
        iflag = 0
        j = 0

        ! main loop
        main: do
            pv = p(1)

            ! evaluate p at s
            qp(1) = pv
            do i = 2, nn
                pv = pv*s + p(i)
                qp(i) = pv
            end do
            mp = abs(pv)

            ! compute a rigorous bound on the error in evaluating p
            ms = abs(s)
            ee = (mre/(are + mre))*abs(qp(1))
            do i = 2, nn
                ee = ee*ms + abs(qp(i))
            end do

            ! iteration has converged sufficiently if the
            ! polynomial value is less than 20 times this bound
            if (mp <= 20.0_wp*((are + mre)*ee - mre*mp)) then
                nz = 1
                szr = s
                szi = 0.0_wp
                return
            end if
            j = j + 1

            ! stop iteration after 10 steps
            if (j > 10) return
            if (j >= 2) then
                if (abs(t) <= 0.001_wp*abs(s - t) .and. mp > omp) then
                    ! a cluster of zeros near the real axis has been encountered,
                    ! return with iflag set to initiate a quadratic iteration
                    iflag = 1
                    sss = s
                    return
                end if
            end if

            ! return if the polynomial value has increased significantly
            omp = mp

            ! compute t, the next polynomial, and the new iterate
            kv = k(1)
            qk(1) = kv
            do i = 2, n
                kv = kv*s + k(i)
                qk(i) = kv
            end do
            if (abs(kv) > abs(k(n))*10.0_wp*eta) then
                ! use the scaled form of the recurrence if the value of k at s is nonzero
                t = -pv/kv
                k(1) = qp(1)
                do i = 2, n
                    k(i) = t*qk(i - 1) + qp(i)
                end do
            else
                ! use unscaled form
                k(1) = 0.0_wp
                do i = 2, n
                    k(i) = qk(i - 1)
                end do
            end if
            kv = k(1)
            do i = 2, n
                kv = kv*s + k(i)
            end do
            t = 0.0_wp
            if (abs(kv) > abs(k(n))*10.*eta) t = -pv/kv
            s = s + t

        end do main

    end subroutine realit

    subroutine calcsc(type)

         !! this routine calculates scalar quantities used to
         !! compute the next k polynomial and new estimates of
         !! the quadratic coefficients.

        integer, intent(out) :: type !! integer variable set here indicating how the
                                     !! calculations are normalized to avoid overflow

        ! synthetic division of k by the quadratic 1,u,v
        call quadsd(n, u, v, k, qk, c, d)
        if (abs(c) <= abs(k(n))*100.0_wp*eta) then
            if (abs(d) <= abs(k(n - 1))*100.0_wp*eta) then
                type = 3
                ! type=3 indicates the quadratic is almost a factor of k
                return
            end if
        end if

        if (abs(d) >= abs(c)) then
            type = 2
            ! type=2 indicates that all formulas are divided by d
            e = a/d
            f = c/d
            g = u*b
            h = v*b
            a3 = (a + g)*e + h*(b/d)
            a1 = b*f - a
            a7 = (f + u)*a + h
        else
            type = 1
            ! type=1 indicates that all formulas are divided by c
            e = a/c
            f = d/c
            g = u*e
            h = v*b
            a3 = a*e + (h/c + g)*b
            a1 = b - a*(d/c)
            a7 = a + g*d + h*f
        end if

    end subroutine calcsc

    subroutine nextk(type)

         !! computes the next k polynomials using scalars computed in calcsc.

        integer, intent(in) :: type

        real(wp) :: temp
        integer :: i

        if (type /= 3) then
            temp = a
            if (type == 1) temp = b
            if (abs(a1) <= abs(temp)*eta*10.0_wp) then
                ! if a1 is nearly zero then use a special form of the recurrence
                k(1) = 0.0_wp
                k(2) = -a7*qp(1)
                do i = 3, n
                    k(i) = a3*qk(i - 2) - a7*qp(i - 1)
                end do
                return
            end if

            ! use scaled form of the recurrence
            a7 = a7/a1
            a3 = a3/a1
            k(1) = qp(1)
            k(2) = qp(2) - a7*qp(1)
            do i = 3, n
                k(i) = a3*qk(i - 2) - a7*qp(i - 1) + qp(i)
            end do

        else
            ! use unscaled form of the recurrence if type is 3
            k(1) = 0.0_wp
            k(2) = 0.0_wp
            do i = 3, n
                k(i) = qk(i - 2)
            end do
        end if

    end subroutine nextk

    subroutine newest(type, uu, vv)

        ! compute new estimates of the quadratic coefficients
        ! using the scalars computed in calcsc.

        integer, intent(in) :: type
        real(wp), intent(out) :: uu
        real(wp), intent(out) :: vv

        real(wp) :: a4, a5, b1, b2, c1, c2, c3, c4, temp

        ! use formulas appropriate to setting of type.
        if (type /= 3) then
            if (type /= 2) then
                a4 = a + u*b + h*f
                a5 = c + (u + v*f)*d
            else
                a4 = (a + g)*f + h
                a5 = (f + u)*c + v*d
            end if

            ! evaluate new quadratic coefficients.
            b1 = -k(n)/p(nn)
            b2 = -(k(n - 1) + b1*p(n))/p(nn)
            c1 = v*b2*a1
            c2 = b1*a7
            c3 = b1*b1*a3
            c4 = c1 - c2 - c3
            temp = a5 + b1*a4 - c4
            if (temp /= 0.0_wp) then
                uu = u - (u*(c3 + c2) + v*(b1*a1 + b2*a7))/temp
                vv = v*(1.0_wp + c4/temp)
                return
            end if
        end if

        ! if type=3 the quadratic is zeroed
        uu = 0.0_wp
        vv = 0.0_wp

    end subroutine newest

    subroutine quadsd(nn, u, v, p, q, a, b)

        ! divides `p` by the quadratic `1,u,v` placing the
        ! quotient in `q` and the remainder in `a,b`.

        integer, intent(in) :: nn
        real(wp), intent(in) :: u, v, p(nn)
        real(wp), intent(out) :: q(nn), a, b

        real(wp) :: c
        integer :: i

        b = p(1)
        q(1) = b
        a = p(2) - u*b
        q(2) = a
        do i = 3, nn
            c = p(i) - u*a - v*b
            q(i) = c
            b = a
            a = c
        end do

    end subroutine quadsd

    subroutine quad(a, b1, c, sr, si, lr, li)

         !! calculate the zeros of the quadratic a*z**2+b1*z+c.
         !! the quadratic formula, modified to avoid overflow, is used to find the
         !! larger zero if the zeros are real and both zeros are complex.
         !! the smaller real zero is found directly from the product of the zeros c/a.

        real(wp), intent(in) :: a, b1, c
        real(wp), intent(out) :: sr, si, lr, li

        real(wp) :: b, d, e

        if (a == 0.0_wp) then
            sr = 0.0_wp
            if (b1 /= 0.0_wp) sr = -c/b1
            lr = 0.0_wp
            si = 0.0_wp
            li = 0.0_wp
            return
        end if

        if (c == 0.0_wp) then
            sr = 0.0_wp
            lr = -b1/a
            si = 0.0_wp
            li = 0.0_wp
            return
        end if

        ! compute discriminant avoiding overflow
        b = b1/2.0_wp
        if (abs(b) >= abs(c)) then
            e = 1.0_wp - (a/b)*(c/b)
            d = sqrt(abs(e))*abs(b)
        else
            e = a
            if (c < 0.0_wp) e = -a
            e = b*(b/abs(c)) - e
            d = sqrt(abs(e))*sqrt(abs(c))
        end if

        if (e >= 0.0_wp) then
            ! real zeros
            if (b >= 0.0_wp) d = -d
            lr = (-b + d)/a
            sr = 0.0_wp
            if (lr /= 0.0_wp) sr = (c/lr)/a
            si = 0.0_wp
            li = 0.0_wp
            return
        end if

        ! complex conjugate zeros
        sr = -b/a
        lr = sr
        si = abs(d/a)
        li = -si

    end subroutine quad

end subroutine rpoly
!*****************************************************************************************

END MODULE rpoly_m
