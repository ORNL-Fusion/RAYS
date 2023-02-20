 MODULE bin_to_uniform_grid_m

! generic procedure: uniform grid binner to accumulate an extensive quantity, Q, into a
! set of uniformly spaced bins in x.  For example to discretize the power absorbed from
! a wave as it propagates in x, giving the local profile of power deposition.

! Input Q(:): Extensive quantity to be distributed over x bins.
! Input x(:) Strictly increasing parameter vector.  size(x) must equal size(Q))
! Input xmin: Scalar value of left edge of lowest bin
! Input xmax: Scalar value of right edge of highest bin

! Output binned_Q(:): intent(inout).  number of bins taken from size(binned_Q)

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    interface bin_to_uniform_grid 
        module procedure binner_single, binner_real
    end interface
    
    CONTAINS
 !*********************************************************************************

    SUBROUTINE binner_single(Q, x, xmin, xmax, binned_Q, ierr)
    
        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: Q(:), x(:)
        REAL(kind = skind), INTENT(IN) :: xmin, xmax
        REAL(kind = skind), INTENT(INOUT) :: binned_Q(:)
        INTEGER, INTENT(OUT) :: ierr
        
        INTEGER :: nx, n_bins, is
        REAL(kind = skind) :: x_low, x_high, x_range, bin_width
        
        INTEGER :: index_low, index_high, index_temp, i, delta_i
        REAL(kind = skind) :: delta_Q, Q_density, Q_incrL, Q_incrH
        REAL(kind = skind) :: Rindex_low, Rindex_high, Rdelta_index
        REAL(kind = skind) :: fraction_low, fraction_high
        
        ierr = 0
        
! Check that s(:), Q(:), x(:) are all the same size
		if (size(x) .ne. size(Q)) then
			write(*,*) 'uniform_grid_binner: Input arrays sizes inconsistent ', &
			&   '  size(Q) = ', size(Q), '   size(x) = ', size(x)
			ierr = 1
			return
		end if
		
		nx = size(x)
		
! Check that xmin <= all x <= xmax
		x_low = minval(x)
		x_high = maxval(x)
		
		write(*,*) 'x_low, x_high', x_low, x_high
		
		if (x_low < xmin .or. x_high > xmax) then
			write(*,*) 'uniform_grid_binner: x outside range xmin, xmax'
			ierr = 2
			return
		end if
		
		n_bins = size(binned_Q)
        binned_Q = 0.  	! Empty all bins

		x_range = xmax - xmin
		bin_width = x_range/n_bins

		s_loop: do is = 2, nx
			x_low = x(is-1)
			x_high = x(is)
			
			! find bin indices: xmin .le. x < xmin+bin_width -> bin 1
			! xmax-bin_width < x .le. xmax -> bin n_bins
			index_low = floor((x_low - xmin)/bin_width) + 1
			index_high = floor((x_high - xmin)/bin_width) + 1
			if (x_high .ge. xmax) index_high = n_bins ! special case x_high = xmax
			delta_i = index_high - index_low

			if (index_high < index_low) then ! Switch
				index_temp = index_high
				index_high = index_low
				index_low = index_temp
			end if
			
			delta_Q = Q(is) - Q(is-1)

!  Distribute Q
			if (delta_i == 0) then  ! Put all delta_Q in that one bin
				binned_Q(index_low) = binned_Q(index_low) + delta_Q

			else if (delta_i > 0) then  ! Divide between bins
                Rindex_low = (x_low - xmin)/bin_width ! real number index of bottom point
                Rindex_high = (x_high - xmin)/bin_width ! real number index of top point
                Rdelta_index = Rindex_high - Rindex_low ! interval in index space
				
				fraction_low = (real(index_low) - Rindex_low)/Rdelta_index
				Q_incrL = delta_Q*fraction_low ! contribution to bin index_low
				binned_Q(index_low) = binned_Q(index_low) + Q_incrL

				fraction_high = (Rindex_high - real(index_high - 1))/Rdelta_index
				Q_incrH = delta_q*fraction_high ! contribution to bin index_high
				binned_Q(index_high) = binned_Q(index_high) + Q_incrH
			
				if (delta_i > 1) then ! increment middle bins by Q_density
					Q_density = delta_Q/(Rindex_high - Rindex_low)
					do i = index_low + 1, index_high - 1
						binned_Q(i) = binned_Q(i) + Q_density
					end do
				end if
			end if
						
		end do s_loop
    
    RETURN
    END SUBROUTINE binner_single
 !*********************************************************************************
    

    SUBROUTINE binner_real(Q, x, xmin, xmax, binned_Q, ierr)
    
        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: Q(:), x(:)
        REAL(kind = rkind), INTENT(IN) :: xmin, xmax
        REAL(kind = rkind), INTENT(INOUT) :: binned_Q(:)
        INTEGER, INTENT(OUT) :: ierr
        
        INTEGER :: nx, n_bins, is
        REAL(kind = rkind) :: x_low, x_high, x_range, bin_width
        
        INTEGER :: index_low, index_high, index_temp, i, delta_i
        REAL(kind = rkind) :: delta_Q, Q_density, Q_incrL, Q_incrH
        REAL(kind = rkind) :: Rindex_low, Rindex_high, Rdelta_index
        REAL(kind = rkind) :: fraction_low, fraction_high
        
        ierr = 0
        
! Check that s(:), Q(:), x(:) are all the same size
		if (size(x) .ne. size(Q)) then
			write(*,*) 'uniform_grid_binner: Input arrays sizes inconsistent ', &
			&   '  size(Q) = ', size(Q), '   size(x) = ', size(x)
			ierr = 1
			return
		end if
		
		nx = size(x)
		
! Check that xmin <= all x <= xmax
		x_low = minval(x)
		x_high = maxval(x)
		
		write(*,*) 'x_low, x_high', x_low, x_high
		
		if (x_low < xmin .or. x_high > xmax) then
			write(*,*) 'uniform_grid_binner: x outside range xmin, xmax'
			ierr = 2
			return
		end if
		
		n_bins = size(binned_Q)
        binned_Q = 0.  	! Empty all bins

		x_range = xmax - xmin
		bin_width = x_range/n_bins

		s_loop: do is = 2, nx
			x_low = x(is-1)
			x_high = x(is)
			
			! find bin indices: xmin .le. x < xmin+bin_width -> bin 1
			! xmax-bin_width < x .le. xmax -> bin n_bins
			index_low = floor((x_low - xmin)/bin_width) + 1
			index_high = floor((x_high - xmin)/bin_width) + 1
			if (x_high .ge. xmax) index_high = n_bins ! special case x_high = xmax
			delta_i = index_high - index_low

			if (index_high < index_low) then ! Switch
				index_temp = index_high
				index_high = index_low
				index_low = index_temp
			end if
			
			delta_Q = Q(is) - Q(is-1)

!  Distribute Q
			if (delta_i == 0) then  ! Put all delta_Q in that one bin
				binned_Q(index_low) = binned_Q(index_low) + delta_Q

			else if (delta_i > 0) then  ! Divide between bins
                Rindex_low = (x_low - xmin)/bin_width ! real number index of bottom point
                Rindex_high = (x_high - xmin)/bin_width ! real number index of top point
                Rdelta_index = Rindex_high - Rindex_low ! interval in index space
				
				fraction_low = (real(index_low) - Rindex_low)/Rdelta_index
				Q_incrL = delta_Q*fraction_low ! contribution to bin index_low
				binned_Q(index_low) = binned_Q(index_low) + Q_incrL

				fraction_high = (Rindex_high - real(index_high - 1))/Rdelta_index
				Q_incrH = delta_q*fraction_high ! contribution to bin index_high
				binned_Q(index_high) = binned_Q(index_high) + Q_incrH
			
				if (delta_i > 1) then ! increment middle bins by Q_density
					Q_density = delta_Q/bin_width
					do i = index_low + 1, index_high - 1
						binned_Q(i) = binned_Q(i) + Q_density
					end do
				end if
			end if
						
		end do s_loop
    
    RETURN
    END SUBROUTINE binner_real
 !*********************************************************************************
    

        
END MODULE bin_to_uniform_grid_m