 subroutine deallocate_array
!   allocates arrays which are never deallocated in the main loop.

    use equilibrium_m, only : ns, gradns, ts, gradts, omgc, omgp2, alpha, gamma
    use suscep_m, only : chis_cold, chis

    deallocate( ns, gradns )
    deallocate( ts, gradts )
    deallocate( omgc, omgp2 )
    deallocate( alpha, gamma )

    deallocate( chis_cold, chis )

    return
 end subroutine deallocate_array
