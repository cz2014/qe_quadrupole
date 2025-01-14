!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine add_zstar_ue (imode0, npe)
  !-----------------------------------------------------------------------
  ! add the contribution of the modes imode0+1 -> imode+npe
  ! to the effective charges Z(Us,E) (Us=scf,E=bare)
  !
  ! trans =.true. is needed for this calculation to be meaningful
  !
  USE kinds, only : DP
  USE klist, ONLY : xk, wk, ngk, igk_k
  USE uspp,  ONLY : vkb
  USE wvfct, ONLY : npwx
  USE wavefunctions,  ONLY: evc
  USE noncollin_module,      ONLY: noncolin
  USE buffers,  ONLY : get_buffer
  USE qpoint,   ONLY: nksq, ikks, ikqs
  USE eqv,      ONLY: dpsi, dvpsi
  USE efield_mod, ONLY: zstarue0_rec
  USE units_ph,   ONLY : iudwf, lrdwf
  USE units_lr,   ONLY : iuwfc, lrwfc
  USE control_lr, ONLY : nbnd_occ

  implicit none

  integer, intent(in) :: imode0, npe

  integer :: ibnd, jpol, ipert, nrec, mode, ik, ikk, ikq ! -cz ikq
  ! counter on bands
  ! counter on polarization
  ! counter on pertubations
  ! counter on records
  ! counter on modes
  ! counter on k points
  INTEGER :: npw, npwq

  real(DP) :: weight

  call start_clock('add_zstar_ue')
  zstarue0_rec=(0.0_DP,0.0_DP)
  do ik = 1, nksq
     ikk=ikks(ik)
     ikq=ikqs(ik)
     npw = ngk(ikk)
     npwq = ngk(ikq)
     weight = wk (ikk)
     if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ikk)
     call init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)
     do jpol = 1, 3
        !
        ! read/compute DeltaV*psi(bare) for electric field
        !
        call dvpsi_e (ik, jpol)
        !
        do ipert = 1, npe
           mode = imode0 + ipert
           nrec = (ipert - 1) * nksq + ik
           !
           ! read dpsi(scf)/du for phonon mode # mode
           !
           call get_buffer (dpsi, lrdwf, iudwf, nrec)
           do ibnd = 1, nbnd_occ(ik)
              zstarue0_rec (mode, jpol) = zstarue0_rec (mode, jpol) - 2.d0 * weight * &
                   dot_product (dpsi(:,ibnd), dvpsi(:,ibnd))
              IF (noncolin) &
                 zstarue0_rec(mode,jpol)=zstarue0_rec (mode, jpol) - 2.d0 * weight * &
                   dot_product (dpsi(1+npwx:npw+npwx,ibnd), dvpsi(1+npwx:npw+npwx,ibnd))

           enddo
        enddo
     enddo

  enddo
  call stop_clock('add_zstar_ue')
  return
end subroutine add_zstar_ue

!-----------------------------------------------------------------------
subroutine add_zstar_ue2 (imode0, npe)
  !-----------------------------------------------------------------------
  ! add the contribution of the modes imode0+1 -> imode+npe
  ! to the charge repsonce Q(q, \kappa \beta)
  !
  ! modified from add_zstar_ue -cz
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi
  USE cell_base, ONLY : at, omega
  USE fft_base, ONLY : dfftp
  USE fft_types, ONLY : fft_index_to_3d
  USE lsda_mod, ONLY : nspin
  USE noncollin_module, ONLY : nspin_mag
  USE efield_mod, ONLY : zstarue0_rec
  USE qpoint, ONLY : xq
  USE units_ph, ONLY : lrdrho, iudrho
  USE output, ONLY : fildrho

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: imode0, npe

  LOGICAL :: offrange
  INTEGER :: ir_end, ir, i, j, k, ipert, mode
  REAL(DP) :: arg, xq_cry(3)
  COMPLEX(DP) :: dres
  COMPLEX(DP), ALLOCATABLE :: aux(:,:,:), eiqr(:)

  IF (fildrho == ' ') CALL errore('FIXME', 'fildrho should be set up', 1)
  IF (nspin_mag /= 1) CALL errore('FIXME', 'magnetism not implemented', 1)
  IF (nspin == 2) CALL errore('FIXME', 'LSDA magnetism not implemented', 1)

  ALLOCATE(aux(dfftp%nnr,nspin_mag,npe))
  ! drho for each mode
  ALLOCATE(eiqr(dfftp%nnr))
  ! phase factor

  call start_clock('add_zstar_ue')
! xq_cry = xq
! CALL cryst_to_cart (1, xq_cry, at, -1)
  ! obtain xq in crystal basis


  !! here we read drho
  DO ipert = 1, npe
      mode = imode0 + ipert
      CALL davcio_drho (aux (1,1,ipert), lrdrho, iudrho, mode, -1)
  END DO

  !! here we try to set up eiqr
  eiqr = (1.d0, 0.d0)
  !CALL multiply_iqr(dfftp, -xq, eiqr)
#if defined (__MPI)
  ir_end = MIN(dfftp%nnr, dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)
#else
  ir_end = dfftp%nnr
#endif
  !
! DO ir = 1, ir_end
!   !
!   CALL fft_index_to_3d(ir, dfftp, i, j, k, offrange)
!   IF ( offrange ) CYCLE
!   !
!   ! (i,j,k) is the zero-based coordinate of the real-space grid
!   arg = -tpi * (  xq_cry(1) * REAL(i, DP) / REAL(dfftp%nr1, DP) &
!                + xq_cry(2) * REAL(j, DP) / REAL(dfftp%nr2, DP) &
!                + xq_cry(3) * REAL(k, DP) / REAL(dfftp%nr3, DP)  )
!  !eiqr(ir) = CMPLX( COS(arg), SIN(arg), kind=DP )
!   eiqr(ir) = (1.d0, 0.d0)
!   !
! END DO

  !! here we multiply phase factor and integrate in unit cell
  zstarue0_rec=(0.0_DP,0.0_DP)
  DO ipert = 1, npe
      mode = imode0 + ipert
     !CALL davcio_drho (aux (1,1), lrdrho, iudrho, mode, -1)
      ! read # mode from drho files
      dres = dot_product (eiqr(1:dfftp%nnr), aux(1:dfftp%nnr,1,ipert))
      dres = dres*omega / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
      zstarue0_rec(mode,1) = -1.d0*dres
      ! here we only multiply phase factor for ispin = 1
  END DO
  !
  DEALLOCATE(aux)
  DEALLOCATE(eiqr)
  !
  call stop_clock('add_zstar_ue')
  !
  RETURN
  !
END SUBROUTINE add_zstar_ue2

