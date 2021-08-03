subroutine fr_speel_z(natm,invlx,invly,invlz, alpha,r4pie,chge,rz, kmax, &
                   &  fz_io,eng_io,enc,ens)
  implicit none
  integer:: natm, kmax
  real(8):: invlx,invly,invlz,alpha,r4pie
  real(8),dimension(1:natm):: chge,rz

  real(8),dimension(1:natm):: fz_io
  real(8):: eng_io

  real(8),dimension(1:natm,1:kmax):: enc,ens
  !------------------------------------------------------------------------------
  !
  ! fr_speel_z():
  !
  ! subroutine to compute forces for symmetry-preserving  erf-electrostatics
  !                       f r_       s        p           e   el_
  !
  !                       in the z direction
  !                              z
  !
  ! fr_speel_z() take enc, ens as dummy input for use of -Ofast option
  !
  ! fr_speel_z() is called after the real space erfc-electrostatics and lj
  !
  ! there is no initialization of force in fr_speel_z()
  !
  ! j runs over 1 to natm (number of charges)
  !
  ! k runs over 1 to kmax (number of kz vectors)
  !
  !
  ! this is the SPMF code, which is much better than LMF/LRT or LMF/EXP
  !                                 fast, straightfoward, and no iteration
  !
  ! Ref: Chem. Commun. (2014), 50, 14397
  !
  !------------------------------------------------------------------------------
  !
  ! derived for inverse length
  real(8),parameter:: tpi=6.283185307179586d0
  real(8),dimension(1:natm):: ckc,cks
  real(8):: r4alpsq, ssz,rkz,rkzsq,ckcs,ckss,akk
  real(8):: engcpe,engcoeffi, engc2
  integer:: j,k
  !
  ! two large and too small alpha are not acceptable
  if(alpha.lt.1.d-6 .or. alpha.gt.100.0d0) stop 'too large or too small alpha'
  !
  ! set up 1/(4 \alpha^2) and engcoefficient
  r4alpsq=-0.25d0/alpha**2; engcoeffi = 2.0d0 * tpi*(invlx*invly*invlz)*r4pie
  engc2 = 2.0d0 * engcoeffi
  !
  ! calculate and store sin(kz) and cos(kz) for each atom, each k
  do j=1, natm
    ssz = rz(j)*invlz; enc(j,1)=cos(tpi*ssz); ens(j,1)=sin(tpi*ssz)
  end do
  do k = 2, kmax
    do j = 1, natm
      enc(j,k)=enc(j,k-1)*enc(j,1)-ens(j,k-1)*ens(j,1)
      ens(j,k)=ens(j,k-1)*enc(j,1)+enc(j,k-1)*ens(j,1)
    end do
  end do
  !
  ! sum over k vectors for energy and force
  engcpe = 0.0d0
  do k=1,kmax
    rkz=tpi*dble(k)*invlz; rkzsq = rkz ** 2
    ! calculate exp(ikr) terms and product with charges
    do j = 1, natm
      ckc(j)=chge(j) * enc(j,k)
      cks(j)=chge(j) * ens(j,k)
    end do
    ckcs = sum(ckc(1:natm)); ckss = sum(cks(1:natm))
    akk=exp(r4alpsq*rkzsq)/rkzsq
    !
    ! accumulate potential energy, force and virial terms (not now)
    engcpe=engcpe+akk*(ckcs ** 2 +ckss ** 2)
    do j = 1, natm
      fz_io(j) = fz_io(j) + engc2 * akk*(cks(j)*ckcs-ckc(j)*ckss)*rkz
    end do
  end do
  !
  ! calculate final energy
  eng_io=eng_io + engcoeffi*engcpe
  return
end subroutine fr_speel_z
