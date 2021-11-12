!----------------------------------------------------------------------------------------
!     READS the transition dipole matrix elements and atomic projections.
!     calculates the Dielectric tensor
!      
!     
!     Drude-Lorentz part of the dielectric tensor has not been added
!     This part is included in QE and is accessible via the flag
!     "intrasmear"
!
!     TDMatrix.dat : Transition Dipole Matrix elements written
!                    by epsilon.x
!     nk  ibnd jbnd  Dx_real  Dx_im  Dy_real  Dy_im  Dz_real  Dz_im   
! 
!-----------------------------------------------------------------------------------------

      program main
      implicit none
      INTEGER :: iL, imu, ik, iq, ipol, jpol, ic, iv, ibnd, jbnd, ir
      INTEGER :: i_rank, i_fi, n_rank, n_fi, ikf, ikr, nkf, nkr
      INTEGER :: nat, nmodes, n_L, nk, nq, nbnd, nl, nr
      INTEGER :: nbnd_min, nbnd_max
      REAL*8 :: e_L, e_L_min, de_L 
      REAL*8 :: vb_min,vb_max,cb_min,cb_max
      REAL*8 :: et, delta, efermi
      REAL*8 :: factor, const, vol
      REAL*8, allocatable :: eigval(:,:), focc (:,:)
      INTEGER, allocatable :: nvb(:),vb_i(:),vb_a(:),cb_i(:),cb_a(:)


      REAL*8 :: full_occ
  
      COMPLEX*8 :: fac
      COMPLEX*8, allocatable :: TDMat (:,:,:,:)
      COMPLEX*8, allocatable :: eps(:,:,:) 

      COMPLEX*8, allocatable :: TDMat_ib (:,:,:)
      COMPLEX*8, allocatable :: eps_ib(:,:,:)
      REAL*8 :: sigma, eta, xarg
      REAL*8, allocatable :: dfocc (:,:)


      REAL*8 :: i, j, k, num1, num2, num3, num4, num5, num6, num7, num8
      CHARACTER (len=200) :: line, filout, dir, code_flag
      CHARACTER (len=200) :: filout_xx, filout_yy, filout_zz

      COMPLEX*8 :: G1_1, G2_1, G3_1, G1_2, G2_2, G3_2
      INTEGER :: ibq, jbq, ikd, ibr, nbq, nbr

      INTEGER :: l_a, l_b
!projection
      INTEGER :: natmwfc, natmwfc_v,natmwfc_c,iw,ib
      INTEGER, allocatable :: atmwfc_v(:),atmwfc_c(:)
      REAL*8, allocatable :: atmwfc(:,:,:),  atmwfc_vc(:,:,:)

      
      REAL :: start, finish
      REAL*8, PARAMETER :: Pi = 3.1415927D0
      REAL*8, PARAMETER :: bohr_to_ang = 0.529177D0
      REAL*8, PARAMETER :: thz_to_eV = 0.00413567D0
      REAL*8, PARAMETER :: thz_to_cm = 33.35643D0 
      REAL*8, PARAMETER :: meV_to_eV = 0.001D0
      REAL*8, PARAMETER :: ry_to_eV = 13.605698D0
  
!---------------------------------------------!
      call cpu_time(start)
      open (21, file='input.dat', status='unknown')    
       read (21,*) dir
       read (21,*) nat
       read (21,*) code_flag
       if ( code_flag .eq. 'QE' ) then
        read (21,*) nk
       else if ( code_flag .eq. 'EPW' ) then
        read (21,*) nk, n_rank, n_fi
       else
        write (6,*) "ONLY QE AND EPW CODES SUPPORTED"
        CALL EXIT
       endif
       read (21,*) nbnd
       read (21,*) nbnd_min, nbnd_max
       read (21,*) nq
       read (21,*) efermi
       read (21,*) vol
       read (21,*) delta, eta, sigma 
       read (21,*) e_L_min, de_L, n_L
       read (21,*) natmwfc
       read (21,*) natmwfc_v
      allocate (atmwfc_v(natmwfc_v))
       read (21,*) (atmwfc_v(iv), iv =1,natmwfc_v)
       read (21,*) natmwfc_c
      allocate (atmwfc_c(natmwfc_c))
       read (21,*) (atmwfc_c(ic), ic =1,natmwfc_c)

      close (21)
      write (6,*) "INPUT.DAT READ"

!-----------------------------------------------
!     COARSE GRID
!-----------------------------------------------
      full_occ = 2.0D0
!     if (nspin .eq. 2) then 
!       full_occ = 1.0D0
!     else 
!       full_occ = 2.0D0
!     endif
!--------------------------
!--------------------------
!       Read projwfc_up.dat 
!--------------------------
      open (22, file=TRIM(ADJUSTL(dir))//"projwfc_up.dat", &
            status='unknown')    ! projection file
      allocate (atmwfc(natmwfc, nk,nbnd))

      do iw = 1, natmwfc
        read (22,*)
        do ik = 1,nk
          do ib = 1,nbnd
            read (22,*) i,j, atmwfc(iw,ik,ib)
          enddo
        enddo
      enddo
      write(6,*) "PROJWFC_UP.DAT READ"

      allocate (atmwfc_vc(2, nk,nbnd))
      atmwfc_vc(:, :,:) = 0.0d0
      do iv = 1,natmwfc_v
        atmwfc_vc(1, :,:) = atmwfc_vc(1, :,:) + atmwfc(atmwfc_v(iv),:,:)
      enddo
      do ic = 1,natmwfc_c
        atmwfc_vc(2, :,:) = atmwfc_vc(2, :,:) + atmwfc(atmwfc_c(ic),:,:)
      enddo

!-----------------------------------------------
!     TDMAT (COARSE GRID)
!-----------------------------------------------
      !
      allocate (TDMat(nk,nbnd,nbnd,3))
      allocate (TDMat_ib(nk,nbnd,3))
      TDMat(:,:,:,:) = CMPLX(0.0D0, 0.0D0)
      !
      open (21, file=TRIM(ADJUSTL(dir))//"TDMatrix.dat", &
            status='unknown')    ! TransitionDipole Matrix elements
      !
      !
      if ( code_flag .eq. 'QE' ) then ! for QE
       !
!--------------------------------------------------------------
!     adapted QE format
!--------------------------------------------------------------

      allocate (eigval(nk,nbnd))
      allocate (focc(nk,nbnd))

!--------------------------------------------------------------
        open (22, file=TRIM(ADJUSTL(dir))//"Eig_Occ.dat", &
               status='old')
        do ik = 1, nk
          do ibnd = 1, nbnd
           read (22,*) num1,num2, eigval(ik,ibnd), focc(ik,ibnd)
          enddo
        enddo
        close (22)
        eigval(:,:) = eigval(:,:)*ry_to_eV   !*2.0d0
        focc(:,:) = focc(:,:)/2.0d0
 
!--------------------------------------------------------------
       ! 
       do ik = 1, nk
        !
        do jbnd = 1, nbnd
         !
         do ibnd = 1, nbnd
          read (21,*) i, j, k, num1, num2, num3, num4, num5, num6
          !
          if ( ibnd .eq. jbnd ) then
           TDMat_ib(ik,ibnd,1) = CMPLX( num1, num2) * 2.0D0
           TDMat_ib(ik,ibnd,2) = CMPLX( num3, num4) * 2.0D0
           TDMat_ib(ik,ibnd,3) = CMPLX( num5, num6) * 2.0D0
           !-----------------------------------------------------------
           ! in units of velocity (p/m)
           ! in units of bohr-1, factor  of 2.0D0 comes from m_e = 1/2
           !  TDMAT  elements read from QE/6.3 and EPW
           !-----------------------------------------------------------
          endif 
          if ( abs (eigval(ik,ibnd) - eigval(ik,jbnd)) & 
              .ge. 1.0D-04 ) then  
           TDMat(ik,ibnd,jbnd,1) = CMPLX( num1, num2)/ & 
                   (eigval(ik,ibnd)-eigval(ik,jbnd)) &
                   * 14.3996 / bohr_to_ang  ! conversion to bohr
           TDMat(ik,ibnd,jbnd,2) = CMPLX( num3, num4)/ &
                   (eigval(ik,ibnd)-eigval(ik,jbnd)) &
                   * 14.3996 / bohr_to_ang  ! conversion to bohr
           TDMat(ik,ibnd,jbnd,3) = CMPLX( num5, num6)/ &
                   (eigval(ik,ibnd)-eigval(ik,jbnd)) &
                   * 14.3996 / bohr_to_ang  ! conversion to bohr
          endif
          !
         enddo ! ibnd
         !
        enddo ! jbnd
        !
       enddo ! ik
       close (21)
       write (6,*) "TDM READ"
       !
      else if ( code_flag .eq. 'EPW') then
       allocate (eigval(nk,nbnd))     ! allocate eigval for EPW
       allocate (focc(nk,nbnd))       ! allocate focc for EPW
       ! 
       nkf = nk/n_fi
       ik = 0
       do i_fi = 1, n_fi
       !
       read (21,'(a100)') line
       do iq = 1, nq
        !
        read (21,'(a100)') line
        do ikf = 1, nkf
         !
         ik = ik + 1
         read (21,'(a100)') line
         read (21,'(a100)') line
         read (21,'(a100)') line
         do ibnd = 1, nbnd
          !
          do jbnd = 1, nbnd
           !
           read (21,*) i, j, num1, num2, & 
             num3, num4, num5, num6, num7, num8
           if (ibnd .eq. 1) then
             !-------------------
             !read eigenvalue and occupation
             !------------------
             eigval(ik,jbnd) = num2
             !
             if (num2 .gt. efermi) then  !-------------------  
               focc(ik,jbnd) = 0         ! For semiconducting 
             else                        ! system only
               focc(ik,jbnd) = 1         !-------------------
             endif
             !
           endif
           !
           if ( ibnd .eq. jbnd ) then
            TDMat_ib(ik,ibnd,1) = CMPLX( num3, num4) 
            TDMat_ib(ik,ibnd,2) = CMPLX( num5, num6) 
            TDMat_ib(ik,ibnd,3) = CMPLX( num7, num8) 
            ! in units of velocity (p/m)
            ! In EPW, the factor of 2.0D0 is already multiplied 
           endif
           !--------------------------------------------
           ! Ry*bohr to eV*angstrom. Than divided by eV. 
           !--------------------------------------------
           if ( abs (num1 - num2) .ge. 1.0D-04 ) then
             TDMat(ik,ibnd,jbnd,1) = CMPLX( num3, num4)/ & 
               (num2-num1) * 7.1998           ! 0.529177249 * 13.6056980659           !
             TDMat(ik,ibnd,jbnd,2) = CMPLX( num5, num6)/ &
               (num2-num1) * 7.1998           ! conversion to Angstrom
             TDMat(ik,ibnd,jbnd,3) = CMPLX( num7, num8)/ &
               (num2-num1) * 7.1998           ! conversion to Angstrom
           endif
           !
          enddo   !jbnd
          !
         enddo    !ibnd
         !
         !-------------------------------------------
         ! To deal with EPW dipole matrix formatting
         ! ignoring k+q elements (written twice)
         !-------------------------------------------
         l_a = 1 + (nbnd*nbnd+4)  
         do l_b = 1, l_a
           read (21,'(a100)') line
         enddo
         !
        enddo ! ikf
        !
       enddo ! iq
       !
       enddo ! i_fi
       !
       close (21)
       write (6,*) "TDM READ"
       write (6,*) "ENERGY & OCCUPATIONS ON COARSE GRID READ"
       !
      else 
       write (6,*) "ONLY QE AND EPW CODES SUPPORTED"
       CALL EXIT      
      endif
      !      


!-----------------------------------------------
!     VBM (COARSE GRID) 
!-----------------------------------------------
      !

      allocate (dfocc(nk,nbnd))
      do ik = 1, nk
       do ibnd = 1, nbnd
        xarg=((efermi - eigval(ik,ibnd))/sigma)
        if ( xarg .lt. -36.0d0 ) then
          dfocc(ik,ibnd) = 0.d0
          focc(ik,ibnd) = 0.d0
        elseif (xarg .gt. 36.0d0 ) then
          dfocc(ik,ibnd) = 0.d0
          focc(ik,ibnd) = 1.d0
        else
          dfocc(ik,ibnd) = 1.0d0 * &
          exp(-xarg)/(sigma*(1 + exp(-xarg))**2)
          focc(ik,ibnd) = 1.0d0/(1 + exp(-xarg))
        endif
        !write (6,*) eigval(ik,ibnd), focc(ik,ibnd)
       enddo
      enddo
!---------------------------------------------!
! 
!  READING PART DONE
!
!---------------------------------------------!
      !---------------------------------------------!
      !  Intraband term
      !---------------------------------------------!

      !
      allocate (eps_ib(3,3,n_L))
      eps_ib (:,:,:) = CMPLX(0.0d0,0.0d0)
      !
      do iL = 1, n_L
        !
        e_L = e_L_min + dE_L * (iL - 1)         ! eV
        !write (6,'("LASER ENERGY ",f10.2)') e_L
        do ik = 1, nk
         !
         do ibnd = 1, nbnd_max
          !
          IF ( focc(ik,ibnd) <   1.0D0 ) THEN
          IF ( focc(ik,ibnd) >=  1.0D-04) THEN
          !
          fac = dfocc(ik,ibnd)/(e_L**2 + (0,1)*eta*e_L) & 
              * ry_to_eV**3
          ! conversion from eV to Ry
          do ipol = 1, 3
           !
           do jpol = 1, 3
            !
            eps_ib(ipol,jpol,iL) = eps_ib(ipol,jpol,iL) &
            - (fac*TDMat_ib(ik,ibnd,ipol)*CONJG(TDMat_ib(ik,ibnd,jpol)))
            !
           enddo ! jpol
           !
          enddo ! ipol
          !
          ENDIF
          ENDIF
         enddo
         !
        enddo
        !
      enddo
      !
      ! const = 2*4*PI*e^2/(vol*nk) = 16*PI/(vol*nk) ( if tdmat [bohr]
      ! and
      ! energy [Ry])
      const = 16 * Pi / (vol * nk) ! all in Ry natural units 
      !
      eps_ib (:,:,:) = const * eps_ib (:,:,:)
      ! 
      open (21,file='OUT_Xi_intraband.dat',action='WRITE')
!      write (21,*) "# EL(eV) ReXi_xx ImXi_xx ReXi_yy ImXi_yy ReXi_zz" &
!      " ImXi_zz ReXi_xy ImXi_xy ReXi_yz ImXi_yz ReXi_zx ImXi_zx #"
      do iL = 1, n_L       ! loop over laser frequency
        e_L =  e_L_min + dE_L * (iL - 1)
        write (21,1000) e_L, &
         eps_ib(1,1,iL), eps_ib(1,2,iL), eps_ib(1,3,iL), &
         eps_ib(2,1,iL), eps_ib(2,2,iL), eps_ib(2,3,iL), &
         eps_ib(3,1,iL), eps_ib(3,2,iL), eps_ib(3,3,iL)
      enddo
      close (21)
      write (6,*) " Intraband Xi_ij matrix calculated"
      !
      !---------------------------------------------!
      !  Interband term
      !---------------------------------------------!
      !
      allocate (eps(3,3,n_L))
      eps (:,:,:) = CMPLX(0.0d0,0.0d0)
      !
      do iL = 1, n_L           ! laser frequency 
        !
        e_L = e_L_min + dE_L * (iL - 1)         ! eV
        !write (6,'("LASER ENERGY ",f10.2)') e_L
        do ik = 1, nk 
         !
         do iv = 1, nbnd_max
          IF ( focc(ik,iv) >= 1.0E-04 ) THEN
          !
          do ic = iv+1, nbnd_max
           IF ( focc(ik,ic) <  1.0D0) THEN
           !
           et = eigval(ik,ic) - eigval(ik,iv)
           ! 
           IF (abs(focc(ik,ic)-focc(ik,iv))< 1e-4) CYCLE
           !
           fac = (focc(ik,iv)-focc(ik,ic)) * atmwfc_vc(1, ik,iv)*atmwfc_vc(2, ik,ic)*(2*et) / & 
!           (et**2 - e_L**2 - (0,1)*2*delta*e_L) &
           (et**2 - e_L**2 + delta**2 - (0,1)*2*delta*e_L) &
           * ry_to_eV
           ! conversion from eV to ry
           !

           do ipol = 1, 3
            !
            do jpol = 1, 3
             !
             eps (ipol,jpol,iL) = eps (ipol,jpol,iL) + & 
             (fac*TDMat(ik,ic,iv,ipol)*CONJG(TDMat(ik,ic,iv,jpol)))
             !
            enddo ! jpol
            !
           enddo ! ipol
           !
           ENDIF
          enddo  ! iv 
          !
          ENDIF
         enddo   ! ic
         !
        enddo   ! ik
        !
      enddo  ! iL
      !
      ! const = 2*4*PI*e^2/(vol*nk) = 16*PI/(vol*nk) ( if tdmat [bohr] and
      ! energy [Ry])
      ! 
      ! e^2 = 2
      ! factor of 2 : 2 electrons per state
      !
      ! Results for epsilon match with QE
      ! In QE const = 64*PI/(vol*nk)
      !
      const = 16 * Pi / (vol * nk) ! all in Ry natural units 
      !
      eps (:,:,:) = const * eps (:,:,:)
      !
      open (21,file='OUT_Xi.dat',action='WRITE')
!      write (21,*) "# EL(eV) ReXi_xx ImXi_xx ReXi_yy ImXi_yy ReXi_zz" &
!      " ImXi_zz ReXi_xy ImXi_xy ReXi_yz ImXi_yz ReXi_zx ImXi_zx #"
      do iL = 1, n_L       ! loop over laser frequency
        e_L =  e_L_min + dE_L * (iL - 1) 
        write (21,1000) e_L, &
         eps(1,1,iL), eps(1,2,iL), eps(1,3,iL), &
         eps(2,1,iL), eps(2,2,iL), eps(2,3,iL), &
         eps(3,1,iL), eps(3,2,iL), eps(3,3,iL)
      enddo
      close (21) 
      write (6,*) " Xi_ij matrix calculated"
      1000 format(f12.4, 1X, 9(e16.8,1X,e16.8))
      !
      call cpu_time(finish)
      print '("Time = ",f10.2," minutes.")',(finish-start)/60
      !
      end program main
!---------------------------------------------!
