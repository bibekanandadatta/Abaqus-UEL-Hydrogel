      include 'uel_hydrogel_dev.for'

! **********************************************************************
! **********************************************************************

!     using ifort compiler compile the code from intel oneAPI terminal:
!     ifort /Qmkl -o gel test_hydrogel.for
!     this command will generate an executable call gel
!     execute the pegel on the terminal using: .\gel

!     the purpose of this code is to compare the the element tangent
!     matrix for different element types during development
!     currently it tests JTYPE=4 (QUAD4-PE)
!     it prints out the results in a .dat file based on the JTYPE

! **********************************************************************
! **********************************************************************

      PROGRAM HYDROGEL_DEV

      use global_parameters, only: wp, zero, one, two, half
      implicit none

      real(wp), allocatable :: RHS(:,:), AMATRX(:,:), SVARS(:)
      real(wp), allocatable :: PROPS(:), COORDS(:,:), DUall(:,:)
      real(wp), allocatable :: Uall(:), Vel(:), Accn(:)
      real(wp), allocatable :: ADLMAG(:,:), PREDEF(:,:,:), DDLMAG(:,:)

      integer, allocatable  :: JPROPS(:), JDLTYP(:,:)

      real(wp)  :: ENERGY(8), PNEWDT, TIME(2), DTIME, PARAMS(3), PERIOD

      integer   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE
      integer   :: KSTEP, KINC, JELEM, NDLOAD, NPREDF
      integer   :: MLVARX, MDLOAD, NJPROPS, LFLAGS(5), nlSdv, ngSdv


      integer   :: nDim, nStress, nOrder
      integer   :: nInt, fbarFlag, matID, nPostVars
      real(wp)  :: Rgas, theta, phi0, rho, Gshear, Kappa, lam_L, Vp
      real(wp)  :: mu0, Vw, chi, Dw, initMU

      real(wp)  :: elemLen, lam_1_App, lam_2_App, u_1_App, u_2_App
      real(wp)  :: mu_App_Fac

      integer   :: i, j, k, l

      ! additional vectors and matrices for numerical derivatives
      real(wp), allocatable :: Uall_1(:), Uall_2(:)
      real(wp), allocatable :: RHS_1(:,:), AMATRX_1(:,:)
      real(wp), allocatable :: RHS_2(:,:), AMATRX_2(:,:)
      real(wp), allocatable :: AMATRX_h(:,:)
      real(wp)              :: delta_h

      integer, parameter    :: case     = 1
      integer, parameter    :: fileUnit = 15
      character(len=256)    :: fileName



      !! element type, element no, dimension, etc.
      JELEM     = 1           ! no of elements in patch test
      JTYPE     = 4           ! type of element
      elemLen   = 1.0e-3_wp   ! element side length

      ! nInt : no of volume integration point, nPostVars: no of post-processed variables
      if (JTYPE .eq. 4) then
        nDim      = 2               ! spatial dimension of the problem
        nStress   = 3               ! independent stress quantities
        nNode     = 4               ! no of nodes of the element
        nOrder    = 1               ! polynomial order of the lagrangian element
        nInt      = 4               ! no of integration points
        nPostVars = 2*nStress + 1
        fileName  = 'hydrogel_debug_QUAD4PE.dat'
      else
        write(*,'(A)')  'Element unavailable for debugging: ', jtype
        stop
      endif

      fbarFlag    = 0
      matID       = 1

      !! open the debugging file
      open(unit=fileUnit, file=fileName, status='unknown')

      !! Abaqus variables
      NRHS    = 1               ! right hand side residual column dimension
      NSVARS  = nInt            ! no of state variables
      NPROPS  = 12              ! total no of material properties
      MCRD    = nDim            ! no of dimension (ABAQUS does it differently)
      NJPROPS = 4               ! no of integer properties
      NPREDF  = 1               ! no of predefined field
      NDOFEL  = nNode*(nDim+1)  ! no of total degrees of freedom
      MLVARX  = NDOFEL          ! maximum variables (assuming same as NDOFEL)

      ! no distributed load
      NDLOAD  = 1
      MDLOAD  = 1

      allocate( RHS(MLVARX,NRHS), AMATRX(NDOFEL,NDOFEL), SVARS(NSVARS),
     &      PROPS(NPROPS), COORDS(MCRD,NNODE), DUALL(MLVARX,1),
     &      UALL(NDOFEL), Vel(NDOFEL), Accn(NDOFEL),
     &      JDLTYP(MDLOAD,1), ADLMAG(MDLOAD,1), DDLMAG(MDLOAD,1),
     &      PREDEF(2,NPREDF,NNODE), JPROPS(NJPROPS) )


      ! initializing the output of UEL subroutine
      RHS     = zero            ! residual vector
      AMATRX  = zero            ! stiffness matrix
      SVARS   = zero            ! array containing state variables
      ENERGY  = zero            ! array containing energy
      PNEWDT  = zero            ! time stepping flag

      ! initialize some other not-so-necessary (for statics) input to UEL
      LFLAGS = [72, 1, 0, 0, 0] ! Abaqus step flag
      VEL     = zero            ! velocity
      ACCN    = zero            ! acceleration
      PARAMS  = zero            ! time parameters (irrelevant now)
      PREDEF  = zero            ! no predefined field

      KSTEP   = 1
      KINC    = 1
      PERIOD  = 1.0e-3_wp
      TIME    = 1.0e-3_wp
      DTIME   = 1.0e-3_wp


      if (JTYPE .EQ. 4) then
        COORDS(1,:) = [0, 1, 1, 0]
        COORDS(2,:) = [0, 0, 1, 1]
      end if

      COORDS  = COORDS*elemLen


      !! material properties
      Rgas        = 8.3145_wp         ! Universal gas constant
      theta       = 298.0_wp          ! Absolute temperature (K)
      phi0        = 0.99_wp           ! Initial polymer volume fraction
      rho         = 1100.0_wp         ! Density of the gel
      Gshear      = 100.0e3_wp        ! Shear modulus
      Kappa       = 50.0_wp*Gshear    ! Bulk modulus
      lam_L       = 0.0_wp            ! Locking stretch (only for AB model)
      Vp          = 8.92e-3_wp        ! Molar volume of polymer
      mu0         = 0.0_wp            ! Chemical potential of pure solvent
      Vw          = 1.8e-5_wp         ! Molar volume of the solvent
      chi         = 0.485_wp          ! Flory-Huggins parameter
      Dw          = 5.0e-9_wp         ! Diffusion coefficient of the solvent

      matID       = 1

      ! allocate the material and element properties
      PROPS(1:NPROPS)   = [ Rgas, theta, phi0, rho, Gshear, Kappa,
     &                      lam_L, Vp, mu0, Vw, chi, Dw ]
      JPROPS(1:NJPROPS) = [ nInt, fbarFlag, matID, nPostVars ]

      initMU      = mu0
     &              + Rgas*theta*(phi0 + log(1-phi0) + chi*phi0**two)
     &              - Kappa*Vw*log(phi0)
     &              + half*Kappa*Vw*(log(phi0))**two



      !! define the "simulated" nodal solutions (or boundary )
      Uall(1:nDOFEL) = [zero,   zero,   initMU,
     &                  zero,   zero,   initMU,
     &                  zero,   zero,   initMU,
     &                  zero,   zero,   initMU]



      !! call the UEL subroutine
      call UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)



      !! numerical tangent matrix calculation
      allocate( UALL_1(NDOFEL), RHS_1(MLVARX,NRHS),
     &          UALL_2(NDOFEL), RHS_2(MLVARX,NRHS),
     &          AMATRX_1(NDOFEL,NDOFEL), AMATRX_2(NDOFEL,NDOFEL),
     &          AMATRX_h(NDOFEL,NDOFEL) )


      !! perturb the degrees of freedom and calculate the tangent matrix
      delta_h   = sqrt(epsilon(one))

      do i = 1, ndofel
        do j = 1, ndofel

          ! perturb in +ve direction
          RHS_1     = zero
          AMATRX_1  = zero
          SVARS     = zero
          Uall_1    = Uall
          Uall_1(j) = Uall_1(j) + delta_h

          call UEL(RHS_1,AMATRX_1,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &      PROPS,NPROPS,COORDS,MCRD,NNODE,Uall_1,DUall,Vel,Accn,JTYPE,
     &      TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &      NJPROPS,PERIOD)

          ! perturb in -ve direction
          RHS_2     = zero
          AMATRX_2  = zero
          SVARS     = zero
          Uall_2    = Uall
          Uall_2(j) = Uall_2(j) - delta_h

          call UEL(RHS_2,AMATRX_2,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &      PROPS,NPROPS,COORDS,MCRD,NNODE,Uall_2,DUall,Vel,Accn,JTYPE,
     &      TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &      NJPROPS,PERIOD)

          AMATRX_h(i,j) = - (RHS_1(i,1) - RHS_2(i,1))/(2*delta_h)
        end do
      end do


      write(fileUnit,'(A)') 'RHS: '
      do i = 1,NDOFEL
        write(fileUnit,*) RHS(i,1)
      enddo

      write(fileUnit,'(A)') 'AMATRX: '
      do i = 1, NDOFEL
        write(fileUnit, '(*(F24.14))') AMATRX(i,:)
      end do

      write(fileUnit,'(A)') 'AMATRX_h: '
      do i = 1, NDOFEL
        write(fileUnit, '(*(F24.14))') AMATRX_h(i,:)
      end do


      !! close the debugging file
      close(unit=fileUnit)


      END PROGRAM HYDROGEL_DEV

! **********************************************************************

      ! this subroutine emulates the XIT subroutine in ABAQUS
      SUBROUTINE XIT
        stop
      END SUBROUTINE XIT

! **********************************************************************

      SUBROUTINE PRINT_MAT( Kuu, Kum, Kmu, Kmm, Ru, Rm, Amatrx, RHS,
     &                      nDim, nNode, uDOFEL, mDOFEL, nDOFEl )

      use global_parameters, only: wp

      implicit none

      integer,    intent(in)  :: nDim, nNode, uDOFEL, mDOFEL, nDOFEl
      real(wp),   intent(in)  :: Kuu(uDOFEL,UDOFEL)
      real(wp),   intent(in)  :: Kum(uDOFEL,mDOFEL)
      real(wp),   intent(in)  :: Kmu(mDOFEL,uDOFEL), Kmm(mDOFEL,mDOFEL)
      real(wp),   intent(in)  :: Amatrx(nDOFEL,nDOFEL)
      real(wp),   intent(in)  :: Ru(uDOFEL,1), Rm(mDOFEL,1)
      real(wp),   intent(in)  :: RHS(nDOFEL,1)

      integer                 :: i

      integer, parameter      :: fileUnit = 15


      write(fileUnit,*) 'Kuu: '
      do i = 1,nDim*nNode
        write(fileUnit,'(100(F16.6,2x))') Kuu(i,:)
      enddo

      write(fileUnit,*) 'Kum: '
      do i = 1,nDim*nNode
        write(fileUnit,'(100(F18.14,2x))') Kum(i,:)
      enddo

      write(fileUnit,*) 'Kmu: '
      do i = 1,nNode
        write(fileUnit,'(100(F18.12,1x))') Kmu(i,:)
      enddo

      write(fileUnit,*) 'Kmm: '
      do i = 1,nNode
        write(fileUnit,'(100(F18.14,2x))') Kmm(i,:)
      enddo

      write(fileUnit,'(A)') 'Ru: '
      do i = 1,UDOFEL
        write(fileUnit,*) Ru(i,1)
      enddo


      write(fileUnit,'(A)') 'Rm: '
      do i = 1,MDOFEL
        write(fileUnit,*) Rm(i,1)
      enddo


      write(fileUnit,'(A)') 'AMATRX: '
      do i = 1, NDOFEL
        write(fileUnit,'(100(g0,2x))') AMATRX(i,:)
      end do


      write(fileUnit,'(A)') 'RHS: '
      do i = 1,NDOFEL
        write(fileUnit,*) RHS(i,1)
      enddo

      END SUBROUTINE PRINT_MAT

! **********************************************************************