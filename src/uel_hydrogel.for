! **********************************************************************
! *********** Abaqus/ STANDARD USER ELEMENT SUBROUTINE (UEL) ***********
! **********************************************************************
!           Fully coupled chemo-mechanics of non-ionic hydrogels
!   the formulation uses PK-II stress based total Lagrangian framework
! with F-bar modification for fully-integrated HEX8 and QUAD4-PE element
! currently supports first-order elements for 3D, 2D-AX, and 2D-PE cases
! **********************************************************************
!                     BIBEKANANDA DATTA (C) MAY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
!  This subroutine and dependencies are shared under 3-clause BSD license
! **********************************************************************
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL HEX8 ELEMENT
!     U3                AXISYMMETRIC TRI3 ELEMENT
!     U4                AXISYMMETRIC QUAD4 ELEMENT
!     U5                PLANE STRAIN TRI3 ELEMENT
!     U6                PLANE STRAIN QUAD4 ELEMENT
! **********************************************************************
!                       LIST OF MATERIAL PROPERTIES
!
!     props(1)-props(12): universal parameters and gel properties
!
!     Rgas        = props(1)        Universal gas constant
!     theta       = props(2)        Absolute temperature (K)
!     phi0        = props(3)        Initial polymer volume fraction
!     rho         = props(4)        Density of the gel
!     Gshear      = props(5)        Shear modulus
!     Kappa       = props(6)        Bulk modulus
!     lam_L       = props(7)        Locking stretch (only for AB model)
!     Vp          = props(8)        Molar volume of polymer
!     mu0         = props(9)        Chemical potential of pure solvent
!     Vw          = props(10)       Molar volume of the solvent
!     chi         = props(11)       Flory-Huggins parameter
!     Dw          = props(12)       Diffusion coefficient of the solvent
! **********************************************************************
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = fbarFlag        flag to use F-bar modification or not
!     jprops(3)   = matID           constitutive model for elastomeric network
!     jprops(4)   = nPostVars       no of local (int pt) post-processing variables
! **********************************************************************
!          VOIGT NOTATION CONVENTION FOR STRESS/ STRAIN TENSORS
!
!       In this subroutine we adopted the following convention for
!       symmetric stress and strain tensor following Voigt notation
!       This is different than what is followed by Abaqus/ Standard
!
!          sigma11, sigma22, sigma33, sigma23, sigma13, sigma12
!       strain11, strain22, strain33, strain23, strain13, strain12
! **********************************************************************
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)               Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)     Euler-Almansi strain tensor components
!     uvar(2*nStress+1)             Polymer volume fraction (phi)
! **********************************************************************
!               VARIABLES TO BE UPDATED WITHIN THE SUBROUTINE
!
!     RHS(i,NRHS)                   Right hand side vector
!     AMATRX(i,j)                   Stiffness matrix (NDOFEL x NDOFEL)
!     SVARS(1:NSVARS)               Element state variables.  Must be updated in this routine
!     ENERGY(1:8)                   Energy(1) Kinetic Energy
!                                   Energy(2) Elastic Strain Energy
!                                   Energy(3) Creep Dissipation
!                                   Energy(4) Plastic Dissipation
!                                   Energy(5) Viscous Dissipation
!                                   Energy(6) Artificial strain energy
!                                   Energy(7) Electrostatic energy
!                                   Energy(8) Incremental work done by loads applied to the element
!     PNEWDT                        Allows user to control ABAQUS time increments.
!                                   If PNEWDT<1 then time step is abandoned and computation is restarted with
!                                   a time increment equal to PNEWDT*DTIME
!                                   If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
!
!                       VARIABLES PROVIDED FOR INFORMATION
!
!     NDOFEL                        Total # DOF for the element
!     NRHS                          Dimension variable
!     NSVARS                        Total # element state variables
!     PROPS(1:NPROPS)               User-specified properties of the element
!     NPROPS                        No. properties
!     JPROPS(1:NJPROPS)             Integer valued user specified properties for the element
!     NJPROPS                       No. integer valued properties
!     COORDS(i,N)                   ith coordinate of Nth node on element
!     MCRD                          Maximum of (# coords,minimum of (3,#DOF)) on any node
!     Uall                          Vector of DOF at the end of the increment
!     DUall                         Vector of DOF increments
!     Vel                           Vector of velocities (defined only for implicit dynamics)
!     Accn                          Vector of accelerations (defined only for implicit dynamics)
!     JTYPE                         Integer identifying element type (the number n in the Un specification in the input file)
!     TIME(1:2)                     TIME(1)   Current value of step time
!                                   TIME(2)   Total time
!     DTIME                         Time increment
!     KSTEP                         Current step number
!     KINC                          Increment number
!     JELEM                         User assigned element number in ABAQUS
!     PARAMS(1:3)                   Time increment parameters alpha, beta, gamma for implicit dynamics
!     NDLOAD                        Number of user-defined distributed loads defined for this element
!     JDLTYP(1:NDLOAD)              Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
!     ADLMAG(1:NDLOAD)              Distributed load magnitudes
!     DDLMAG(1:NDLOAD)              Increment in distributed load magnitudes
!     PREDEF(1:2,1:NPREDF,1:NNODE)  Predefined fields.
!     PREDEF(1,...)                 Value of predefined field
!     PREDEF(2,...)                 Increment in predefined field
!     PREDEF(1:2,1,k)               Value of temperature/temperature increment at kth node
!     PREDEF(1:2,2:NPREDF,k)        Value of user defined field/field increment at kth node
!     NPREDF                        Number of predefined fields
!     LFLAGS                        Load type control variable
!     MLVARX                        Dimension variable
!     MDLOAD                        Total number of distributed loads and/or fluxes defined on this element
!     PERIOD                        Time period of the current step
! **********************************************************************
! **********************************************************************

      ! make sure the relative directory is correct
      include 'global_parameters.for'     ! global parameters module
      include 'error_logging.for'         ! error/ debugging module
      include 'linear_algebra.for'        ! linear algebra module
      include 'nonlinear_solver.for'      ! nonlinear solver module
      include 'lagrange_element.for'      ! Lagrange element module
      include 'gauss_quadrature.for'      ! Guassian quadrature module
      include 'surface_integration.for'   ! surface integration module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module
      include 'hydrogel_material.for'     ! material model module
      include 'hydrogel_element.for'      ! element formulation module

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

! **********************************************************************
!     This subroutine is called by Abaqus with above arguments
!     for each user elements defined in an Abaqus model. Users are
!     responsible for programming the element tangent/ stiffness
!     matrix and residual vectors which will be then assembled and
!     solved by Abaqus after applying the boundary conditions.
! **********************************************************************

      use global_parameters
      use error_logging
      use post_processing
      use hydrogel_element


      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT


      character(len=2)      :: analysis
      character(len=8)      :: abqProcedure
      logical               :: nlgeom
      integer               :: nInt, nPostVars
      integer               :: nDim, nStress
      integer               :: uDOF, uDOFEL, mDOF, mDOFEL


      ! initialize primary output variable to be zero
      amatrx        = zero
      rhs(:,nrhs)   = zero
      energy        = zero

      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1) .eq. 72) .or. (lflags(1) .eq. 73)) then
        abqProcedure = 'COUPLED'
      else
        call msg%ferror(flag=error, src='UEL',
     &            msg='Incorrect Abaqus procedure.', ia=lflags(1))
        call xit
      end if

      ! check if the procedure is linear or nonlinear
      ! it does not matter except for post processing
      if (lflags(2) .eq. 0) then
        NLGEOM = .false.
      else if (lflags(2) .eq. 1) then
        NLGEOM = .true.
      end if


      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4) .eq. 1) then
        call msg%ferror(flag=error, src='UEL',
     &        msg='The step should be a GENERAL step.', ia=lflags(4))
        call xit
      end if

      ! assign parameter specific to analysis and element types
      ! define different element parameters
      if ( (jtype .eq. 1) .or. (jtype .eq. 2) ) then
        nDim      = 3
        analysis  = '3D'          ! three-dimensional analysis
        nStress   = 6
      else if ( (jtype .eq. 3) .or. jtype .eq. 4) then
        nDim      = 2
        analysis  = 'AX'          ! plane axisymmetric analysis
        nStress   = 4
      else if ( (jtype .eq. 5) .or. (jtype .eq. 6) ) then
        nDim      = 2
        analysis  = 'PE'          ! 2D plane-strain analysis
        nStress   = 3
      else if ( (jtype .eq. 7) .or. (jtype .eq. 8) ) then
        nDim      = 2
        analysis  = 'PS'          ! 2D plane-stress analysis (currently unavailable)
        nStress   = 3
      else
        call msg%ferror(error,src='uel',
     &            msg='Element type is unavailable: ', ia=jtype)
        call xit
      end if


      uDOF      = nDim
      uDOFEl    = uDOF*nNode
      mDOF      = 1
      mDOFEL    = mDOF*nNODE

      nInt      = jprops(1)
      matID     = jprops(3)
      nPostVars = jprops(4)

      ! array containing variables for post-processing
      if ( .not. allocated(globalPostVars) ) then

        allocate(globalPostVars(numElem,nInt,nPostVars))

        call msg%finfo('---------------------------------------')
        call msg%finfo('--------- ABAQUS HYDROGEL UEL ---------')
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- PROCEDURE       = ', ch=abqProcedure)
        call msg%finfo('------- ANALYSIS TYPE   = ', ch=analysis)
        call msg%finfo('---------- NLGEOM       = ', la=nlgeom)
        call msg%finfo('------- MODEL DIMENSION = ', ia=nDim)
        call msg%finfo('------- ELEMENT NODES   = ', ia=nNode)
        call msg%finfo('---------------------------------------')
        call msg%finfo('-------- INTEGRATION SCHEME -----------')
        call msg%finfo('------------ NINT  = ', ia=nInt)
        call msg%finfo('---------------------------------------')
        call msg%finfo('---------- POST-PROCESSING ------------')
        call msg%finfo('--- NO OF ELEMENTS            = ',ia=numElem)
        call msg%finfo('--- OVERLAY ELEMENT OFFSET    = ',ia=elemOffset)
        call msg%finfo('--- NO OF VARIABLES AT INT PT = ',ia=nPostVars)
        call msg%finfo('---------------------------------------')

      end if

      ! return when Abaqus performs dummy step calculation with dt = 0
      if( dtime .eq. zero ) return

      ! call the first-order polyelectrolyte gel element subroutine
      select case (jtype)
      case (1, 2, 5, 6)
      ! (1) TET4, (2) HEX8, (5) TRI3-PE, and (6) QUAD4-PE
      call hydrogel_general(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &    PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     &    DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     &    NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,
     &    PERIOD,NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,MDOF,MDOFEL)

      case (3, 4)
      ! axisymmetric subroutine for (3) TRI3-AX and (4) QUAD4-AX
      call hydrogel_axisymmetric(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     &    NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,
     &    ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,
     &    JPROPS,NJPROPS,PERIOD,NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,
     &    MDOF,MDOFEL)

      case default
        call msg%ferror(flag=error, src='UEL',
     &                  msg='Wrong element type: ', ia=jtype)
        call xit
      end select

      RETURN

      END SUBROUTINE UEL

! **********************************************************************
! ************** ABAQUS USER DATABASE FILE I/O SUBROUTINE **************
! **********************************************************************

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

! **********************************************************************
!     This Abaqus/Standard user subroutine file is used for external
!     file i/o operation. User can open, close, write to external files
!     at different stages of the simulation through this subroutine
! **********************************************************************

      use global_parameters
      use error_logging

      DIMENSION                TIME(2)

      integer, intent(in)   :: lop, lrestart, kstep, kinc
      real(wp), intent(in)  :: time, dtime

      real(wp)              :: currentTime, totalTime
      integer               :: lenJobName,lenOutDir
      character(len=256)    :: jobName, outDir
      character(len=512)    :: errFile


      ! possible LOP argument parameter values
      integer, parameter    :: startAnalysis    = 0
      integer, parameter    :: startIncrement   = 1
      integer, parameter    :: endIncrement     = 2
      integer, parameter    :: endAnalysis      = 3
      integer, parameter    :: restartAnalysis  = 4
      integer, parameter    :: startStep        = 5
      integer, parameter    :: endStep          = 6

      ! possible LRESTART argument parameter values
      integer, parameter    :: restartIgnore    = 0
      integer, parameter    :: restartWrite     = 1
      integer, parameter    :: restartOverwrite = 2


      currentTime           = time(1)
      totalTime             = time(2)

      if (LOP .eq. startAnalysis) then

        call getJobName(jobName, lenJobName)
        call getOutDir(outDir, lenOutDir)
        errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
        call msg%fopen(errfile=errFile)

      else if (LOP .eq. startIncrement) then
        return

      else if (LOP .eq. endIncrement) then
        return

      else if (LOP .eq. endAnalysis) then
        call msg%finfo('Abaqus job completed successfully.')

      else if (LOP .eq. restartAnalysis) then
        return

      else if (LOP .eq. startStep) then
        return

      else if (LOP .eq. endStep) then
        return

      end if

      RETURN

      END SUBROUTINE UEXTERNALDB

! **********************************************************************
! ************** ABAQUS USER OUTPUT VARIABLES SUBROUTINE ***************
! **********************************************************************

       SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

! **********************************************************************
!     This subroutine is called by Abaqus at each material point (int pt)
!     to obtain the user defined output variables for standard Abaqus
!     elements. We used an additional layer of standard Abaqus elements
!     with same topology (same number of nodes and int pts) on top of
!     the user elements with an offset in the numbering between the user
!     elements and standard elements. This number is defined in the
!     post_processing module and should match with Abaqus input file.
! **********************************************************************

      use global_parameters
      use post_processing

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.

      ! explicityly define the type for uvar to avoid issues
      real(wp)        :: uvar

      ! assign the stored global variables to the UVAR for Abaqus to process
      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      RETURN

      END SUBROUTINE UVARM

! **********************************************************************
! **********************************************************************