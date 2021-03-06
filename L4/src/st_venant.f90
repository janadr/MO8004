!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!
!!   A simple ocean model
!!
!!   Simple initial template version for students for labs of MO7007
!!
!!   Based on one-layer, linearized shallow water equations
!!
!!   Time scheme: Fully Explicit, Leapfrog O(dt^2) => Asselin filter
!!
!!   Space:       C-grid, finite-difference, centered O(dx^2)
!!
!!
!!       AUTHORS:  Laurent Brodeau, Summer 2013
!!                 (Edited Henrik Carlson, Autumn 2013)
!!                 (Edited Sara Berglund & Aitor Aldama Campino, Spring 2017)
!!                 (Edited Aitor Aldama Campino, Spring 2018)
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM ST_VENANT

    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! Adhémar Jean Claude Barré de Saint-Venant
    !! (August 23, 1797, Villiers-en-Bière, Seine-et-Marne – January 1886,
    !! Saint-Ouen, (Loire-et-Cher) was a mechanician and mathematician who
    !! contributed to early stress analysis and also developed the unsteady
    !! open channel flow shallow water equations, also known as the
    !! Saint-Venant equations that are a fundamental set of equations used in
    !! modern hydraulic engineering.
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    USE mod_write_nc

    IMPLICIT none

    REAL, PARAMETER :: &
                       &  g = 9.81           ! gravity acceleration        [m s^-2]

    !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !! Parameters which are given a value in the namelist:
    !! - setting a default value anyway
    !! => change these parameters in the namelist file, NOT HERE!
    !! => you can save your different namelist files...

    !! DONT CHANGE THE FOLLOWING PARAMETERS HERE!!!! DO IT IN THE NAMELIST
    REAL ::                      &
         & D = 4000.,            &  ! average depth
         & f0 = 1.e-4,           &  ! Coriolis constant (only used if l_coriolis is set to true)
         & rcfl = 0.3,           &  ! CFL criterion
         & gamma = 0.1,          &  ! Asselin coefficient
         & Lx = 5.*1e+7,         &  ! size of domain (m)
         & Ly = 5.*1e+7,         &  !    //
         & h0 = 5.,              &  ! Initial amplitude
         & Lw = 20.,             &  ! Width of the gaussian (only used if l_gaussian is set to true)
         & Tm_d = 10.,           &  ! length of the run in days
         & rfsave = 1.              ! frequency at which to save fields in hours.
    !!
    INTEGER :: &
         & Nx = 100,             &  ! number of x T-points
         & Ny = 100,             &  ! number of y T-points
         & S = 20,               &  ! Size of the sponge
         & btype = 2                ! Where the sponge should be implemented.
    !!
    LOGICAL :: l_write_uv = .TRUE.         ! should we output U and V? H is always writen!
    LOGICAL :: l_coriolis = .TRUE.         ! should we include coriolis?
    LOGICAL :: l_gaussian = .TRUE.         ! Shound the initial condition be a gaussian?
    LOGICAL :: l_solidbc = .TRUE.          ! Should we use solid bc?
    LOGICAL :: l_periodic = .TRUE.         ! Should we use periodic bc?
    LOGICAL :: l_up = .TRUE.               ! Should bc be periodic in u?
    LOGICAL :: l_vp = .TRUE.               ! Should bc be periodic in v?
    LOGICAL :: l_spongebc = .TRUE.         ! Should we use a sponge bc?
    LOGICAL :: l_write_domain = .TRUE.     ! Should we save the fields with the sponge?
    LOGICAL :: l_alpha = .TRUE.            ! Should we output alpha_u and alpha_v?
    !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    LOGICAL :: lexist

    REAL    :: dx, dy, dt, ra, rb
    REAL(4)    :: tm ! length of the run in seconds...

    ! Define matrices that will be allocated later:
    REAL, DIMENSION(:, :, :), ALLOCATABLE :: u_tmp, v_tmp, h_tmp
    !
    REAL(4), DIMENSION(:, :, :), ALLOCATABLE :: u, v, h ! Storage arrays, no need for heavy double precision => float!
    !
    REAL, DIMENSION(:, :), ALLOCATABLE :: Xdu, Xdv, Xdh, U_t, V_t, U_v, V_u, alpha_u, alpha_v ! lolo
    !!
    REAL, DIMENSION(:), ALLOCATABLE :: vx_t, vy_t, vx_u, vy_v, vtime

    !! Indices for loops
    INTEGER :: nt, nt_save, jt_save = 1, ji, jj, jf, jt, jc, nb_save, sx, sy
    REAL    :: du, dv, dh

    !! Defining namelist sections:
    NAMELIST /ngrid/ Nx, Ny, Lx, Ly, Tm_d, rfsave
    NAMELIST /nphysics/ D, l_coriolis, f0
    NAMELIST /ninitial/ l_gaussian, h0, Lw
    NAMELIST /nnumerics/ rcfl, gamma
    NAMELIST /nboundaries/ l_write_uv, l_solidbc, l_periodic, l_up, l_vp, l_spongebc, S, btype, l_write_domain, l_alpha

    !! Declarations are done, starting the program !

    !! First testing if the namelist file is present in the current directory:
    INQUIRE (FILE='namelist', EXIST=lexist)
    IF (.NOT. lexist) THEN
        PRINT *, 'ERROR: file "namelist" not found!'; STOP
    END IF
    PRINT *, 'Reading namelist..'
    !! Opening and reading the namelist:
    OPEN (UNIT=11, FILE='namelist', FORM='FORMATTED', STATUS='OLD')
    READ (11, ngrid)
    READ (11, nphysics)
    READ (11, ninitial)
    READ (11, nnumerics)
    READ (11, nboundaries)
    CLOSE (11)
    PRINT *, 'Done reading namelist.'
    !! Namelist paramters are now known!

    dx = Lx/REAL(Nx)  !delta x
    dy = Ly/REAL(Ny)  !delta y

    sx = 0
    sy = 0

    !! Redifining the size of the grid if using a sponge bc.
    IF (l_spongebc .and. (btype == 0 .or. btype == 2)) THEN
        Nx = Nx + S
        sx = S/2
    END IF
    IF (l_spongebc .and. (btype == 0 .or. btype == 1)) THEN
        Ny = Ny + S
        sy = S/2
    END IF

    PRINT *, 'Defining time-step of the model'
    !! Time stuff
    !! ==========

    tm = real(Tm_d*24.*3600., 8) ! Length of the run in seconds...

    !! Time step length is determined by grid box size,
    dt = 0.
    dt = rcfl*min(dx, dy)/sqrt(g*D)   !! a CFL-number, and phase speed (c^2 = g*D)
    PRINT *, ' Initial dt =', dt
    IF (dt > 10.) THEN
        IF (MOD(int(dt), 10) /= 0) THEN
            dt = REAL(int(dt) - MOD(int(dt), 10))  !! Want a multiple of 10 seconds:
        END IF
    END IF
    dt = REAL(int(dt))  ! in secons without decimals!
    PRINT *, ' *** Time step to use in seconds:', dt

    nt = tm/dt; ! How many time steps to go:
    PRINT *, 'Number of time steps to go =', nt

    nt_save = (rfsave*3600.)/dt
    PRINT *, ' Will save every', nt_save, ' time steps...'

    nb_save = INT(nt/nt_save)
    IF (MOD(nt, nt_save) /= 0) THEN
        PRINT *, 'Extending run length so it matches the save frequency!'
        nb_save = nb_save + 1
        nt = nt_save*nb_save
    END IF
    nb_save = nb_save + 1 ! include initial condition
    PRINT *, ' Actual number of time steps to go: nt =', nt
    !PRINT *, ' Will save', nb_save, 'time shots!'
    PRINT *, '  => new run length in days is ', nt*dt/(24.*3600.)
    PRINT *, ''

    !! ====================================================
    !! === Start ===
    !! ====================================================

    PRINT *, ''
    PRINT *, ' =============== Shallow water model ============ '
    PRINT *, '  Lx, Nx, dx = ', Lx, Nx, dx
    PRINT *, '  Ly, Ny, dy = ', Ly, Ny, dy
    PRINT *, '  Mean depth = ', D
    PRINT *, '  dt = ', dt
    PRINT *, '  Number of numerical time steps :', nt
    PRINT *, '  Number of time snapshots to save   :', nb_save
    PRINT *, ' =============================================== '
    PRINT *, ' ============== Settings of model ============== '
    IF (l_coriolis) PRINT *, '  Coriolis is on'
    IF (l_periodic) THEN
        IF (l_vp) PRINT *, '  Periodic boundary condition in N-S on'
        IF (l_up) PRINT *, '  Periodic boundary condition in E-W on'
    END IF
    IF (l_spongebc) PRINT *, '  Sponge boundary condition on.'
    IF (l_solidbc) PRINT *, '  Solid boundary condition on'
    PRINT *, ' ================================================ '
    PRINT *, ''

    !! Allocating arrays:
    ALLOCATE (h(Nx, Ny, nb_save), vx_u(0:Nx), vy_v(0:Ny), vx_t(Nx), vy_t(Ny), &
         &     vtime(nb_save), u_tmp(0:Nx, Ny, -1:1), v_tmp(Nx, 0:Ny, -1:1), h_tmp(Nx, Ny, -1:1),   &
         &     Xdu(0:Nx, Ny), Xdv(Nx, 0:Ny), Xdh(Nx, Ny),    &
         &     U_t(Nx, Ny), V_t(Nx, Ny), U_v(1:Nx, 0:Ny), V_u(0:Nx, 1:Ny))

    IF (l_write_uv) ALLOCATE (u(0:Nx, Ny, nb_save), v(Nx, 0:Ny, nb_save))

    !! Building coordinates and calendar vectors:
    vx_u(:) = (/((ji - sx)*dx, ji=0, Nx)/) ! x coordinates at U-points, starts for i=0 ! => size Nx+1
    vy_v(:) = (/((jj - sy)*dy, jj=0, Ny)/) ! y coordinates at V-points, starts for j=0 ! => size Ny+1

    vx_t(:) = vx_u(0:Nx - 1) + dx/2.; ! x coordinates at T-points...
    vy_t(:) = vy_v(0:Ny - 1) + dy/2.; ! y coordinates at T-points...

    vtime(:) = (/((jt - 1)*rfsave, jt=1, nb_save)/) ! in hours !!!

    !! Building sponge.

    ALLOCATE (alpha_u(0:Nx, Ny))
    ALLOCATE (alpha_v(Nx, 0:Ny))

    IF (l_spongebc) THEN
        CALL SPONGE_BC(alpha_u, alpha_v)
    END IF

    PRINT *, 'Defining Initial condition'
    !! ====================================================
    !! === Initial condition ===
    !! ====================================================
    IF (l_gaussian) THEN
        PRINT *, 'Initial condition is set as Guassian'
        !! Gaussian inside:

        ! Width of the gaussian
        Lw = min(Lx, Ly)/Lw
        PRINT *, 'Width of gaussian is set to:', Lw
        PRINT *, 'Hight of gaussian is set to:', h0

        DO jj = 1, Ny
            DO ji = 1, Nx
                h(ji, jj, 1) = h0*exp(-((vx_t(ji) - 0.5*Lx)/Lw)**2.-((vy_t(jj) - 0.5*Ly)/Lw)**2.)
            END DO
        END DO
    END IF

    ! Saving initial condition:
    CALL WRITE_NC_2D(vx_t, vy_t, h(:, :, 1), 'H_init_cond.nc', 'h0')

    h_tmp(:, :, -1) = h(:, :, 1)

    jt_save = 1; jt = 0
    PRINT *, 'Save fields at time', REAL(vtime(jt_save), 4), ' time step =', jt, ' of', nt
    jt_save = jt_save + 1

    PRINT *, 'Done saving initial condition.'
    !! =========================================================================
    !! === Solving with Euler forward on a C-grid for first time step (jt=1) ===
    !! =========================================================================

    jt = jt + 1

    Xdu = 0.; Xdv = 0.; Xdh = 0.

    ! Coriolis
    ! ========

    IF (l_coriolis) THEN
        ! Define V_u, v in a u point
        ! Define U_v, u in a v point
        

    END IF

    ! Height gradient
    ! ===============
    Xdu(1:Nx - 1, :) = Xdu(1:Nx - 1, :) - g/dx*(h_tmp(2:Nx, :, -1) - h_tmp(1:Nx - 1, :, -1))
    Xdv(:, 1:Ny - 1) = Xdv(:, 1:Ny - 1) - g/dy*(h_tmp(:, 2:Ny, -1) - h_tmp(:, 1:Ny - 1, -1))

    IF (l_up .and. l_periodic) THEN
        Xdu(0, :) = Xdu(0, :) - g/dx*(h_tmp(1, :, -1) - h_tmp(Nx, :, -1))
        Xdu(Nx, :) = Xdu(Nx, :) - g/dx*(h_tmp(1, :, -1) - h_tmp(Nx, :, -1))
    END IF
    IF (l_vp .and. l_periodic) THEN
        Xdv(:, 0) = Xdv(:, 0) - g/dx*(h_tmp(:, 1, -1) - h_tmp(:, Ny, -1))
        Xdv(:, Ny) = Xdv(:, Ny) - g/dx*(h_tmp(:, 1, -1) - h_tmp(:, Ny, -1))
    END IF

    ! Adding Coriolis
    ! ===============
    IF (l_coriolis) THEN
        ! Add Coriolis here

    END IF

    ! Convergence/divergence
    ! ======================
    Xdh(:, :) = Xdh(:, :) - D*((u_tmp(1:Nx, :, -1) - u_tmp(0:Nx - 1, :, -1))/dx + (v_tmp(:, 1:Ny, -1) - v_tmp(:, 0:Ny - 1, -1))/dy)

    ! Time step
    ! =========
    u_tmp(:, :, 0) = u_tmp(:, :, -1) + Xdu(:, :)*dt
    v_tmp(:, :, 0) = v_tmp(:, :, -1) + Xdv(:, :)*dt
    h_tmp(:, :, 0) = h_tmp(:, :, -1) + Xdh(:, :)*dt

    IF (l_solidbc) THEN
        CALL APPLY_SOLID_BC(jt_step=-1)
    END IF
    IF (l_periodic) THEN
        CALL APPLY_PERIODIC_BC(jt_step=-1)
    END IF
    IF (l_spongebc) THEN
        CALL APPLY_SPONGE_BC(jt_step=-1)
    END IF

    PRINT *, 'Done with first Euler step.                  time step =', jt

    !****************  M A I N   T I M E   L O O P ******************************

    ! ====================================================
    ! === Solving with Leap frog on a C-grid ===
    ! ====================================================

    ! Time loops
    DO jt = 2, nt

        Xdu = 0.; Xdv = 0.; Xdh = 0.

        ! Coriolis
        ! =============

        IF (l_coriolis) THEN
            ! Define V_u, v in a u point
            ! Define U_v, u in a v point
        END IF

        ! Height gradient
        ! ===============
        Xdu(1:Nx - 1, :) = Xdu(1:Nx - 1, :) - g/dx*(h_tmp(2:Nx, :, 0) - h_tmp(1:Nx - 1, :, 0))
        Xdv(:, 1:Ny - 1) = Xdv(:, 1:Ny - 1) - g/dy*(h_tmp(:, 2:Ny, 0) - h_tmp(:, 1:Ny - 1, 0))

        IF (l_up .and. l_periodic) THEN
            Xdu(0, :) = Xdu(0, :) - g/dx*(h_tmp(1, :, 0) - h_tmp(Nx, :, 0))
            Xdu(Nx, :) = Xdu(Nx, :) - g/dx*(h_tmp(1, :, 0) - h_tmp(Nx, :, 0))
        END IF
        IF (l_vp .and. l_periodic) THEN
            Xdv(:, 0) = Xdv(:, 0) - g/dx*(h_tmp(:, 1, 0) - h_tmp(:, Ny, 0))
            Xdv(:, Ny) = Xdv(:, Ny) - g/dx*(h_tmp(:, 1, 0) - h_tmp(:, Ny, 0))
        END IF

        ! Adding Coriolis
        ! ===============
        IF (l_coriolis) THEN
            ! Add Coriolis here
        END IF

        ! Convergence/divergence
        ! ===============
        Xdh(:, :) = Xdh(:, :) - D*((u_tmp(1:Nx, :, 0) - u_tmp(0:Nx - 1, :, 0))/dx + (v_tmp(:, 1:Ny, 0) - v_tmp(:, 0:Ny - 1, 0))/dy)

        ! Time step
        ! =========
        u_tmp(:, :, 1) = u_tmp(:, :, -1) + Xdu(:, :)*2.*dt
        v_tmp(:, :, 1) = v_tmp(:, :, -1) + Xdv(:, :)*2.*dt
        h_tmp(:, :, 1) = h_tmp(:, :, -1) + Xdh(:, :)*2.*dt

        !! Applying Asselin filter:
        u_tmp(:, :, 0) = (1.-2.*gamma)*u_tmp(:, :, 0) + gamma*(u_tmp(:, :, 1) + u_tmp(:, :, -1))
        v_tmp(:, :, 0) = (1.-2.*gamma)*v_tmp(:, :, 0) + gamma*(v_tmp(:, :, 1) + v_tmp(:, :, -1))
        h_tmp(:, :, 0) = (1.-2.*gamma)*h_tmp(:, :, 0) + gamma*(h_tmp(:, :, 1) + h_tmp(:, :, -1))

        ! === Boundaries ===
        IF (l_solidbc) THEN
            CALL APPLY_SOLID_BC()
        END IF
        IF (l_periodic) THEN
            CALL APPLY_PERIODIC_BC()
        END IF
        IF (l_spongebc) THEN
            CALL APPLY_SPONGE_BC()
        END IF

        ! Re-arrange matrices
        u_tmp(:, :, -1) = u_tmp(:, :, 0)
        v_tmp(:, :, -1) = v_tmp(:, :, 0)
        h_tmp(:, :, -1) = h_tmp(:, :, 0)
        u_tmp(:, :, 0) = u_tmp(:, :, 1)
        v_tmp(:, :, 0) = v_tmp(:, :, 1)
        h_tmp(:, :, 0) = h_tmp(:, :, 1)
        u_tmp(:, :, 1) = 0.
        v_tmp(:, :, 1) = 0.
        h_tmp(:, :, 1) = 0.

        IF (MOD(jt, nt_save) == 0) THEN

            ! Save time t (not t+1) to u,v,h matrices
            IF (l_write_uv) THEN
                u(:, :, jt_save) = REAL(u_tmp(:, :, 0), 4)
                v(:, :, jt_save) = REAL(v_tmp(:, :, 0), 4)
            END IF
            h(:, :, jt_save) = REAL(h_tmp(:, :, 0), 4)

            PRINT *, 'Save fields at time', REAL(vtime(jt_save), 4), ' time step =', jt, ' of', nt

            jt_save = jt_save + 1
        END IF

    END DO   ! END of time LOOP!!!

    PRINT *, 'End of Program.'

    ! Write data to files
    !====================

    ! Choose to write the whole domain (when using a sponge bc) or only the fields.
    IF (l_write_domain .and. l_spongebc) THEN

        IF (btype == 2 .or. btype == 0) THEN
            sx = S/2
        END IF
        IF (btype == 1 .or. btype == 0) THEN
            sy = S/2
        END IF

        CALL WRITE_NC(vx_t(1 + sx:Nx - sx), vy_t(1 + sy:Ny - sy), vtime, h(1 + sx:Nx - sx, 1 + sy:Ny - sy, :), 'h_test.nc', 'h')

        IF (l_write_uv) THEN
            CALL WRITE_NC(vx_u(0 + sx:Nx - sx), vy_t(1 + sy:Ny - sy), vtime, u(0 + sx:Nx - sx, 1 + sy:Ny - sy, :), 'u_test.nc', 'u')
            CALL WRITE_NC(vx_t(1 + sx:Nx - sx), vy_v(0 + sy:Ny - sy), vtime, v(1 + sx:Nx - sx, 0 + sy:Ny - sy, :), 'v_test.nc', 'v')
        END IF

    ELSE

        CALL WRITE_NC(vx_t(1:Nx), vy_t(1:Ny), vtime, h(1:Nx, 1:Ny, :), 'h_test.nc', 'h')

        IF (l_write_uv) THEN
            CALL WRITE_NC(vx_u(0:Nx), vy_t(1:Ny), vtime, u(0:Nx, 1:Ny, :), 'u_test.nc', 'u')
            CALL WRITE_NC(vx_t(1:Nx), vy_v(0:Ny), vtime, v(1:Nx, 0:Ny, :), 'v_test.nc', 'v')
        END IF

    END IF

CONTAINS

    SUBROUTINE APPLY_SOLID_BC(jt_step)

        INTEGER, OPTIONAL, INTENT(in) :: jt_step   !: at what time step are we in *_tmp files?
        !!                                            default => jt_step=2
        !!                                            Only 2 or 1 !

        INTEGER :: jts, jtn

        jts = 0
        IF (present(jt_step)) jts = jt_step

        jtn = jts + 1  !: next time step

        !! Western boundary:
        !! =================

        u_tmp(0, :, jtn) = 0.   ! Solid wall
        u_tmp(0, :, jts) = 0.   ! Solid wall

        !! Eastern boundary:
        !! =================

        u_tmp(Nx, :, jtn) = 0.   ! Solid wall
        u_tmp(Nx, :, jts) = 0.   ! Solid wall

        ! === Southern boundary:
        ! ======================

        v_tmp(:, 0, jtn) = 0.   ! Solid wall
        v_tmp(:, 0, jts) = 0.   ! Solid wall

        ! === Northern boundary:
        ! ======================

        v_tmp(:, Ny, jtn) = 0.   ! Solid wall
        v_tmp(:, Ny, jts) = 0.   ! Solid wall

    END SUBROUTINE APPLY_SOLID_BC

    SUBROUTINE SPONGE_BC(alpha_u, alpha_v)

        REAL, DIMENSION(0:Nx, Ny), INTENT(inout) :: alpha_u
        REAL, DIMENSION(Nx, 0:Ny), INTENT(inout) :: alpha_v

        INTEGER :: Sh, iii

        IF (mod(S, 2) /= 0) THEN
            PRINT *, 'S is not an even number'
            STOP
        ELSE
            Sh = S/2
        END IF

        DO iii = 0, Sh

            IF (btype == 0) THEN
                alpha_u(iii:Nx - iii, 1 + iii:Ny - iii) = FLOAT(Sh - iii)/Sh
                alpha_v(iii + 1:Nx - iii, iii:Ny - iii) = FLOAT(Sh - iii)/Sh
            ELSEIF (btype == 1) THEN
                alpha_u(:, 1 + iii:Ny - iii) = FLOAT(Sh - iii)/Sh
                alpha_v(:, iii:Ny - iii) = FLOAT(Sh - iii)/Sh
            ELSEIF (btype == 2) THEN
                alpha_u(iii:Nx - iii, :) = FLOAT(Sh - iii)/Sh
                alpha_v(iii + 1:Nx - iii, :) = FLOAT(Sh - iii)/Sh
            END IF

        END DO

        IF (l_alpha) THEN
            ! Write alpha_u and alpha_v
            CALL WRITE_NC_2D(vx_u(0:Nx), vy_t(1:Ny), REAL(alpha_u(0:Nx, 1:Ny), 4), 'alpha_u.nc', 'alpha_u')
            CALL WRITE_NC_2D(vx_t(1:Nx), vy_v(0:Ny), REAL(alpha_v(1:Nx, 0:Ny), 4), 'alpha_v.nc', 'alpha_v')
        END IF

    END SUBROUTINE SPONGE_BC

    SUBROUTINE APPLY_SPONGE_BC(jt_step)

        INTEGER, OPTIONAL, INTENT(in) :: jt_step

        INTEGER :: jts, jtn

        jts = 0
        IF (present(jt_step)) jts = jt_step
        jtn = jts + 1

        v_tmp(:, :, jtn) = (1.-alpha_v(:, :))*v_tmp(:, :, jtn)
        u_tmp(:, :, jtn) = (1.-alpha_u(:, :))*u_tmp(:, :, jtn)

    END SUBROUTINE APPLY_SPONGE_BC

    SUBROUTINE APPLY_PERIODIC_BC(jt_step)

        INTEGER, OPTIONAL, INTENT(in) :: jt_step ! At what time step are we in *_tmp files?

        INTEGER :: jts, jtn

        jts = 0
        IF (present(jt_step)) jts = jt_step
        jtn = jts + 1

        IF (l_up) THEN
            u_tmp(Nx, :, jtn) = u_tmp(0, :, jtn)
        END IF

        IF (l_vp) THEN
            v_tmp(:, Ny, jtn) = v_tmp(:, 0, jtn)
        END IF

    END SUBROUTINE APPLY_PERIODIC_BC

END PROGRAM ST_VENANT
