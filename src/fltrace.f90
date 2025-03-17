!*******************************************************************************
PROGRAM fltrace
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code.
!*******************************************************************************
    USE shared_data
    USE current_tools
    USE grid

    IMPLICIT NONE

    CHARACTER(LEN=64):: input_value
    CHARACTER(LEN=64):: parameter_filename
    CHARACTER(LEN=5):: snap_id


    call get_command_argument(1, input_value)
    read(unit=input_value,fmt=*) snap

    write (snap_id,'(I5.5)') snap
    parameter_filename = trim('./fl_data/flparameters'//trim(snap_id)//'.txt')

    !Import parameters from text file (saved out by python in the fltrace directory)
    !Need to update this so it can do different numbers
    open(1,file= parameter_filename)
    read(1, *) flparameters
    close(1)

    !##########################################

    ! Put some of the major variables in here - things that can be changed occasionally but not in a series of runs

    run = int(flparameters(0))
    nx = int(flparameters(1))
    ny = int(flparameters(2))
    nz = int(flparameters(3))
    x0 = flparameters(4)
    x1 = flparameters(5)
    y0 = flparameters(6)
    y1 = flparameters(7)
    z0 = flparameters(8)
    z1 = flparameters(9)

    nx_out = int(flparameters(18))  !Number of CELLS in the output emissivity matrix thing
    ny_out = int(flparameters(19))
    nz_out = int(flparameters(20))

    machine_flag = int(flparameters(17))
    print_flag = int(flparameters(21))
    save_all = int(flparameters(22))

    !DATA ROOT HERE. NEED TO READ IN AS STRING REALLY

    data_root = './tmp/'

    snap = int(flparameters(10))
    nstarts = int(flparameters(11))

    max_line_length = int(flparameters(13))
    ds_factor =  flparameters(14)
    weakness_limit =  flparameters(15)

    data_source = int(flparameters(16))

    data_root = trim(data_root)

    call establish_grid()  !Establish the grid and read in the magnetic field

    !call calculate_current()  !Don't need to do this now as the current is read in directly

    call establish_starts()  !Import the location of the start points

    call bfield_boundary() !Populate ghost points so the interpolator works

    call integrate_fieldlines()

    call export_emissivity()

    contains

    SUBROUTINE establish_starts

    !Read in the start points from the provided 'starts.txt'
    IMPLICIT NONE

    INTEGER:: i
    REAL(num), DIMENSION(:):: starts_import(0:nstarts*3-1)
    CHARACTER(LEN=64):: starts_filename
    CHARACTER(LEN=5):: snap_id


    ALLOCATE(starts(0:nstarts-1,0:2))

    write (snap_id,'(I5.5)') snap

    starts_filename = trim('./fl_data/starts'//trim(snap_id)//'.txt')

    open(1,file= starts_filename)
    read(1, *) starts_import
    close(1)

    do i = 0, nstarts-1
        starts(i,:) = starts_import(3*i:3*i+2)
    end do

    END SUBROUTINE establish_starts

    SUBROUTINE interpolate_bfield(x_point, y_point, z_point)
    !Finds the magnetic field vector at a given point, taking advantage of the staggered grid.
    IMPLICIT NONE

    REAL(num):: x_point, y_point, z_point !Coordinates of point to interpolate
    REAL(num):: xp, yp, zp !Coordinates in the respective dimentions
    INTEGER:: xi, yi, zi !Cell indices in each dimension
    REAL(num):: xf, yf, zf !Distance 'up' each cell

    b1 = 0.0_num
    !Establish ratios
    xp = nx*(x_point - x0)/(x1 - x0)
    yp = ny*(y_point - y0)/(y1 - y0)
    zp = nz*(z_point - z0)/(z1 - z0)

    !Interpolate bx
    xi = int(xp); yi = int(yp + 0.5_num); zi = int(zp + 0.5_num)
    xf = xp - xi; yf = yp + 0.5_num - yi; zf = zp + 0.5_num - zi

    b1(0) = b1(0) + bx(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + bx(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    b1(0) = b1(0) + bx(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + bx(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    b1(0) = b1(0) + bx(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + bx(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    b1(0) = b1(0) + bx(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + bx(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !Interpolate by
    xi = int(xp+0.5_num); yi = int(yp); zi = int(zp + 0.5_num)
    xf = xp + 0.5_num - xi; yf = yp - yi; zf = zp + 0.5_num - zi

    b1(1) = b1(1) + by(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + by(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    b1(1) = b1(1) + by(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + by(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    b1(1) = b1(1) + by(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + by(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    b1(1) = b1(1) + by(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + by(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !Interpolate bz
    xi = int(xp+0.5_num); yi = int(yp+0.5_num); zi = int(zp)
    xf = xp + 0.5_num - xi; yf = yp + 0.5_num - yi; zf = zp - zi

    b1(2) = b1(2) + bz(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + bz(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    b1(2) = b1(2) + bz(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + bz(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    b1(2) = b1(2) + bz(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + bz(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    b1(2) = b1(2) + bz(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + bz(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !print*, bz(xi+1,yi+1,zi), bz(xi+1,yi+1,zi+1)
    !print*, xf, yf, zf,  b1
    !print*, '______________________'
    if (sqrt(sum(b1**2)) < weakness_limit) then
        null_point = 1
    end if

    b1 = b1/sqrt(sum(b1**2))

    END SUBROUTINE interpolate_bfield

    SUBROUTINE obtain_bz(x_point, y_point, z_point)
    REAL(num):: x_point, y_point, z_point !Coordinates of point to interpolate
    INTEGER:: xi, yi, zi !Cell indices in each dimension
    REAL(num):: xp, yp, zp !Coordinates in the respective dimentions
    REAL(num):: xf, yf, zf !Distance 'up' each cell

    bz1 = 0.0_num
    !Interpolate bz

    xp = nx*(x_point - x0)/(x1 - x0)
    yp = ny*(y_point - y0)/(y1 - y0)
    zp = nz*(z_point - z0)/(z1 - z0)

    xi = int(xp+0.5_num); yi = int(yp+0.5_num); zi = int(zp)
    xf = xp + 0.5_num - xi; yf = yp + 0.5_num - yi; zf = zp - zi

    bz1 = bz1 + bz(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + bz(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    bz1 = bz1 + bz(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + bz(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    bz1 = bz1 + bz(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + bz(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    bz1 = bz1 + bz(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + bz(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    END SUBROUTINE obtain_bz

    SUBROUTINE bfield_boundary

    !Populates the ghost points so the interpolator works. Actual values aren't really that important
    IMPLICIT NONE

    by(0:nx+1,0:ny,0) = by(0:nx+1,0:ny,1) - dz*(bz(0:nx+1,1:ny+1,0) - bz(0:nx+1, 0:ny,0))/dy
    bx(0:nx, 0:ny+1,0) = bx(0:nx,0:ny+1,1) - dz*(bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0))/dx

    !UPPER BOUNDARY (Zero Current)
    by(0:nx+1,0:ny,nz+1) = by(0:nx+1,0:ny,nz) + dz*(bz(0:nx+1,1:ny+1,nz) - bz(0:nx+1, 0:ny,nz))/dy
    bx(0:nx, 0:ny+1,nz+1) = bx(0:nx,0:ny+1,nz) + dz*(bz(1:nx+1,0:ny+1,nz) - bz(0:nx,0:ny+1,nz))/dx

    !x boundaries (Zero current, and zero flux)
    bz(0,0:ny+1,0:nz) = bz(1,0:ny+1,0:nz) - dx*(bx(0,0:ny+1,1:nz+1) - bx(0, 0:ny+1,0:nz))/dz
    by(0,0:ny,0:nz+1) = by(1, 0:ny,0:nz+1) - dx*(bx(0,1:ny+1,0:nz+1) - bx(0,0:ny,0:nz+1))/dy

    bz(nx+1,0:ny+1,0:nz) = bz(nx,0:ny+1,0:nz) + dx*(bx(nx,0:ny+1,1:nz+1) - bx(nx, 0:ny+1,0:nz))/dz
    by(nx+1,0:ny,0:nz+1) = by(nx, 0:ny,0:nz+1) + dx*(bx(nx,1:ny+1,0:nz+1) - bx(nx,0:ny,0:nz+1))/dy

    !y boundaries (Zero current, and zero flux)
    bz(0:nx+1,0,0:nz) = bz(0:nx+1, 1,0:nz) - dy*(by(0:nx+1,0,1:nz+1) - by(0:nx+1,0,0:nz))/dz
    bx(0:nx,0,0:nz+1) = bx(0:nx,1,0:nz+1) - dy*(by(1:nx+1,0,0:nz+1) - by(0:nx, 0,0:nz+1))/dx

    bz(0:nx+1,ny+1,0:nz) = bz(0:nx+1, ny,0:nz) + dy*(by(0:nx+1,ny,1:nz+1) - by(0:nx+1,ny,0:nz))/dz
    bx(0:nx,ny+1,0:nz+1) = bx(0:nx,ny,0:nz+1) + dy*(by(1:nx+1,ny,0:nz+1) - by(0:nx, ny,0:nz+1))/dx

    END SUBROUTINE bfield_boundary

    SUBROUTINE integrate_line(start, updown)

    IMPLICIT NONE

    REAL(num), DIMENSION(:):: start(0:2), pt(0:2)
    REAL(num):: avgj, maxs, bdist
    INTEGER:: lcount, updown, i
    lcount = 0
    pt = start; null_point = 0
    maxs = 0.0_num; bdist = 1e6
    current_line = 1e6
    line_j = 0.0_num; line_flh = 0.0_num; line_winding = 0.0_num; line_twist = 0.0_num
    do while (.true.)

        if ((pt(0) < x0 - 1e-6) .or. (pt(0) > x1 + 1e-6) .or. (pt(1) < y0 - 1e-6) .or. (pt(1) > y1 + 1e-6) .or. &
        (pt(2) < z0 - 1e-6) .or. (pt(2) > z1 + 1e-6) .or. (lcount > max_line_length-1)) then
            exit
        end if
        if (null_point > 0.5_num) then
            current_line(1:max_line_length-1,0:2) = 1e6
            exit
        end if
        current_line(lcount,0:2) = pt
        !Interpolate the various fields here
        call interpolate_bfield(pt(0), pt(1), pt(2))
        call interpolate_jfield(pt(0), pt(1), pt(2))

        pt = pt + updown*ds*b1  !B1 is obtained from interpolate_bfield

        line_j(lcount) = jmag

        CALL interpolate_generic_field(flh_density, pt(0), pt(1), pt(2), flh)
        line_flh(lcount) = flh

        CALL interpolate_generic_field(winding_density, pt(0), pt(1), pt(2), winding)
        line_winding(lcount) = winding

        CALL interpolate_generic_field(twist_density, pt(0), pt(1), pt(2), twist)
        line_twist(lcount) = twist

        lcount = lcount + 1
    end do

    if (lcount > 10) then
    !Disregard current outliers
    do i = 0, 9
        line_j(maxloc(line_j,1)-1) = 0.0_num
    end do
    end if
    line_j(0) = 0.0_num !The first one is often spuriously big
    avgj = sum(line_j)/lcount !Average current throughout the line

    current_line_length = lcount
    END SUBROUTINE

    SUBROUTINE update_emissivity(line, linej, lcount)
        !Using the line coordinates and the average current, update the emissivity matrix (established in GRID)
        REAL(num), DIMENSION(:,:):: line(0:max_line_length-1,0:2)
        REAL(num), DIMENSION(0:2):: pt
        REAL(num):: linej
        INTEGER:: lcount, n, xp, yp, zp
        INTEGER:: xprev, yprev, zprev
        INTEGER:: splurge

        xprev = -1; yprev = -1; zprev = -1
        !Run through the lines and update as appropriate. May be slow...
        splurge = 1
        do n = 0, lcount-1
             pt = line(n,:)
             if (maxval(pt) < 1d6) then
                xp = int(nx_out*(pt(0) - x0)/(x1 - x0))
                yp = int(ny_out*(pt(1) - y0)/(y1 - y0))
                zp = int(nz_out*(pt(2) - z0)/(z1 - z0))
                if (zp > splurge .and. zp < nz_out - splurge) then
                     emiss(xp-splurge:xp+splurge,yp-splurge:yp+splurge,zp-splurge:zp+splurge) = emiss(xp-splurge:xp+splurge,yp-splurge:yp+splurge,zp-splurge:zp+splurge) + linej
                end if
                emiss(xp,yp,zp) = emiss(xp,yp,zp) + linej

            end if
        end do

    END SUBROUTINE

    SUBROUTINE integrate_fieldlines()
    !Uses the 'starts' array to output an array of all the fieldlines
    !Traces in both directions
    IMPLICIT NONE

    integer:: start_index, line_number
    real(num), dimension(:,:):: line1(0:max_line_length-1,0:2), line2(0:max_line_length-1,0:2)
    real(num):: linej_total, winding_total, flh_total, twist_total, surface_at_point
    INTEGER:: length1, length2
    allocate(current_line(0:max_line_length-1,0:2))
    allocate(line_j(0:max_line_length-1))
    allocate(line_flh(0:max_line_length-1))
    allocate(line_winding(0:max_line_length-1))
    allocate(line_twist(0:max_line_length-1))
    !allocate(all_lines(0:2*nstarts-1, 0:max_line_length-1, 0:2))
    allocate(flh_array(0:nstarts-1), winding_array(0:nstarts-1), twist_array(0:nstarts-1), surface_array(0:nstarts-1), current_array(0:nstarts-1))
    line_number = 0

    do start_index = 0, nstarts-1
        linej_total = 0.0_num
        flh_total = 0.0_num; winding_total = 0.0_num; twist_total = 0.0_num

         !Integrated quantities should only really count once...

         !Integrate opposite to magnetic field.
         call integrate_line(starts(start_index,:),-1)
         line_number = line_number + 1
         line1 = current_line
         linej_total = linej_total + sum(line_j)
         flh_total = flh_total + sum(line_flh)
         winding_total = winding_total + sum(line_winding)
         twist_total = twist_total + sum(line_twist)

         length1 = current_line_length
         !Integrate with magnetic field
         call integrate_line(starts(start_index,:),1)
         line_number = line_number + 1
         line2 = current_line
         linej_total = linej_total + sum(line_j)
         flh_total = flh_total + sum(line_flh)
         winding_total = winding_total + sum(line_winding)
         twist_total = twist_total + sum(line_twist)

         length2 = current_line_length

         linej_total = linej_total/((length1 + length2))   !Average throughout for this one
         !Sum for the others
         flh_total = flh_total*ds
         winding_total = winding_total*ds
         twist_total = twist_total*ds

         call obtain_bz(starts(start_index,0), starts(start_index,1), starts(start_index,2))
         surface_at_point = bz1
         !Can change what you're plotting here, theoretically
         CALL update_emissivity(line1, linej_total, length1)
         CALL update_emissivity(line2, linej_total, length2)

        if (minval(current_line) < 1e6) then
        if (print_flag > 0.5_num .and. (mod(line_number, 100) == 0)) print*, 'Field line number',  line_number, 'of', nstarts*2, 'integrated', int(50*line_number/nstarts), '%'
        flh_array(start_index) = flh_total
        winding_array(start_index) = winding_total
        twist_array(start_index) = twist_total
        surface_array(start_index) = bz1
        current_array(start_index) = linej_total*(length1 + length2)*ds
    end if

     end do

    if (print_flag > 0.5_num) print*, nstarts*2, 'Field lines integrated'

    END SUBROUTINE integrate_fieldlines


END PROGRAM fltrace







