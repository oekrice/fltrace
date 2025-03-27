!*******************************************************************************
MODULE grid
!*******************************************************************************
! Generates all the shared grid data, as grid3d does in the python code. Maybe additionally put shared grid in arrays in a different file - not sure yet. Most (if not all) of the grid arrays should not depend on the process, but I guess we'll see about that.
!*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

    contains

    subroutine establish_grid()
        !Establishes the grid arrays, which hopefully can be read in mainly from the netcdf. Perhaps.
        !Also read in the magnetic field here because may as well
        IMPLICIT NONE
        INTEGER:: ncid, vid
        INTEGER:: i,j,k
        character(len=5) :: snap_id
        allocate(xs(0:nx),ys(0:ny),zs(0:nz))
        allocate(xc(1:nx),yc(1:ny),zc(1:nz))
        allocate(xs_out(0:nx_out),ys_out(0:ny_out),zs_out(0:nz_out))

        allocate(bx(0:nx,-1:ny+1,-1:nz+1))
        allocate(by(-1:nx+1,0:ny,-1:nz+1))
        allocate(bz(-1:nx+1,-1:ny+1,0:nz))

        allocate(jx(0:nx+1,0:ny  ,0:nz  ))
        allocate(jy(0:nx  ,0:ny+1,0:nz  ))
        allocate(jz(0:nx  ,0:ny  ,0:nz+1))

        allocate(flh_density(0:nx+1,0:ny+1  ,0:nz+1))
        allocate(winding_density(0:nx+1,0:ny+1  ,0:nz+1))
        allocate(twist_density(0:nx+1,0:ny+1  ,0:nz+1))

        allocate(emiss(0:nx_out,0:ny_out,0:nz_out))

        flh_density = 0.0_num; winding_density = 0.0_num; twist_density = 0.0_num
        flh_flag  = .false.; winding_flag = .false.; twist_flag = .false.
        emiss = 0.0_num

        write (snap_id,'(I5.5)') snap
        bfield_filename = trim(trim(data_root)//trim(snap_id)//'.nc')

        if (print_flag > 0.5) print*, 'Importing magnetic field to Fortran, tmp fname: ', bfield_filename

        call try(nf90_open(trim(bfield_filename), nf90_nowrite, ncid))

        call try(nf90_inq_varid(ncid, 'bx', vid))
        call try(nf90_get_var(ncid, vid, bx(0:nx,1:ny,1:nz)))

        call try(nf90_inq_varid(ncid, 'by', vid))
        call try(nf90_get_var(ncid, vid, by(1:nx,0:ny,1:nz)))

        call try(nf90_inq_varid(ncid, 'bz', vid))
        call try(nf90_get_var(ncid, vid, bz(1:nx,1:ny,0:nz)))


        call try(nf90_inq_varid(ncid, 'jx', vid))
        call try(nf90_get_var(ncid, vid, jx(1:nx,0:ny,0:nz)))

        call try(nf90_inq_varid(ncid, 'jy', vid))
        call try(nf90_get_var(ncid, vid, jy(0:nx,1:ny,0:nz)))

        call try(nf90_inq_varid(ncid, 'jz', vid))
        call try(nf90_get_var(ncid, vid, jz(0:nx,0:ny,1:nz)))

        call try(nf90_inq_varid(ncid, 'xs', vid))
        call try(nf90_get_var(ncid, vid, xs(0:nx)))

        call try(nf90_inq_varid(ncid, 'ys', vid))
        call try(nf90_get_var(ncid, vid, ys(0:ny)))

        call try(nf90_inq_varid(ncid, 'zs', vid))
        call try(nf90_get_var(ncid, vid, zs(0:nz)))

        !Read in field line things if possible
        call try2(nf90_inq_varid(ncid, 'flh_density', vid), flh_flag)
        if (flh_flag) call try(nf90_get_var(ncid, vid, flh_density(1:nx,1:ny,1:nz)))

        call try2(nf90_inq_varid(ncid, 'winding_density', vid), winding_flag)
        if (winding_flag) call try(nf90_get_var(ncid, vid, winding_density(1:nx,1:ny,1:nz)))

        call try2(nf90_inq_varid(ncid, 'twist_density', vid), twist_flag)
        if (twist_flag) call try(nf90_get_var(ncid, vid, twist_density(1:nx,1:ny,1:nz)))

        call try(nf90_close(ncid))

        print*, 'Field line quantities imported?', flh_flag, winding_flag, twist_flag
        x0 = xs(0); x1 = xs(nx); y0 = ys(0); y1 = ys(ny); z0 = zs(0); z1 = zs(nz)

        xc(1:nx) = 0.5_num*(xs(0:nx-1) + xs(1:nx))
        yc(1:ny) = 0.5_num*(ys(0:ny-1) + ys(1:ny))
        zc(1:nz) = 0.5_num*(zs(0:nz-1) + zs(1:nz))

        dx = sum((xs(1:nx) - xs(0:nx-1)))/ nx
        dy = sum((ys(1:ny) - ys(0:ny-1)))/ ny
        dz = sum((zs(1:nz) - zs(0:nz-1)))/ nz

        do i = 0, nx_out
            xs_out(i) = xs(0) + (xs(nx) - xs(0))*i/nx_out
        end do

        do j = 0, ny_out
            ys_out(j) = ys(0) + (ys(ny) - ys(0))*j/ny_out
        end do

        do k = 0, nz_out
            zs_out(k) = zs(0) + (zs(nz) - zs(0))*k/nz_out
        end do

        !call check_divergence

        if (print_flag > 0.5) print*, 'Grid established and magnetic field read-in'

        ds = ds_factor*min(dx, dy, dz)  !Tracing 'timestep'

    end subroutine establish_grid

    subroutine check_divergence
    !Checks that the solenoidal condition is fine in every cell
    IMPLICIT none
    real(num):: div(1:nx,1:ny,1:nz)

    div = 0.0_num
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (bx(1:nx,1:ny,1:nz) - bx(0:nx-1,1:ny,1:nz))/dx
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (by(1:nx,1:ny,1:nz) - by(1:nx,0:ny-1,1:nz))/dy
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (bz(1:nx,1:ny,1:nz) - bz(1:nx,1:ny,0:nz-1))/dz

    if (maxval(abs(div)) > 1d-10) then
        print*, 'Read-in magnetic field is not divergence free, stopping'
        STOP
    end if

    end subroutine check_divergence

    SUBROUTINE export_emissivity
    IMPLICIT NONE

    character(len=64):: filename
    integer:: aid, bid, cid, did, vid, ncid
    integer:: xid, yid, zid
    integer:: fid, gid, hid, iid, jid
    CHARACTER(LEN=5):: snap_id

    real(num), dimension(:,:):: xsum(0:ny_out,0:nz_out), ysum(0:nx_out,0:nz_out), zsum(0:nx_out,0:ny_out)

    write (snap_id,'(I5.5)') snap

    filename = trim('./fl_data/emiss'//trim(snap_id)//'.nc')

    call try(nf90_create(trim(filename), nf90_clobber, ncid))

    !Define variables
    call try(nf90_def_dim(ncid, 'nx_out', nx_out+1, aid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'ny_out', ny_out+1, bid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'nz_out', nz_out+1, cid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'nstarts', nstarts, did))  !Make up fake dimensions here

    if (save_all > 0.5_num) call try(nf90_def_var(ncid, 'emiss', nf90_double, (/aid, bid, cid/), vid))
    call try(nf90_def_var(ncid, 'emiss_xsum', nf90_double, (/bid, cid/), xid))
    call try(nf90_def_var(ncid, 'emiss_ysum', nf90_double, (/aid, cid/), yid))
    call try(nf90_def_var(ncid, 'emiss_zsum', nf90_double, (/aid, bid/), zid))

    call try(nf90_def_var(ncid, 'flh_array', nf90_double, (/did/), fid))
    call try(nf90_def_var(ncid, 'winding_array', nf90_double, (/did/), gid))
    call try(nf90_def_var(ncid, 'twist_array', nf90_double, (/did/), hid))
    call try(nf90_def_var(ncid, 'surface_array', nf90_double, (/did/), iid))
    call try(nf90_def_var(ncid, 'current_array', nf90_double, (/did/), jid))

    call try(nf90_enddef(ncid))

    xsum = sum(emiss, dim = 1)
    ysum = sum(emiss, dim = 2)
    zsum = sum(emiss, dim = 3)

    !Write variables
    if (save_all > 0.5_num) call try(nf90_put_var(ncid, vid, emiss))
    call try(nf90_put_var(ncid, xid, xsum))
    call try(nf90_put_var(ncid, yid, ysum))
    call try(nf90_put_var(ncid, zid, zsum))

    call try(nf90_put_var(ncid, fid, flh_array))
    call try(nf90_put_var(ncid, gid, winding_array))
    call try(nf90_put_var(ncid, hid, twist_array))
    call try(nf90_put_var(ncid, iid, surface_array))
    call try(nf90_put_var(ncid, jid, current_array))

    call try(nf90_close(ncid))

    if (print_flag > 0.5_num) print*, 'Emissivity exported to file', filename

    END SUBROUTINE export_emissivity

    subroutine export_fieldlines
    !Output the calculated magnetic field lines as a netcdf file flines.nc
    IMPLICIT NONE

    character(len=64):: filename
    integer:: aid, bid, cid, vid, ncid, i
    CHARACTER(LEN=5):: run_id

    write (run_id,'(I5.5)') snap

    filename = trim('./fl_data/flines'//trim(run_id)//'.nc')

    call try(nf90_create(trim(filename), nf90_clobber, ncid))

    !Define variables
    call try(nf90_def_dim(ncid, 'a', 3, aid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'b', max_line_length, bid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'c', nstarts*2, cid))  !Make up fake dimensions here

    call try(nf90_def_var(ncid, 'lines', nf90_double, (/cid,bid,aid/), vid))
    call try(nf90_enddef(ncid))

    !Write variables
    call try(nf90_put_var(ncid, vid, all_lines))
    call try(nf90_close(ncid))

    if (print_flag > 0.5_num) print*, 'Field lines exported to file ', filename

    return

    end subroutine export_fieldlines


    SUBROUTINE try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        print*, 'Ensure data directory is correct in fltrace.f90'
        stop
    end if

    END SUBROUTINE try

    SUBROUTINE try2(status, flag)
    ! Check if there is flh stuff to be calculated
    INTEGER, INTENT(IN):: status
    LOGICAL:: flag
    if (status /= NF90_noerr) THEN
        flag = .false.
    else
        flag = .true.
    end if
    END SUBROUTINE try2

END MODULE grid
