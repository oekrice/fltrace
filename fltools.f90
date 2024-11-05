function integrate_line(bx,by,bz,start,max_line_length,nx,ny,nz,x0,x1,y0,y1,z0,z1,ds,weakness_limit,updown) result(line)

    REAL(8), intent(in)    :: bx(:,:,:), by(:,:,:), bz(:,:,:)
    REAL(8), intent(in)    :: start(0:2)
    REAL(8):: line(0:max_line_length-1,0:2), pt(0:2)
    REAL(8):: b1(0:2)
    INTEGER:: nx, ny, nz, updown, max_line_length, lcount, null_point
    REAL(8):: x0, x1, y0, y1, z0, z1, mag, ds, weakness_limit

    REAL(8):: xp, yp, zp !Coordinates in the respective dimentions
    INTEGER:: xi, yi, zi !Cell indices in each dimension
    REAL(8):: xf, yf, zf !Distance 'up' each cell

    line = 1e6_8
    pt = start; lcount = 0; null_point = 0
    do while (.true.)
        if ((pt(0) < x0) .or. (pt(0) > x1) .or. (pt(1) < y0) .or. (pt(1) > x1) .or. &
        (pt(2) < z0) .or. (pt(2) > z1) .or. (lcount > max_line_length-1) .or. (null_point > 0.5_8)) then
            exit
        end if
        line(lcount, 0:2) = pt

        b1 = 0.0_8
        !Establish ratios
        xp = nx*(pt(0) - x0)/(x1 - x0)
        yp = ny*(pt(1) - y0)/(y1 - y0)
        zp = nz*(pt(2) - z0)/(z1 - z0)

        !Interpolate bx
        xi = int(xp+1.0_8); yi = int(yp + 1.5_8); zi = int(zp + 1.5_8)
        xf = xp - xi + 1.0_8; yf = yp + 1.5_8 - yi; zf = zp + 1.5_8 - zi

        b1(0) = b1(0) + bx(xi,yi,zi)*(1.0_8-xf)*(1.0_8-yf)*(1.0_8-zf) + bx(xi,yi,zi+1)*(1.0_8-xf)*(1.0_8-yf)*(zf)
        b1(0) = b1(0) + bx(xi,yi+1,zi)*(1.0_8-xf)*(yf)*(1.0_8-zf)       + bx(xi,yi+1,zi+1)*(1.0_8-xf)*(yf)*(zf)
        b1(0) = b1(0) + bx(xi+1,yi,zi)*(xf)*(1.0_8-yf)*(1.0_8-zf)       + bx(xi+1,yi,zi+1)*(xf)*(1.0_8-yf)*(zf)
        b1(0) = b1(0) + bx(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_8-zf)             + bx(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

        !Interpolate by
        xi = int(xp+1.5_8); yi = int(yp+1.0_8); zi = int(zp + 1.5_8)
        xf = xp + 1.5_8 - xi; yf = yp +1.0_8 - yi; zf = zp + 1.5_8 - zi

        b1(1) = b1(1) + by(xi,yi,zi)*(1.0_8-xf)*(1.0_8-yf)*(1.0_8-zf) + by(xi,yi,zi+1)*(1.0_8-xf)*(1.0_8-yf)*(zf)
        b1(1) = b1(1) + by(xi,yi+1,zi)*(1.0_8-xf)*(yf)*(1.0_8-zf)       + by(xi,yi+1,zi+1)*(1.0_8-xf)*(yf)*(zf)
        b1(1) = b1(1) + by(xi+1,yi,zi)*(xf)*(1.0_8-yf)*(1.0_8-zf)       + by(xi+1,yi,zi+1)*(xf)*(1.0_8-yf)*(zf)
        b1(1) = b1(1) + by(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_8-zf)             + by(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

        !Interpolate bz
        xi = int(xp+1.5_8); yi = int(yp+1.5_8); zi = int(zp +1.0_8)
        xf = xp + 1.5_8 - xi; yf = yp + 1.5_8 - yi; zf = zp - zi +1.0_8

        b1(2) = b1(2) + bz(xi,yi,zi)*(1.0_8-xf)*(1.0_8-yf)*(1.0_8-zf) + bz(xi,yi,zi+1)*(1.0_8-xf)*(1.0_8-yf)*(zf)
        b1(2) = b1(2) + bz(xi,yi+1,zi)*(1.0_8-xf)*(yf)*(1.0_8-zf)       + bz(xi,yi+1,zi+1)*(1.0_8-xf)*(yf)*(zf)
        b1(2) = b1(2) + bz(xi+1,yi,zi)*(xf)*(1.0_8-yf)*(1.0_8-zf)       + bz(xi+1,yi,zi+1)*(xf)*(1.0_8-yf)*(zf)
        b1(2) = b1(2) + bz(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_8-zf)             + bz(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

        mag = sqrt(sum(b1**2))
        b1 = b1/mag
        if (mag < weakness_limit) then
            null_point = 1
        end if
        pt = pt + updown*ds*b1
        lcount = lcount + 1
    end do

 END function integrate_line
!

! SUBROUTINE bfield_boundary
!
! !Populates the ghost points so the interpolator works. Actual values aren't really that important
! IMPLICIT NONE
!
! by(0:nx+1,0:ny,0) = by(0:nx+1,0:ny,1) - dz*(bz(0:nx+1,1:ny+1,0) - bz(0:nx+1, 0:ny,0))/dy
! bx(0:nx, 0:ny+1,0) = bx(0:nx,0:ny+1,1) - dz*(bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0))/dx
!
! !UPPER BOUNDARY (Zero Current)
! by(0:nx+1,0:ny,nz+1) = by(0:nx+1,0:ny,nz) + dz*(bz(0:nx+1,1:ny+1,nz) - bz(0:nx+1, 0:ny,nz))/dy
! bx(0:nx, 0:ny+1,nz+1) = bx(0:nx,0:ny+1,nz) + dz*(bz(1:nx+1,0:ny+1,nz) - bz(0:nx,0:ny+1,nz))/dx
!
! !x boundaries (Zero current, and zero flux)
! bz(0,0:ny+1,0:nz) = bz(1,0:ny+1,0:nz) - dx*(bx(0,0:ny+1,1:nz+1) - bx(0, 0:ny+1,0:nz))/dz
! by(0,0:ny,0:nz+1) = by(1, 0:ny,0:nz+1) - dx*(bx(0,1:ny+1,0:nz+1) - bx(0,0:ny,0:nz+1))/dy
!
! bz(nx+1,0:ny+1,0:nz) = bz(nx,0:ny+1,0:nz) + dx*(bx(nx,0:ny+1,1:nz+1) - bx(nx, 0:ny+1,0:nz))/dz
! by(nx+1,0:ny,0:nz+1) = by(nx, 0:ny,0:nz+1) + dx*(bx(nx,1:ny+1,0:nz+1) - bx(nx,0:ny,0:nz+1))/dy
!
! !y boundaries (Zero current, and zero flux)
! bz(0:nx+1,0,0:nz) = bz(0:nx+1, 1,0:nz) - dy*(by(0:nx+1,0,1:nz+1) - by(0:nx+1,0,0:nz))/dz
! bx(0:nx,0,0:nz+1) = bx(0:nx,1,0:nz+1) - dy*(by(1:nx+1,0,0:nz+1) - by(0:nx, 0,0:nz+1))/dx
!
! bz(0:nx+1,ny+1,0:nz) = bz(0:nx+1, ny,0:nz) + dy*(by(0:nx+1,ny,1:nz+1) - by(0:nx+1,ny,0:nz))/dz
! bx(0:nx,ny+1,0:nz+1) = bx(0:nx,ny,0:nz+1) + dy*(by(1:nx+1,ny,0:nz+1) - by(0:nx, ny,0:nz+1))/dx
!
! END SUBROUTINE bfield_boundary
!
! SUBROUTINE integrate_line(start, updown, line_8ber)
!
! IMPLICIT NONE
!
! REAL(8), DIMENSION(:):: start(0:2), pt(0:2)
! REAL(8), DIMENSION(:,:):: line(0:max_line_length-1,0:2)
! INTEGER:: lcount, updown, line_8ber
! line = 1e6; lcount = 0
!
! pt = start; null_point = 0
! do while (.true.)
!     if ((pt(0) < x0) .or. (pt(0) > x1) .or. (pt(1) < y0) .or. (pt(1) > x1) .or. &
!     (pt(2) < z0) .or. (pt(2) > z1) .or. (lcount > max_line_length-1) .or. (null_point > 0.5_8)) then
!         exit
!     end if
!     line(lcount,0:2) = pt
!     call interpolate_bfield(pt(0), pt(1), pt(2))
!     pt = pt + updown*ds*b1
!     lcount = lcount + 1
!
! end do
!
! END SUBROUTINE
!
!
!
!










