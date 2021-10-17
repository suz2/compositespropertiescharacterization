SUBROUTINE umat (stress,  statev,  ddsdde,  sse,     spd,    &
                scd,     rpl,     ddsddt,  drplde,  drpldt,  &
                stran,  dstran, time,    dtime,   temp,      &
                dtemp,   predef,  dpred,   cmname,  ndi,     &
                nshr,    ntens,   nstatv,  props,   nprops,  &
                coords,  drot,    pnewdt,  celent,  dfgrd0,  &
                dfgrd1,  noel,    npt,     layer,   kspt,    &
                kstep,   kinc )

    implicit none                        
                     
!   Variables passed into the UMAT sub

    character*80, intent(in)  :: cmname
    integer(kind=8), intent(in)  :: ntens, nstatv, nprops,         &
                                        ndi, nshr, noel, npt,          &
                                        kspt, kstep, kinc, layer                                       
    real(kind=8) :: sse, spd, scd, rpl, drpldt,                    &
                        dtime, temp, dtemp, pnewdt, celent

!   Dimension arrays passed into the UMAT sub
    real(kind=8)        &
    stress(ntens),          &! Cauchy stress (vector form)
    statev(nstatv),         &! State variables
    ddsdde(ntens,ntens),    &! Tangent Stiffness Matrix
    ddsddt(ntens),          &! Change in stress per change in temperature
    drplde(ntens),          &! Change in heat generation per change in strain
    stran(ntens),           &! Strain tensor (vector form)
    dstran(ntens),          &! Strain increment tensor (vector form)
    time(2),                &! Time Step and Total Time
    predef(1),              &! Predefined state vars dependent on field variables
    dpred(1),               &! Change in predefined state variables
    props(nprops),          &! Material properties
    coords(3),              &! Coordinates of Gauss pt. being evaluated
    drot(3,3),              &! Incremental rotation matrix
    dfgrd0(3,3),            &! Deformation gradient at t_n
    dfgrd1(3,3)              ! Deformation gradient at t_(n+1)

    integer i, j, STAT, fibernum
    CHARACTER*200  FNAME,dir
    real(kind=8) :: Em, num, Ef, nuf, compliance_matrix(3,3), ddsdde_noshear(3,3), matrixprop(2)
    real(kind=8) :: Em_bar, Em_inter, alpha, ll, fiber_radius
    real(KIND=8), ALLOCATABLE :: dist(:), fiberlocation(:,:)
    real(kind=8), parameter :: pi = 3.141592653

!*******************
    Ef = 19500/1.d6
    nuf= 0.36

    Em_bar = 0.00506
    Em_inter = 0.0075426
    alpha = -0.23465
    num = 0.34

    fiber_radius = 5.d0
    dir = '/home/zsu/NASAIRAD/opt_stratafiedsampling_estimationobjfunc/'

!*******************

    
    stran = stran + dstran
    ddsdde = 0.d0

    IF (CMNAME.EQ.'2') THEN

        FNAME = TRIM(ADJUSTL(dir))//'fibercoord.dat'

        OPEN(555, FILE = FNAME, ACTION='READ', IOSTAT=STAT)
        READ(555,*) fibernum
        ALLOCATE (fiberlocation(2,fibernum),dist(fibernum))
        DO I = 1,fibernum
            READ(555,*) fiberlocation(1:2,i)
        END DO
        CLOSE(555, STATUS='KEEP') 

        DO I = 1,fibernum
            dist(i) = dsqrt((coords(1) - fiberlocation(1,i))**2.d0 + (coords(2) - fiberlocation(2,i))**2.d0)
        END DO
        ll = minval(dist(1:fibernum)) - fiber_radius

        Em = Em_bar + (Em_inter - Em_bar) * exp(alpha * ll)


!        FNAME = TRIM(ADJUSTL(dir))//'Em_gp.dat'
!        OPEN(555, FILE = FNAME, ACTION='WRITE', position="append", IOSTAT=STAT)
!            WRITE(555,*) Em,ll,coords(1),coords(2)
!        CLOSE(555, STATUS='KEEP')

        compliance_matrix(1,:) = (/1/Em,        -num/Em,   -num/Em/)
        compliance_matrix(2,:) = (/-num/Em,        1/Em,   -num/Em/)
        compliance_matrix(3,:) = (/-num/Em,     -num/Em,      1/Em/)
        ddsdde(4,4) = Em/(2*(1+num))
        ddsdde(5,5) = Em/(2*(1+num))
        ddsdde(6,6) = Em/(2*(1+num))

    ELSEIF (CMNAME.EQ.'1') THEN

        compliance_matrix(1,:) = (/1/Ef,        -nuf/Ef,   -nuf/Ef/)
        compliance_matrix(2,:) = (/-nuf/Ef,        1/Ef,   -nuf/Ef/)
        compliance_matrix(3,:) = (/-nuf/Ef,     -nuf/Ef,      1/Ef/)
        ddsdde(4,4) = Ef/(2*(1+nuf))
        ddsdde(5,5) = Ef/(2*(1+nuf))
        ddsdde(6,6) = Ef/(2*(1+nuf))

    ENDIF

 
    CALL INVERSION(compliance_matrix,ddsdde_noshear)

    ddsdde(1,1:3) = ddsdde_noshear(1,:)
    ddsdde(2,1:3) = ddsdde_noshear(2,:)
    ddsdde(3,1:3) = ddsdde_noshear(3,:)

    stress = MATMUL(ddsdde, stran)

    RETURN 
END

SUBROUTINE INVERSION(A,B)

!! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(kind=8), intent(in) :: A(3,3)   !! Matrix
    real(kind=8)             :: B(3,3)   !! Inverse matrix
    real(kind=8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

END SUBROUTINE


