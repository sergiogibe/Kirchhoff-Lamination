module readingData
    
    implicit none
    integer, public :: numberOfPlies
    real, allocatable, dimension(:,:) :: plyParameters
    real, allocatable, dimension(:,:) :: plyResistances
    real, dimension(6,1) :: loadingOfAPoint
    real, dimension(2,1) :: loadingOfAPointQ
    real,dimension(6,6) :: laminateStiffness
    real,dimension(6,6) :: laminateFlexibility
    real,dimension(6,1) :: genDeform
    
    
contains
    subroutine reading()
    integer :: i,j
    
    !Reading Plies Data 
    read(1,*) numberOfPlies
    allocate(plyParameters(numberOfPlies,8))
    allocate(plyResistances(numberOfPlies,5))
    do i=1,numberOfPlies
            read(1,*) plyParameters(i,1),plyParameters(i,2),plyParameters(i,3),plyParameters(i,4),plyParameters(i,5),plyParameters(i,6),plyParameters(i,7),plyParameters(i,8)
    end do
    do i=1,numberOfPlies
            read(1,*) plyResistances(i,1),plyResistances(i,2),plyResistances(i,3),plyResistances(i,4),plyResistances(i,5)
    end do
    
    !Reading Loading of a Point
    read(1,*) loadingOfAPoint(1,1),loadingOfAPoint(2,1),loadingOfAPoint(3,1),loadingOfAPoint(4,1),loadingOfAPoint(5,1),loadingOfAPoint(6,1)
    read(1,*) loadingOfAPointQ(2,1),loadingOfAPointQ(2,1)
    
    !Writing Plies Data
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '%%%%%%%%%%%%%%%%%%%%%%       RUNNING PROGRAM     %%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '______________________PARAMETERS_____________________________'
    write(2,*) ' '
    write(2,*) 'Thickness[m] | Orientation[degrees] | E1[GPa] | E2[GPa] | G12[GPa] | v12 | G23[GPa] | v23'
    write(2,*) ' '
    do i=1,numberOfPlies
        write(2,*) 'Ply ',i
        write(2,*) ' '
        write(2,*) '---------------------------------------------------------- ' 
        write(2,*) ' '
        do j=1,8
            write(2,'(F10.3,X)',advance='no') plyParameters(i,j)
        end do
        write(2,*) ' '
    end do
    write(2,*) ' '
    write(2,*) '________________________RESISTANCES____________________________'
    write(2,*) ' '
    write(2,*) 'Xt[MPa] | Xc[MPa] | Yt[MPa] | Yc[MPa] | S[MPa]'
    write(2,*) ' '
    do i=1,numberOfPlies
        write(2,*) ' '
        write(2,*) 'Ply ' , i
        write(2,*) ' '
        write(2,*) '----------------------------------------------------------- ' 
        write(2,*) ' '
        do j=1,5
            write(2,'(F10.3,X)',advance='no') plyResistances(i,j)
        end do
        write(2,*) ' '
    end do
    
    !Writing Loads
    write(2,*) ' '
    write(2,*) '_____________________________LOADING_____________________________'
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) 'Nx[N/m] | Ny[N/m] | Nxy[N/m] | Mx[N] | My[N] | Mxy[N]'
    write(2,*) ' '
    do i=1,6
        write(2,'(F10.3,X)',advance='no') loadingOfAPoint(i,1)
    end do 
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '----------------------------------------------------------------- '
    write(2,*) ' '
    write(2,*) 'Qx[N/m]  |  Qy[N/m]'
    write(2,*) ' '
    do i=1,2
        write(2,'(F10.3,X)',advance='no') loadingOfAPointQ(i,1)    
    end do
    write(2,*) ' '
    
    
    end subroutine reading
    
    subroutine getQPlyMatrix(QPlyMatrix,plyNumber)
    
    !Declaring Variables
    integer :: i,j,plyNumber
    real :: E1,E2,G12,V12
    real,dimension(3,3) :: QPlyMatrix
    
    !Zeros Q Ply Matrix
    do i = 1,3
        do j = 1,3
            QPlyMatrix(i,j) = 0
        end do
    end do
    
    !Setting Engineering Constants 
    E1 = plyParameters(plyNumber,3)
    E2 = plyParameters(plyNumber,4)
    G12 = plyParameters(plyNumber,5)
    v12 = plyParameters(plyNumber,6)
    
    !Setting Q Ply Matrix
    QPlyMatrix(1,1) = (E1**2)/(E1-(v12**2)*E2)
    QPlyMatrix(1,2) = (v12*E1*E2)/(E1-(v12**2)*E2)
    QPlyMatrix(2,1) = (v12*E1*E2)/(E1-(v12**2)*E2)
    QPlyMatrix(2,2) = (E1*E2)/(E1-(v12**2)*E2)
    QPlyMatrix(3,3) = G12
    
    end subroutine getQPlyMatrix
    
    subroutine getRotationMatrix(rotT,invertedRotT,theta)
    
    !Declaring Variables
    real :: theta
    real,dimension(3,3) :: rotT,invertedRotT
    
    !Setting the rotation T
    rotT(1,1) = (cosd(theta))**2
    rotT(1,2) = (sind(theta))**2
    rotT(1,3) = 2*sind(theta)*cosd(theta)
    rotT(2,1) = (sind(theta))**2
    rotT(2,2) = (cosd(theta))**2
    rotT(2,3) = -2*sind(theta)*cosd(theta)
    rotT(3,1) = -1*sind(theta)*cosd(theta)
    rotT(3,2) = sind(theta)*cosd(theta)
    rotT(3,3) = (cosd(theta))**2 - (sind(theta))**2
    
    !Setting the iverted rotation T
    invertedRotT(1,1) = (cosd(theta))**2
    invertedRotT(1,2) = (sind(theta))**2
    invertedRotT(1,3) = -2*sind(theta)*cosd(theta)
    invertedRotT(2,1) = (sind(theta))**2
    invertedRotT(2,2) = (cosd(theta))**2
    invertedRotT(2,3) = 2*sind(theta)*cosd(theta)
    invertedRotT(3,1) = sind(theta)*cosd(theta)
    invertedRotT(3,2) = -1*sind(theta)*cosd(theta)
    invertedRotT(3,3) = (cosd(theta))**2 - (sind(theta))**2
    
    end subroutine getRotationMatrix
    
    subroutine getQuotasVector(zVector)
    
    !Declaring Variables
    integer :: i
    real :: laminateThickness
    real,dimension(numberOfPlies+1,1) :: zVector
    
    !Getting Total Thickness
    laminateThickness = 0.0
    do i = 1,numberOfPlies
        laminateThickness = laminateThickness + plyParameters(i,1)
    end do
    
    !Getting Quotas Vector 
    zVector(1,1) = -1*(laminateThickness/2.0)
    
    do i = 1,numberOfPlies
        zVector(i+1,1) = zVector(i,1) + plyParameters(i,1)    
    end do
    
    end subroutine getQuotasVector
    
    subroutine getLaminateStiffnessMatrix()
    
    !Declaring Variables
    integer :: i,j,k
    real :: theta
    real,dimension(3,3) :: Q,barQ,rotT,invertedRotT
    real,dimension(3,3) :: A,B,D
    real,dimension(numberOfPlies+1,1) :: z
    
    !Calling Quotas Z Vector
    
    call getQuotasVector(z)
    
    !Setting zeros A,B and D Matrix
    do i = 1,3
        do j = 1,3
            A(i,j) = 0
            B(i,j) = 0
            D(i,j) = 0
        end do
    end do
    
    !Sweeping Laminate to get A, B and D
    do k = 1,numberOfPlies  
        
        !Getting barQ
        call getQPlyMatrix(Q,k)
        theta = plyParameters(k,2)
        call getRotationMatrix(rotT,invertedRotT,theta) 
        barQ = matmul(matmul(invertedRotT,Q),transpose(invertedRotT))
        
        !Writting barQ
        write(2,*) ' '
        write(2,*) ' ________________________Rotated Q______________________________'
        write(2,*) ' '
        write(2,*) 'Ply ',k
        write(2,*) ' '
        write(2,*) ' --------------------------------------------------------------'
        write(2,*) ' '
        do i=1,3
        do j=1,3
            write(2,'(F15.3,X)',advance='no') barQ(i,j)
        end do
        write(2,*) ' '
        end do
        
        !Sweeping barQ
        do i = 1,3
            do j = 1,3
                A(i,j) = A(i,j) + barQ(i,j) * (z(k+1,1) - z(k,1))
                B(i,j) = B(i,j) + 0.5 * barQ(i,j) * ((z(k+1,1))**2 - (z(k,1))**2)
                D(i,j) = D(i,j) + (1.0/3.0) * barQ(i,j) * ((z(k+1,1))**3 - (z(k,1))**3)
            end do
        end do 
        
    end do
    
    !Assembly of Laminate Stiffness Matrix
    do i = 1,3
        do j = 1,3
            laminateStiffness(i,j) = A(i,j)
            laminateStiffness(i,j+3) = B(i,j)
            laminateStiffness(i+3,j) = B(i,j)
            laminateStiffness(i+3,j+3) = D(i,j)
        end do
    end do
    
    !Writting on Document
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '_____________________________Stiffness Matrix________________________________'
    write(2,*) ' '
    do i=1,6
        write(2,*) ' '
        do j=1,6
            write(2,'(F20.5,X)',advance='no') laminateStiffness(i,j)
        end do
        write(2,*) ' '
    end do
    
    end subroutine getLaminateStiffnessMatrix
    
    subroutine inv(a,c,n)
!============================================================
! Inverse matrix
! Method: Based LU decomposition for Ax=b
! Marouf Abderahmane - Strasbourg University 2016
!-----------------------------------------------------------
implicit none
integer :: n,i, j, k
real,dimension(n,n) :: a,c,U,L
real,dimension(n) :: b,d,x
real :: coeff

L=0.0
U=0.0
b=0.0
do k=1, n-1
do i=k+1,n
coeff=a(i,k)/a(k,k)
L(i,k) = coeff
do j=k+1,n
a(i,j) = a(i,j)-coeff*a(k,j)
end do
end do
end do
do i=1,n
L(i,i) = 1.0
end do
do j=1,n
do i=1,j
U(i,j) = a(i,j)
end do
end do
do k=1,n
b(k)=1.0
d(1) = b(1)
do i=2,n
d(i)=b(i)
do j=1,i-1
d(i) = d(i) - L(i,j)*d(j)
end do
end do
x(n)=d(n)/U(n,n)
do i = n-1,1,-1
x(i) = d(i)
do j=n,i+1,-1
x(i)=x(i)-U(i,j)*x(j)
end do
x(i) = x(i)/u(i,i)
end do
do i=1,n
c(i,k) = x(i)
end do
b(k)=0.0
end do
end subroutine inv
    
    subroutine getGeneralizedDeform()
    
    integer :: i,j
    real,dimension(6,1) :: genDPrint
    
    !Calling invertion
    call inv(laminateStiffness,laminateFlexibility,6)
    
    !Getting Generalized Deformation of Laminate
    genDPrint = matmul(laminateFlexibility,loadingOfAPoint)
    
    do i = 1,6
        genDeform(i,1) = genDPrint(i,1)/1000
    end do 
    
    !Writting on Document
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '__________________________Flexibility Matrix________________________________'
    write(2,*) ' '
    do i=1,6
        do j=1,6
            write(2,'(F20.5,X)',advance='no') laminateFlexibility(i,j)
        end do
        write(2,*) ' '
    end do
    
    !Writting Generalized Deformation
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) ' _________________________Generalized Deformation (e-3)______________________'
    write(2,*) ' '
    do i=1,6
            write(2,'(F20.4,X)',advance='no') genDPrint(i,1)
    end do
    write(2,*) ' '

    
    end subroutine getGeneralizedDeform
    
    end module readingData