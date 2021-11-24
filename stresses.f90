module plyStress
    
    use readingData
    
    implicit none
    real,allocatable,dimension(:,:) :: plyGenDeformSup,plyGenDeformInf
    real,allocatable,dimension(:,:) :: plyDeformSup,plyDeformInf
    real,allocatable,dimension(:,:) :: plyStressesSup,plyStressesInf
    
    contains
    
    subroutine allocatingPlyDeformAndStresses()
    
    allocate(plyGenDeformSup(3,numberOfPlies))
    allocate(plyGenDeformInf(3,numberOfPlies))
    allocate(plyDeformSup(3,numberOfPlies))
    allocate(plyDeformInf(3,numberOfPlies))
    allocate(plyStressesSup(3,numberOfPlies))
    allocate(plyStressesInf(3,numberOfPlies))
    
    end subroutine allocatingPlyDeformAndStresses
    
    subroutine getPlyGenDeform()
    
    !Declaring Variables
    integer :: i,j,k
    real,dimension(numberOfPlies+1,1) :: z
    
    !Calling Plies Quotas z Vector
    call getQuotasVector(z)
    
    !Sweeping Plies
    do k = 1,numberOfPlies
        
        !Superior Interface
        plyGenDeformSup(1,k) = genDeform(1,1) + z(k+1,1)*genDeform(4,1) 
        plyGenDeformSup(2,k) = genDeform(2,1) + z(k+1,1)*genDeform(5,1) 
        plyGenDeformSup(3,k) = genDeform(3,1) + z(k+1,1)*genDeform(6,1) 
        
        !Inferior Interface
        plyGenDeformInf(1,k) = genDeform(1,1) + z(k,1)*genDeform(4,1) 
        plyGenDeformInf(2,k) = genDeform(2,1) + z(k,1)*genDeform(5,1) 
        plyGenDeformInf(3,k) = genDeform(3,1) + z(k,1)*genDeform(6,1) 
        
    end do
    
    !Writting on Document
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '______________Ply GenDeformation (Superior Interface) ___________________'

    do j = 1,numberOfPlies
        write(2,*) ' '
        write(2,*) 'Ply :',j,'      Deform:    E1    E2    Y12'
        write(2,*) '------------------------------------------------------------------ '
        write(2,*) ' '
        do i = 1,3
            write(2,'(F22.4,X)',advance='no') plyGenDeformSup(i,j)
        end do
        write(2,*) ' '
    end do
    
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '________________Ply GenDeformation (Inferior Interface) __________________'

    do j = 1,numberOfPlies
        write(2,*) ' '
        write(2,*) 'Ply :',j,'     Deform:    E1    E2    Y12'
        write(2,*) '------------------------------------------------------------------ '
        write(2,*) ' '
        do i = 1,3
            write(2,'(F22.4,X)',advance='no') plyGenDeformInf(i,j)
        end do
        write(2,*) ' '
    end do
    
    end subroutine getPlyGenDeform
    
    subroutine getPlyDeform()
    
    !Declaring Variables
    integer :: k
    real :: theta
    real,dimension(3,3) :: rotT,invertedRotT
    
    !Sweeping Plies
    do k = 1,numberOfPlies
        
        theta = plyParameters(k,2)
        call getRotationMatrix(rotT,invertedRotT,theta)
        
        !Superior Interface
        plyDeformSup(:,k) = matmul(transpose(invertedRotT),plyGenDeformSup(:,k)) 
        
        !Inferior Interface
        plyDeformInf(:,k) = matmul(transpose(invertedRotT),plyGenDeformInf(:,k)) 
        
    end do
      
    end subroutine getPlyDeform
    
    subroutine getPlyStress()
    
    !Declaring Variables
    integer :: i,j,k
    real,dimension(3,3) :: Q
    
    !Sweeping Plies
    do k = 1,numberOfPlies
        
        call getQPlyMatrix(Q,k)    
        
        !Superior Interface
        plyStressesSup(:,k) = 1000000000*matmul(Q,plyDeformSup(:,k))
        
        !Inferior Interface
        plyStressesInf(:,k) = 1000000000*matmul(Q,plyDeformInf(:,k))
        
    end do
    
    !Writting on Document
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '________________Ply Stresses (Superior Interface) [MPa]___________________'

    do j = 1,numberOfPlies
        write(2,*) ' '
        write(2,*) 'Ply :',j,'     Stresses:    S1    S2    T12'
        write(2,*) '------------------------------------------------------------------ '
        write(2,*) ' '
        do i = 1,3
            write(2,'(F22.4,X)',advance='no') plyStressesSup(i,j)/1000000
        end do
        write(2,*) ' '
    end do
    
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '________________Ply Stresses (Inferior Interface) [MPa]___________________'

    do j = 1,numberOfPlies
        write(2,*) ' '
        write(2,*) 'Ply :',j,'     Stresses:    S1    S2    T12'
        write(2,*) '------------------------------------------------------------------ '
        write(2,*) ' '
        do i = 1,3
            write(2,'(F22.4,X)',advance='no') plyStressesInf(i,j)/1000000
        end do
        write(2,*) ' '
    end do
    
    end subroutine getPlyStress
    
    
    end module plyStress