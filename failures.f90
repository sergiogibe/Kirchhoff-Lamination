module failureAnalysis
    
    use readingData
    use plyStress
    
    implicit none
    integer :: failureCounter
    real :: smallerC
    real, dimension(6,1) :: maxLoading
    real,allocatable,dimension(:,:) :: safetyCArray
    real,allocatable,dimension(:,:) :: safetyCMatrix
    
    contains
    
    subroutine allocatingFailureCriteria()
    
    allocate(safetyCArray(1,numberOfPlies))
    allocate(safetyCMatrix(2,numberOfPlies))
    
    end subroutine allocatingFailureCriteria
    
    subroutine failureCriteriaTsaiHill()
    
    !Declaring Variables
    integer :: i,j,k
    real :: X,Y,S

    !Sweeping Plies
    do k = 1,numberOfPlies
        
        !Setting X,Y and S (Superior Interface)
        S = plyResistances(k,5)*1000000
        
        if (plyStressesSup(1,k) <= 0) then
            X = plyResistances(k,2)*1000000
        end if
        
        if (plyStressesSup(1,k) > 0) then
            X = plyResistances(k,1)*1000000
        end if
    
        if (plyStressesSup(2,k) <= 0) then
            Y = plyResistances(k,4)*1000000
        end if
    
        if (plyStressesSup(2,k) > 0) then
            Y = plyResistances(k,3)*1000000
        end if
        
        !Setting X,Y and S (Inferior Interface)
        if (plyStressesSup(1,k) <= 0) then
            X = plyResistances(k,2)*1000000
        end if
        
        if (plyStressesSup(1,k) > 0) then
            X = plyResistances(k,1)*1000000
        end if
    
        if (plyStressesSup(2,k) <= 0) then
            Y = plyResistances(k,4)*1000000
        end if
    
        if (plyStressesSup(2,k) > 0) then
            Y = plyResistances(k,3)*1000000
        end if
        
        !Getting safety C
        safetyCMatrix(1,k) = (1/(((plyStressesSup(1,k)/X)**2) - ((plyStressesSup(1,k)*plyStressesSup(2,k))/(X**2)) + ((plyStressesSup(2,k)/Y)**2) + ((plyStressesSup(3,k)/S)**2)))**(0.5)
        safetyCMatrix(2,k) = (1/(((plyStressesInf(1,k)/X)**2) - ((plyStressesInf(1,k)*plyStressesInf(2,k))/(X**2)) + ((plyStressesInf(2,k)/Y)**2) + ((plyStressesInf(3,k)/S)**2)))**(0.5)
        
    end do
    
    !Getting Safety C Array
    do i = 1,numberOfPlies
        
        if (safetyCMatrix(1,i) <= safetyCMatrix(2,i)) then
            safetyCArray(1,i) = safetyCMatrix(1,i)
        end if
        
        if (safetyCMatrix(1,i) > safetyCMatrix(2,i)) then
            safetyCArray(1,i) = safetyCMatrix(2,i)
        end if
        
    end do
    
    !Getting Plies Failures
    write(2,*) ' '
    write(2,*) '____________________________Failures________________________________ '
    write(2,*) ' '
    failureCounter = 0
    do i = 1,numberOfPlies  
        
        if (safetyCArray(1,i) <= 1) then
            write(2,*) 'Ply ', i, '        FAILURE'
            write(2,*) ' ------------------------------------------------------'
            write(2,*) ' '
            failureCounter = failureCounter + 1
        end if
        
        if (safetyCArray(1,i) > 1) then
            write(2,*) 'Ply ', i, '        does not fail.'
            write(2,*) ' ------------------------------------------------------'
            write(2,*) ' '
        end if
     
    end do
    
    !Getting Smaller C of the Array Considering NO Failure
    if (failureCounter < 1) then
        
        smallerC = safetyCArray(1,1)        
        do i = 1,numberOfPlies      
            if (smallerC >= safetyCArray(1,i)) then
                smallerC = safetyCArray(1,i)
            end if     
        end do
        
        !Getting the Max.Loading Until First Failure
        do i = 1,6
            maxLoading(i,1) = smallerC*loadingOfAPoint(i,1)
        end do
        
    end if
    
    end subroutine failureCriteriaTsaiHill
    
    subroutine failureModeMaximumStress()
    
    !Declaring Variables
    integer :: k
    
    write(2,*) ' '
    write(2,*) '____________________Failure Modes Maximum Stress___________________________'
    write(2,*) ' '
    
    !Sweeping Plies
    do k = 1,numberOfPlies
        
        if (plyResistances(k,2)*(-1000000) < plyStressesSup(1,k) .AND. plyStressesSup(1,k) < plyResistances(k,1)*1000000) then
            write(2,*) ''
        else
            write(2,*) '-------------------------------------------------------------------------'
            write(2,*) 'Ply ',k, 'FAILURE (MODE 1).'
        end if
        
        if (plyResistances(k,4)*(-1000000) < plyStressesSup(2,k) .AND. plyStressesSup(2,k) < plyResistances(k,3)*1000000) then
            write(2,*) ''
        else
            write(2,*) '--------------------------------------------------------------------------'
            write(2,*) 'Ply ',k, 'FAILURE (MODE 2).'
        end if
        
        if (plyResistances(k,5)*(-1000000) < plyStressesSup(3,k) .AND. plyStressesSup(3,k) < plyResistances(k,5)*1000000) then
            write(2,*) ''
        else
            write(2,*) '--------------------------------------------------------------------------'
            write(2,*) 'Ply ',k, 'FAILURE (MODE 3).'
        end if
           
    end do
    
    end subroutine failureModeMaximumStress
    
    subroutine stepLoadToFailure()
        
    !Declaring Variables
    integer :: i
    
    if (failureCounter < 1) then
    
    !Updating Loads
    loadingOfAPoint = maxLoading
    
    !Writing on Document
    write(2,*) ' '
    write(2,*) ' __________________________First Failure Loading____________________________'
    write(2,*) ' '
    do i = 1,6
        write(2,'(F20.4,X)',advance = 'no') loadingOfAPoint(i,1)
    end do
    
    !Writing on Document
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '%%%%%%%%%%%%%%%%%%%%       RE-RUNNING PROGRAM     %%%%%%%%%%%%%%%%%%%%%%%%%'
    write(2,*) ' '
    write(2,*) ' '
    write(2,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    write(2,*) ' '
    write(2,*) ' '
    
    !Re-running Program
    call getGeneralizedDeform()
    call getPlyGenDeform()
    call getPlyDeform()
    call getPlyStress()
    call failureCriteriaTsaiHill()
    call failureModeMaximumStress()
    
    end if
        
    end subroutine stepLoadToFailure
    
    end module failureAnalysis