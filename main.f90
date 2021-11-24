program LF
    
    !Starting Modules
    use readingData
    use plyStress
    use failureAnalysis
    
    !Declaring Main Variables
    implicit none
    
    !Opening Laminate Data
    open(1,file='dataIn.txt',status='old')
    open(2,file='dataOut.txt',status='unknown')
    
    !Starting Reading and Getting Generalized Laminate Deformation 
    call reading()
    call getLaminateStiffnessMatrix()
    call getGeneralizedDeform()
    
    !Getting Stresses and Deformations of a Ply
    call allocatingPlyDeformAndStresses()
    call getPlyGenDeform()
    call getPlyDeform()
    call getPlyStress()
    
    !Analising First Failure
    call allocatingFailureCriteria()
    call failureCriteriaTsaiHill()
    call failureModeMaximumStress()
    call stepLoadToFailure()
    
    
    !Closing Program
    call closingProgram()
    pause
    
    end program LF
    
    
    subroutine closingProgram()
    
    write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(*,*) '           The program ended successfully.          '
    write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(*,*) '            Data stored in dataOut.txt              '
    write(*,*) '             in the project directory.              '
    write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    
    end subroutine closingProgram