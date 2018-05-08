!================================================================
!	Efield v1.0
!
!   - Arbitrary field values
!   - X-direction only
!
! Tom Halverson
! 1-29-2018
!================================================================


Program Main

  !----Modules-----
  Use Inputs
  Use Globals
  !----------------

  Implicit None
  Double Precision :: t1,t2,tot_time
  Character(len = 5) :: units

  Print *,'========================================================'
  Print *,'|     Grizzly Bear v1.0                                |'
  Print *,'|     Created by: Thomas Halverson                     |'
  Print *,'|     Date: 1-29-18                                    |'
  Print *,'|                                                      |'
  Print *,'|  Computes the eigenstates of a system of             |'
  Print *,'|  dipole-dipole coupled rotors in a time varying      |'
  Print *,'|  electric field using the eigenstates of the         |'
  Print *,'|  free rotor: | LM >                                  |'
  Print *,'========================================================'
  Print *

  Call CPU_time(t1)

  Call GetInputs() 
  Call Map_Indices()
  Call MakeR()
  Call Propagate()

  Call CPU_time(t2)

  If(t2-t1 > 18000.0d0) then
    tot_time = (t2-t1)/(60.0d0*60.0d0)
    units = 'hrs.'
  Else if (t2-t1 > 300) then
    tot_time = (t2-t1)/(60.0d0)
    units = 'mins.'
  Else 
    tot_time = (t2-t1)
    units = 'sec.'
  End If

  Print *,'*****************************************************************'
  Print *,'      SUMMARY'
  Print *,'  ------------------------------------'
  Print *
  Print *, "   ",Dmax,' Coupled Rotors'
  Print 1001, LMax
  Print 1003, Big_Len
  Print 1004, Trim(Vals_FN)
  Print 1005, Trim(Vecs_FN)
  Print *
  Print *
   
  Print 9999, tot_time, units
  Print *,'*****************************************************************'
  Print *
  Print *
  
  9999 Format ('  Total time taken: ', f10.5, 2x, a)
  1001 Format ('  Max L value:                 ', 3x, i5)
  1002 Format ('  Max Nu value:                ', 3x, i5)
  1003 Format ('  Total number of eigenvalues: ', 3x, i7)
  1004 Format ('  Eigenvalues read to file:    ', a)
  1005 Format ('  Eigenvectors read to file:   ', a)



End Program Main


!==========================================================================================
  Subroutine GetInputs()  
  
   !----Modules-----
    Use Inputs
    Use Globals
   !----------------

    Implicit None
    Character(len = 128) :: input_file

    Integer :: i, j, n
    
    Call GETARG(1, input_file)
    input_file = trim(input_file)

    !input_file = "testinput.dat"

    Open(UNIT = 15, FILE = input_file, STATUS='OLD',  ACTION ='Read')

      Print *,'*****************************************************************'
      Print *,'STEP 1:'
      Print '(1x,a,2x,a20)',"Reading data from: ", input_file

      Read(15,'(A)')Big_FN
      Big_FN = trim(Big_FN)
      Read(15,'(A)')Geos_FN
      Geos_FN = trim(Geos_FN)
      Read(15,'(A)')Params_FN
      Params_FN = trim(Params_FN)
      Read(15,'(A)')E_FN
      E_FN = trim(E_FN)
      Read(15,'(A)') Vals_FN
      Vals_FN = trim(Vals_FN)
      Read(15,'(A)') Vecs_FN
      Vecs_FN = trim(Vecs_FN)
 
    Close(15)   

    !-----------------------------------------
    !  Read in ND basis fucntions from file
    Open(UNIT = 15, FILE = big_FN, STATUS='OLD',  ACTION ='Read')
    Read (15,"(i3,2x,i3,2x,i6)") Dmax, LMax, big_Len
    !Print *, Dmax, LMax, big_Len
    Allocate(big_LM(big_Len,2*DMax))
    Do i = 1, Big_Len
      Do j = 1, 2*Dmax
        Read (15,'(i3,2x)',Advance='NO')   big_LM(i, j)
      End Do 
      Read (15,*)   
    End Do
    Close(15)
    !-----------------------------------------
    

    !-----------------------------------------
    !  Create V12 for generic pair
    Call MakeV()
    !-----------------------------------------


    !-----------------------------------------
    !  Read in gemotries 
    !(in carteisan coordinates)
    Open(UNIT = 15, FILE = Geos_FN, STATUS='OLD',  ACTION ='Read')
    Allocate(Geos(3*Dmax))
    
    Do i = 1, 3*DMax
      Read(15,*) Geos(i)
    End Do
        
    Close(15) 
    !-----------------------------------------


     
    !-----------------------------------------
    !  Read mu and B from file
    Open(UNIT = 15, FILE = params_FN, STATUS='OLD',  ACTION ='Read')
    Read(15,*) B_Const
    Read(15,*) mu_Const	
    !-----------------------------------------

    
    !-----------------------------------------
    !  Read mu and B from file
    Open(UNIT = 15, FILE = E_FN, STATUS='OLD',  ACTION ='Read')
    
    Read(15,*) E_Len
    Allocate(E(2,E_len))
 
    Do i = 1, E_len
      Read (15,*) E(1,i), E(2,i)
    End Do
    !-----------------------------------------
    

    Print *, "Data read: SUCCESSFULL"
    Print *
    Print 1001, "Number of Rotors:                 ", Dmax
    Print 1001, "Maximum L value:                  ", LMax
    Print 1001, 'Number ND basis functions:        ', big_Len
    Print 1001, 'Number of time step grid points:  ', E_Len
    Print 1004, 'Rotational constant (in Hartree): ', B_const
    Print 1004, 'Dipole moment (in atomic units):  ', mu_const
    Print *
    Print *, '  -----------------------------------------------'
    Print *, '        Center of Mass Geometries (in Bohr)        '
    n = 1
    Do i = 1, Dmax
      Write(*,'(5x,a1,i2,a1,3x)',Advance='NO') 'P',i,':'
      Do j = 1, 3
        Write(*,'(f8.5,2x)',Advance='NO') Geos(n)
        n = n + 1
      End Do
      Print *
    End Do
    Print *, '  -----------------------------------------------'
    Print *
    Print 1002, "Eigenvals filename: ", Vals_FN
    Print 1002, "Eigenvals filename: ", Vecs_FN
    Print *,'*****************************************************************'
    Print *
    Print *
    Print *

    1000 Format(2x, a, i13)
    1001 Format(2x, a, i13)
    1002 Format(2x, a, a30)
    1003 Format(10x, a, a1)
    1004 Format(2x, a, f12.10)
    1005 Format(2x, a, f12.6)
    
  End Subroutine GetInputs
!==========================================================================================


!==========================================================================================
  Subroutine MakeV()

   !----Modules-----
    Use Inputs
    Use Globals
   !----------------
  
  	Implicit None

    Integer :: count, LEN
    Integer :: I,J
    Integer :: l1,m1,l2,m2
    Integer :: l1p,m1p,l2p,m2p

    Double Precision :: t1, t2
    Double Precision :: x_term, z_term, dummy
    Double Complex :: y_term, dummy_I

    Double Precision :: XELEM
    Double Complex   :: YELEM
    Double Precision :: ZELEM

    count = 0          
    Do l1 = 0, Lmax
      Do m1 = -l1, l1
        Do l2 = 0, Lmax !- l1
          Do m2 = -l2, l2
            count = count + 1
          End Do
        End Do
      End Do
    End Do

    Small_Len = count

    Allocate(Small_LM(Small_Len,4))
    Allocate(V12(Small_Len,Small_Len))

    Small_LM = 0
    V12 = 0.0d0

    count = 1
    Do l1 = 0, Lmax
      Do m1 = -l1, l1
        Do l2 = 0, Lmax !- l1
          Do m2 = -l2, l2
          
            Small_LM(count,1) = l1
            Small_LM(count,2) = m1
            Small_LM(count,3) = l2
            Small_LM(count,4) = m2

            count = count + 1
         
          End Do
        End Do
      End Do
    End Do

    Do I = 1, Small_Len

      l1p = Small_LM(I, 1)
      m1p = Small_LM(I, 2)
      l2p = Small_LM(I, 3)
      m2p = Small_LM(I, 4)
    
      Do J = 1, Small_Len

        l1 = Small_LM(J, 1)
        m1 = Small_LM(J, 2)
        l2 = Small_LM(J, 3)
        m2 = Small_LM(J, 4)

        x_term = XELEM(l1p,m1p,l1,m1)*XELEM(l2p,m2p,l2,m2)
        y_term = YELEM(l1p,m1p,l1,m1)*YELEM(l2p,m2p,l2,m2)
        z_term = ZELEM(l1p,m1p,l1,m1)*ZELEM(l2p,m2p,l2,m2)

        V12(I,J) = x_term + Real(y_term) - 2.0d0*z_term
         
      End Do
    End Do
    
  End Subroutine
!==========================================================================================



!========================================================================================== 
  Subroutine Map_Indices()

   !----Modules-----
    Use Inputs
    Use Globals
   !----------------

    Implicit NONE
    Integer :: i

    Allocate(IndexMap(0:lmax,-lmax:lmax,0:lmax,-lmax:lmax))

    IndexMap = 0

    Do i = 1, Small_Len
      IndexMap(Small_LM(i,1),Small_LM(i,2),Small_LM(i,3),Small_LM(i,4)) = i
    End Do
  
  End Subroutine Map_Indices
!========================================================================================== 



!==========================================================================================
  Subroutine MakeR()
  
   !----Modules-----
    Use Inputs
    Use Globals
   !----------------

    Implicit NONE
    Integer :: i,j,k
    Integer :: index,count
    Integer :: l,m,mp

    Double Precision :: PI = 3.14159265358979d0

    Double Precision :: ex1,ex2
    Double Precision :: why1,why2
    Double Precision :: zee1,zee2
    Double Precision :: tol, dummyT, dummyP
    Double Precision :: wig_im, wig_re

    Double Precision :: WignerD
    Double Precision :: Dot
    
    Integer, Dimension(20) :: Choose

    Double Precision, Dimension(3) :: rij
    Double Precision, Dimension(3) :: xhat
    Double Precision, Dimension(3) :: yhat
    Double Precision, Dimension(3) :: zhat
    
    tol = 0.000000000000001d0
        
    Choose = (/0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190/)

    xhat = (/1,0,0/)
    yhat = (/0,1,0/)
    zhat = (/0,0,1/)
    
    If (Dmax > 20) Then
      Print *
      Print *, '    ============================'
      Print *, '                ERROR:'
      Print *, '         Add More Factorials'
      Print *, '    ============================'
      Print *
    End if 

    count = Choose(Dmax)
    Allocate(Rijs(count,3))
    Allocate(Wigners(count,0:Lmax,-Lmax:Lmax,-Lmax:Lmax))
    Allocate(Wigners_Star(count,0:Lmax,-Lmax:Lmax,-Lmax:Lmax))
    
    !---------------------------------------------------------------------------------
    index = 1
    Do j = 1, Dmax-1
      Do k = j+1,Dmax
        
        ex1 = Geos(3*j-2)
        ex2 = Geos(3*k-2)
        why1 = Geos(3*j-1)
        why2 = Geos(3*k-1)
        zee1 = Geos(3*j)
        zee2 = Geos(3*k)

        rij(1) = ex2 - ex1
        rij(2) = why2 - why1
        rij(3) = zee2 - zee1
            
        Rijs(index,1) = DSQRT((ex1-ex2)**2+(why1-why2)**2+(zee1-zee2)**2)

        Do i = 1, 3
          rij(i) = rij(i)/Rijs(index,1)
        End Do

        Rijs(index,2) = ACos(Dot(rij,zhat,3))
        Rijs(index,3) = Atan2(Dot(yhat,rij - (Dot(rij,zhat,3)*zhat),3),Dot(xhat,rij - (Dot(rij,zhat,3)*zhat),3))
        
        index = index + 1
          
        End Do
      End Do
      
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      Wigners = Dcmplx(0.0d0,0.0d0)
      Wigners_Star = Dcmplx(0.0d0,0.0d0)
      
      Do index = 1, count
        Do l = 0, Lmax
          Do mp = -l, l
            Do m = -l, l
  
              wig_re = Dcos(-m*Rijs(index,3))*WignerD(l,mp,m,Rijs(index,2))
              wig_im = Dsin(-m*Rijs(index,3))*WignerD(l,mp,m,Rijs(index,2))

              Wigners(index,l,mp,m) = Dcmplx(wig_re,wig_im)
              Wigners_Star(index,l,mp,m) = Dcmplx(wig_re,-wig_im)
          
            End Do
          End Do
        End Do
      End Do

     !---------------------------------------------------------------------------------

    
  End Subroutine MakeR
!==========================================================================================



!==========================================================================================
  Subroutine Propagate()
   !----Modules-----
    Use Inputs
    Use Globals
    Use Outputs
   !----------------  
 

    Implicit NONE
    Integer :: n
    
    Double Precision :: t1, t2, tot_time 
    
    Character(len = 5) :: units
    
    Print *,'*****************************************************************'
    Print *,'STEP: 2'
    Print *,'  Computing eigenspectrum'

    Call CPU_time(t1)

    Do n = 1, E_len
      Allocate(H(Big_Len, Big_Len))
      Allocate(Eigenvals(Big_Len))
      
      Call MakeHam(E(2,n))
      Call Diagonalize()
      Call WriteVals(n)
      
      Deallocate(H)
      Deallocate(Eigenvals)
    End Do

    Call CPU_time(t2)

    If(t2-t1 > 18000.0d0) then
      tot_time = (t2-t1)/(60.0d0*60.0d0)
      units = 'hrs.'
    Else if (t2-t1 > 300) then
      tot_time = (t2-t1)/(60.0d0)
      units = 'mins.'
    Else 
      tot_time = (t2-t1)
      units = 'sec.'
    End If

    Print *, '  Spectrum creation: COMPLETE'
    Print '(a, f10.5, 2x, a)', '  Total time taken: ', tot_time, units
    Print *,'*****************************************************************'  
    Print *
    Print *
    Print *
 
  End Subroutine Propagate
!==========================================================================================



    
!==========================================================================================
  Subroutine MakeHam(gamma)
   !----Modules-----
    Use Inputs
    Use Globals
    Use Outputs
   !----------------  

    Implicit NONE
    Double precision, Intent(IN) :: gamma 
  
    Integer :: i, j, k, count, index
    Integer :: index1,index2
    Integer :: eye, jay, kay
    Integer :: lp, mp, l, m 
    Integer :: m1p, m1, m2p, m2
    Integer :: U1p, U1, U2p, U2
    Integer :: V1p, V1, V2p, V2
 
    Double Complex :: Zdummy,dummy3
    Double Complex :: d1p, d1, d2p ,d2

    Double Precision :: dummy1, dummy2
    Double Precision :: dist, tol = 0.0000000000001d0
    Double Precision :: t1, t2, tot_time
    Character(len = 5) :: units

    
    Integer, Dimension(Dmax) :: ELPs
    Integer, Dimension(Dmax) :: EMPs
    Integer, Dimension(Dmax) :: ELs
    Integer, Dimension(Dmax) :: EMs

    Double Precision, Allocatable, Dimension(:,:,:,:) :: Overlap
    Double Precision, Allocatable, Dimension(:,:,:,:) :: E_MAT
   
    Double Precision :: XELEM
    
    !Print *,'*****************************************************************'
    !Print *,'STEP: 2'
    !Print *,'  Creating Hamiltonian'
    !Print 1001, ' Array size:          ', Big_Len, ' x ', Big_Len
    !Print *
    !Print '(2x,a,f6.3)', 'gamma = ', gamma
    !Print *, '  -----------------------------------------------'
    !Print *, '            RIJ Components (R, th, phi)         '
    !count = 1
    !Do i = 1, Dmax-1
    !  Do j = i+1, Dmax
    !    Write(*,'(5x,a1,i2,i2,a1,3x)',Advance='NO') 'R',i,j,':'
    !    Do k = 1, 3
    !      Write(*,'(f10.7,2x)',Advance='NO') Rijs(count,k)
    !    End Do
    !    Print *
    !    count = count + 1
    !  End Do
      
    !End Do
    !Print *, '  -----------------------------------------------'
    !Print *

    Call CPU_time(t1)
    
    Allocate(Overlap(0:Lmax,-Lmax:Lmax,0:Lmax,-Lmax:Lmax))
    Allocate(E_MAT(0:Lmax,-Lmax:Lmax,0:Lmax,-Lmax:Lmax))

    H = (0.0d0,0.0d0)
    Overlap = 0.0d0
    E_MAT = 0.0d0

  

    Do lp = 0, Lmax
      Do mp = -Lmax, Lmax
        Do l = 0, Lmax
          Do m = -Lmax, Lmax

            If(lp==l.AND.mp==m)Then
              Overlap(lp,mp,l,m) = 1.0d0
            End If
            E_MAT(lp,mp,l,m) = XELEM(lp,mp,l,m)
          End Do
        End Do
      End Do
    End Do

   
  
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Calculates Kenetic Energy  >>>>>>>>>>>>>>>>>>>>>>>>>
    Do i = 1, Big_len
      dummy1 = 0.0d0
      Do eye = 1, Dmax
        dummy1 = dummy1 + Big_LM(i,2*eye-1)*(Big_LM(i,2*eye-1) + 1)
      End Do
      H(i,i) = DCmplx(B_const*dummy1,0.0d0)
    End Do
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  

  

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Calculates Potential Energy  >>>>>>>>>>>>>>>>>>>>>>>>>
    Do i = 1, Big_Len

      !NOTE: These arrays are unnecsary, and can be removed. 
      !      They just make book keeping easier  
      Do eye = 1, Dmax
         ELPs(eye) = Big_LM(i,2*eye-1)
         EMPs(eye) = Big_LM(i,2*eye)  
      End Do  

      
      Do j = i, Big_Len
        
        !NOTE: These arrays are unnecsary, and can be removed. 
        !      They just make book keeping easier
        Do eye = 1, Dmax
          ELs(eye) = Big_LM(j,2*eye-1)
          EMs(eye) = Big_LM(j,2*eye)
        End Do 
        
        dummy3 = DCMPLX(0.0d0,0.0d0)
        count = 1
        
        Do eye = 1, Dmax-1

              
          U1p = ELPs(eye)
          U1  = ELs(eye)
          V1p = EMPs(eye)
          V1  = EMs(eye)
         
          Do jay = eye + 1, Dmax

            U2p = ELPs(jay)
            U2  = ELs(jay)
            V2p = EMPs(jay)
            V2  = EMs(jay)

            dist = 1/(Rijs(count,1))**3

                    !If(i==2.AND.j==7) then
                    !  Print*,-2.0d0*YELEM(U1p,V1p,U1,V1)*YELEM(U2p,V2p,U2,V2)
                    !  Print *
                    !End If
          
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !            Computing the VIJ term for RIJ config
            Zdummy = (0.0d0,0.0d0)
            Do m1p = -U1p, U1p
              d1p = Wigners_Star(count,U1p,m1p,V1p)
              Do m2p = -U2p, U2p
                d2p = Wigners_Star(count,U2p,m2p,V2p)
                Do m1 = -U1, U1
                  d1 = Wigners(count,U1,m1,V1)
                  Do m2 = -U2, U2
                    d2 = Wigners(count,U2,m2,V2)
                          
                    index1 = IndexMap(U1p,m1p,U2p,m2p) 
                    index2 = IndexMap(U1,m1,U2,m2)                

                    Zdummy = Zdummy + d1p*d2p*d1*d2*DCMPLX(V12(index1,index2),0.0d0)

                    !If(i==2.AND.j==7) then
                    !  If(V12(index1,index2).NE.0.0d0) then
                    !    Print'(2x,f5.2,a3,f5.2,a1,&
                    !    &      3x,f5.2,a3,f5.2,a1,&
                    !    &      3x,f5.2,a3,f5.2,a1,&
                    !    &      3x,f5.2,a3,f5.2,a1,&
                    !    &      3x,f5.2,a3,f5.2,a1,&
                    !    &      5x,f5.2,&
                    !    &      5x,f5.2,a3,f5.2,a1)', &
                    !    &  Real(d1p),' + ',aimag(d1p),'i', &
                    !    &  Real(d2p),' + ',aimag(d2p),'i', &
                    !    &  Real(d1),' + ',aimag(d1),'i', &
                    !    &  Real(d2),' + ',aimag(d2),'i', &
                    !    &  Real(d1*d2p),' + ',aimag(d1*d2p),'i', &
                    !    &  V12(index1,index2), &
                    !    &  Real(Zdummy),' + ',aimag(Zdummy),'i'
                    !  End If 
                    !End If
                                                                   
                    End Do 
                  End Do
                End Do 
              End Do
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        	  !Print*, Aimag(Zdummy)
              Do kay = 1, eye - 1
                Zdummy = Zdummy * Dcmplx(Overlap(ELPs(kay),EMPs(kay),ELs(kay),EMs(kay)),0.0d0)
              End Do
              Do kay = eye + 1, jay - 1
                Zdummy = Zdummy * Dcmplx(Overlap(ELPs(kay),EMPs(kay),ELs(kay),EMs(kay)),0.0d0)
              End Do
              Do kay = jay + 1, Dmax
                Zdummy = Zdummy * Dcmplx(Overlap(ELPs(kay),EMPs(kay),ELs(kay),EMs(kay)),0.0d0)
              End Do
              
              !If(ABS(Aimag(Zdummy))>0.00000000001) then
              !  Print*, '******PROBLEM*******  ', i,j
              !End If 
      
              dummy3 = dummy3 + dist*Zdummy
              count = count + 1
              
            End Do
          End Do


        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !        Compute the field
        dummy2 = 0.0d0 
        Do eye = 1, Dmax

          dummy1 = E_Mat(ELPs(eye),EMPs(eye),ELs(eye),EMs(eye))

          Do jay = 1, eye-1
            dummy1 = dummy1 * Overlap(ELPs(jay),EMPs(jay),ELs(jay),EMs(jay))
          End Do

          Do jay = eye+ 1, Dmax
            dummy1 = dummy1 * Overlap(ELPs(jay),EMPs(jay),ELs(jay),EMs(jay))
          End Do
          
          dummy2 = dummy2 + dummy1

        End Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

        H(i,j) =  H(i,j) + DCMPLX(mu_const,0.0d0)*DCMPLX(mu_const,0.0d0)*dummy3 &
                  & + DCMPLX(gamma*dummy2,0.0d0)
        
      End Do
    End Do

    Do i = 1, Big_Len
      Do j = 1, Big_len
        H(j,i) = DCMPLX(Real(H(i,j)),-Aimag(H(i,j)))
      End Do
    End Do

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    Call CPU_time(t2)

    If(t2-t1 > 18000.0d0) then
      tot_time = (t2-t1)/(60.0d0*60.0d0)
      units = 'hrs.'
    Else if (t2-t1 > 300) then
      tot_time = (t2-t1)/(60.0d0)
      units = 'mins.'
    Else 
      tot_time = (t2-t1)
      units = 'sec.'
    End If

    !Print *, '  Hamiltonian creation: COMPLETE'
    !Print 9999, t2-t1, units
    !Print *,'*****************************************************************'  
    !Print *
    !Print *
    !Print *

    Deallocate(Overlap)
    Deallocate(E_MAT)
   

    !Print *
    !Print *
    !Do i = 1, Big_Len
    !  Do j = 1, Big_Len
    !    Write(*,'(f5.2,a3,f5.2,a1,2x)',Advance='NO') Real(H(i,j)),' + ',Aimag(H(i,j)),'i'
    !  End Do
    !  Print *
    !End Do
    !Print *
    !Print *
    !Print *
    !Do i = 1, Length
    !  Print *,i, H(i,2) 
    !End Do
    !Print *

    !l = 0
    !Do i = 1, Big_Len
    !  Do j = i, Big_Len
    !    !Abs(Real(H(i,j)))>tol.OR.
    !    If(Abs(Aimag(H(i,j)))>tol) then
    !      Print '(2x,i3,2x,i3,5x,f13.8,5x,f12.8)', i,j, Real(H(i,j)),Aimag(H(i,j))
    !       !l = l+1
    !    End If
    !  End Do
    !End Do
    !Print *
    !Print *,l
    !Print *

    
 
    1001 Format(2x, a, i10, 5x, a, i10)
    9999 Format ('   Total time taken: ', f10.5, 2x, a)
    
  End Subroutine MakeHam
!==========================================================================================


!==========================================================================================
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine Diagonalize()
!   - Diagonalize subroutine calculates the eigenvalues of the Hamiltonian
!    using the Scalapack routine ZHEEV
!
!  ZHEEV(JOBZ, UPLO, N, A, LDA, W, Work, LWork, RWork, Info)
!
!------------------------------------------------------------------------------
  SubRoutine Diagonalize()

     !----Modules-----
      Use Inputs
      Use Globals
      Use Outputs
     !----------------

	  Implicit NONE
   
      Integer           :: LWORK, INFO, LDA, i, j
      Character*1       :: JOBZ, UPLO
      Double Precision  :: t1, t2, tot_time
      
      Character(len = 5) :: units

      Double Precision, Allocatable, Dimension(:) :: RWork
      Double Complex,   Allocatable, Dimension(:) :: Work

      JOBZ = 'V'
      UPLO = 'U'
      LDA = 1
   
      Call CPU_time(t1)

      
      !Print *,'*****************************************************************'
      !Print *,'STEP: 3'
      !Print *,'  Diagonalizing Hamiltonian'
      !Print *,'  Lapack routine: ZHEEV'
      !Print *  
      
      LWORK = -1
      Allocate(WORK(1))
      Allocate(RWork((3*Big_Len)-2))
      
      Call ZHEEV(JOBZ, UPLO, Big_Len, H, Big_Len, EigenVals, Work, LWork, RWork, Info)
                        
      LWORK = WORK(1)
      
      Deallocate(WORK)
      Allocate(WORK(LWORK))
      
      !Print *, '       Workspace Sizes:'
      !Print *, '------------------------------------'
      !Print *, '  Work array length: ',LWORK
      !Print *, '------------------------------------'
      !Print *
    
      Call ZHEEV(JOBZ, UPLO, Big_Len, H, Big_Len, EigenVals, Work, LWork, RWork, Info)

      Deallocate(WORK)
      Deallocate(RWork)

      Call CPU_time(t2)

      If(t2-t1 > 18000.0d0) then
        tot_time = (t2-t1)/(60.0d0*60.0d0)
        units = 'hrs.'
      Else if (t2-t1 > 300) then
        tot_time = (t2-t1)/(60.0d0)
        units = 'mins.'
      Else 
        tot_time = (t2-t1)
        units = 'sec.'
      End If

   
      !Print *,'  Diagonalization: COMPLETE'
      !Print 9999, t2-t1, units
      !Print *,'*****************************************************************'
      !Print *
      !Print *
      !Print *
   
     
        
    
      9999 Format ('  Diagonalization time: ',f10.5,2x,a)

End SubRoutine
!==========================================================================================





!==========================================================================================
  Subroutine WriteVals(time)
    
   !----Modules-----
    Use Inputs
    Use Globals
    Use Outputs
   !----------------   

    Implicit NONE
    Integer, INTENT(IN) :: time
    Integer	:: i,j 

    Character(len = 128) :: FNAME
    Character(len = 128) :: FNAME_vals
    Character(len = 128) :: FNAME_vecs  
    Character(len = 10) :: blank
 
    If(time < 10) then
      Write(blank,'(i1)') time
      FNAME = '-t'//Trim(blank)//'.dat'
      FNAME=Trim(FNAME)
    Else If(time >= 10 .AND. time < 100) then
      Write(blank,'(i2)') time
      FNAME = '-t'//Trim(blank)//'.dat'
      FNAME=Trim(FNAME)
    Else If(time >= 100 .AND. time < 1000) then
      Write(blank,'(i3)') time
      FNAME = '-t'//Trim(blank)//'.dat'
      FNAME=Trim(FNAME)
    Else If(time >= 1000 .AND. time < 10000) then
      Write(blank,'(i4)') time
      FNAME = '-t'//Trim(blank)//'.dat'
      FNAME=Trim(FNAME)
    Else If(time >= 10000 .AND. time < 100000) then
      Write(blank,'(i5)') time
      FNAME = '-t'//Trim(blank)//'.dat'
      FNAME=Trim(FNAME)
    End If

    FNAME_vals = Trim(vals_FN)//FNAME
    FNAME_vecs = Trim(vecs_FN)//FNAME  
      
    Open(UNIT = 15, FILE = FNAME_vals, ACTION='WRITE')    
    
    Do i = 1, Big_Len
      !Write (15,'(f22.16)') 219474.63*EigenVals(i)
      !Write (15,'(f22.16,a2)') EigenVals(i),'d0'
      Write (15,'(f22.16)') EigenVals(i)
    End Do


    Close(UNIT = 15)


    !Print *
    !Print *, "==========================================="
    !Do i = 1, Big_Len
    !  Print '(2x, f6.4)', EigenVals(i) 
    !  Do j = 1, Big_len
    !    Write (*,'(i2,2x,i2,3x,f22.16)') i,j,H(i,j)
    !  End Do
    !Print *
    !End Do
    !Print *, "==========================================="
    !Print *



    Open(UNIT = 15, FILE = FNAME_vecs, ACTION='WRITE')    
    
    Do i = 1, Big_Len
      Do j = 1, Big_Len
        !Write (15,'(f22.16,a2,2x,f22.16,a2)') Real(H(j,i)),'d0',Aimag(H(j,i)),'d0'
        Write (15,'(f22.16,2x,f22.16)') Real(H(j,i)),Aimag(H(j,i))
      End Do
    End Do

    !Do i = 1, Big_Len
    !  Write (15,'(f22.16,2x,f22.16)') Real(H(i,2)),Aimag(H(i,2))
    !End Do

    !Do i = 1, Big_Len
    !  Write (15,'(f22.16)') H(i,2)
    !End Do

    !Do i = 1, Big_Len
    !  Write (15,'(f22.16)') H(i,3)
    !End Do

    !Do i = 1, Big_Len
    !  Write (15,'(f22.16)') H(i,4)
    !End Do



    Close(UNIT = 15)
   

End Subroutine WriteVals
!==========================================================================================









!==========================================================================================
Double Precision Function XELEM(lp,mp,l,m)

  Implicit NONE

  Integer, INTENT(IN) :: lp
  Integer, INTENT(IN) :: mp
  Integer, INTENT(IN) :: l
  Integer, INTENT(IN) :: m

  Double Precision :: A, B, C, D
  Double precision :: top, bot
  Double Precision :: out

  Double Precision :: Del

  top = DBLE((l+m+1)*(l+m+2))
  bot = DBLE(((2*l)+1))*DBLE((2*l)+3)
  A = 0.50d0*SQRT(top/bot)*Del(lp,l+1)*Del(mp,m+1)

  top = DBLE((l-m-1)*(l-m))
  bot = DBLE(((2*l)-1))*DBLE((2*l)+1)
  B = 0.5d0*SQRT(top/bot)*Del(lp,l-1)*Del(mp,m+1)

  top = DBLE((l-m+1)*(l-m+2))
  bot = DBLE(((2*l)+1))*DBLE((2*l)+3)
  C = 0.5d0*SQRT(top/bot)*Del(lp,l+1)*Del(mp,m-1)
  
  top = DBLE((l+m-1)*(l+m))
  bot = DBLE(((2*l)-1))*DBLE((2*l)+1)
  D = 0.5d0*SQRT(top/bot)*Del(lp,l-1)*Del(mp,m-1)

  out = -A + B + C - D

  XELEM = out

End Function XELEM
!==========================================================================================
Double Complex Function YELEM(lp,mp,l,m)

  Implicit NONE
  
  Integer, INTENT(IN) :: lp
  Integer, INTENT(IN) :: mp
  Integer, INTENT(IN) :: l
  Integer, INTENT(IN) :: m

  Double Precision :: A, B, C, D
  Double precision :: top, bot
  Double Precision :: out

  Double Precision :: Del

  top = DBLE((l+m+1)*(l+m+2))
  bot = DBLE(((2*l)+1))*DBLE((2*l)+3)
  A = 0.50d0*SQRT(top/bot)*Del(lp,l+1)*Del(mp,m+1)

  top = DBLE((l-m-1)*(l-m))
  bot = DBLE(((2*l)-1))*DBLE((2*l)+1)
  B = 0.5d0*SQRT(top/bot)*Del(lp,l-1)*Del(mp,m+1)

  top = DBLE((l-m+1)*(l-m+2))
  bot = DBLE(((2*l)+1))*DBLE((2*l)+3)
  C = 0.5d0*SQRT(top/bot)*Del(lp,l+1)*Del(mp,m-1)
  
  top = DBLE((l+m-1)*(l+m))
  bot = DBLE(((2*l)-1))*DBLE((2*l)+1)
  D = 0.5d0*SQRT(top/bot)*Del(lp,l-1)*Del(mp,m-1)

  out = A - B + C - D


  YELEM = DCMPLX(0.0d0,out)


End Function YELEM
!==========================================================================================
Double Precision Function ZELEM(lp,mp,l,m)

  Implicit NONE

  Integer, INTENT(IN) :: lp
  Integer, INTENT(IN) :: mp
  Integer, INTENT(IN) :: l
  Integer, INTENT(IN) :: m

  Double Precision :: A, B
  Double precision :: top, bot
  Double Precision :: out

  Double Precision :: Del

  top = DBLE((l-m+1)*(l+m+1))
  bot = DBLE(((2*l)+1))*DBLE((2*l)+3)
  A = SQRT(top/bot)*Del(lp,l+1)*Del(mp,m)

  top = DBLE((l-m)*(l+m))
  bot = DBLE(((2*l)-1))*DBLE((2*l)+1)
  B = SQRT(top/bot)*Del(lp,l-1)*Del(mp,m)

  ZELEM = A + B

End Function ZELEM
!==========================================================================================
Double Precision Function Del(a,b)

  Integer, INTENT(IN) :: a
  Integer, INTENT(IN) :: b

  Double Precision :: out

  If(a==b) Then
    out = 1.0d0
  Else
    out = 0.0d0
  End If

  Del = out

End Function Del
!==========================================================================================
!==========================================================================================
Double Precision Function WignerD(l,mp,m,theta)

  Implicit NONE
  Integer, Intent(IN) :: l
  Integer, Intent(IN) :: mp
  Integer, Intent(IN) :: m

  Double Precision, Intent(IN) :: theta

  Integer :: n
  Integer :: start, finish
  Integer :: exp1, exp2
  Double Precision :: A, B, C
  Double Precision :: cosine, sine
  Double Precision :: tot

  Double Precision :: Factorial

  cosine = Dcos(theta/2.0d0)
  sine = Dsin(theta/2.0d0)

  A = DSQRT(Factorial(l+m)*Factorial(l-m)*Factorial(l+mp)*Factorial(l-mp))

  start = Max(0,m - mp)
  finish = Min(l-mp,l+m)
  
  tot = 0.0d0

  Do n = start, finish
    
    B = ((-1)**n)!*((-1)**(mp-m+2*n))
    C = B/(Factorial(l-mp-n)*Factorial(l+m-n)*Factorial(mp-m+n)*Factorial(n))

    exp1 = 2*l + m - mp - 2*n
    exp2 = mp - m + 2*n

    tot = tot + C * (cosine**exp1) * (sine**exp2)    

  End Do

  WignerD = A * tot

End Function WignerD
!==========================================================================================
!==========================================================================================
Double Precision Function Factorial(a)

  Implicit None
  Integer, Intent(In) :: a
  
  Double Precision, Dimension(0:16) :: Fact

  If (a > 16) then
    Print *
    Print *
    Print *, '     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    Print *, '     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    Print *, '     !!!!!!!!  FACTORIAL TOO LARGE  !!!!!!!'
    Print *, '     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    Print *, '     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    Print *
    Print *
  End If

  Fact(0) = 1.0d0
  Fact(1) = 1.0d0
  Fact(2) = 2.0d0
  Fact(3) = 6.0d0
  Fact(4) = 24.0d0
  Fact(5) = 120.0d0
  Fact(6) = 720.0d0
  Fact(7) = 5040.0d0
  Fact(8) = 40320.0d0
  Fact(9) = 362880.0d0
  Fact(10) = 3628800.0d0 
  Fact(11) = 39916800.0d0
  Fact(12) = 479001600.0d0
  Fact(13) = 62270208000.d0
  Fact(14) = 87178291200.0d0
  Fact(15) = 1307674368000.0d0
  Fact(16) = 20922789888000.0d0

  Factorial = fact(a)
  
End Function Factorial
!==========================================================================================
!==========================================================================================
Double Precision Function Dot(v1,v2,len)
  
  Implicit NONE
  Integer, Intent(In) :: len
  
  Double Precision, Intent(In), Dimension(len) :: v1
  Double Precision, Intent(In), Dimension(len) :: v2

  Integer :: i 
  Double Precision :: sum

  sum = 0.0d0
  Do i = 1, len
    sum = sum + v1(i)*v2(i)
  End Do

  Dot = sum

End Function Dot 
!==========================================================================================
