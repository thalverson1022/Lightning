Module Inputs
  Implicit NONE

  Character(len = 128) :: Big_FN
  Character(len = 128) :: Geos_FN
  Character(len = 128) :: Params_FN
  Character(len = 128) :: E_FN
  Character(len = 128) :: Vals_FN
  Character(len = 128) :: Vecs_FN
  
 
  Integer :: LMax         
  Integer :: Dmax

  Double Precision :: B_const
  Double Precision :: mu_const


  Character(len = 1) :: trunc_type
  Character(len = 1) :: parity_type

End Module Inputs

Module Globals
  Implicit NONE
  
  Integer :: Big_Len
  Integer :: Small_Len
  Integer :: E_Len
  
  Integer, Allocatable, Dimension(:,:) :: Big_LM
  Integer, Allocatable, Dimension(:,:) :: Small_LM

  Integer, Allocatable, Dimension(:,:,:,:) :: IndexMap

  Double Precision, Allocatable, Dimension(:) :: Geos
  Double Precision, Allocatable, Dimension(:,:) :: Rijs
  Double Precision, Allocatable, Dimension(:,:) :: V12
  Double Precision, Allocatable, Dimension(:,:) :: E

  Double Complex, Allocatable, Dimension(:,:,:,:) :: Wigners
  Double Complex, Allocatable, Dimension(:,:,:,:) :: Wigners_Star
  
End Module Globals

Module Outputs
  Implicit NONE
  Double Complex, Allocatable, Dimension(:,:)  :: H
  Double Precision, Allocatable, Dimension(:)  :: EigenVals
End Module Outputs