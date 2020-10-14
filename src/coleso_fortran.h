       MODULE COLESO_FORTRAN
       INTERFACE 

       SUBROUTINE coleso_add_function(FUNCNAME, ID) 
     *   BIND(C, NAME='coleso_add_function') 
         USE, INTRINSIC :: ISO_C_BINDING  
         IMPLICIT NONE 
         CHARACTER, DIMENSION(*) :: FUNCNAME
         INTEGER (C_INT) :: ID
       END SUBROUTINE coleso_add_function 

       SUBROUTINE coleso_read_file(ID, FILENAME) 
     *   BIND(C, NAME='coleso_read_file') 
         USE, INTRINSIC :: ISO_C_BINDING  
         IMPLICIT NONE 
         INTEGER (C_INT), VALUE :: ID
         CHARACTER, DIMENSION(*) :: FILENAME    
       END SUBROUTINE coleso_read_file

       SUBROUTINE coleso_set_parameter(PARAMNAME, PARAMVALUE) 
     *   BIND(C, NAME='coleso_set_parameter') 
         USE, INTRINSIC :: ISO_C_BINDING  
         IMPLICIT NONE 
         CHARACTER, DIMENSION(*) :: PARAMNAME, PARAMVALUE
       END SUBROUTINE coleso_set_parameter

       SUBROUTINE coleso_read_set(ID) 
     *   BIND(C, NAME='coleso_read_set') 
         USE, INTRINSIC :: ISO_C_BINDING  
         IMPLICIT NONE 
         INTEGER (C_INT), VALUE :: ID
       END SUBROUTINE coleso_read_set

       SUBROUTINE coleso_init(ID)
     *   BIND(C, NAME='coleso_init') 
         USE, INTRINSIC :: ISO_C_BINDING  
         IMPLICIT NONE 
         INTEGER (C_INT), VALUE :: ID
       END SUBROUTINE coleso_init

       SUBROUTINE coleso_pointvalue(I, T, C, V)
     *   BIND(C, NAME='coleso_pointvalue') 
         USE, INTRINSIC :: ISO_C_BINDING  
         IMPLICIT NONE 
         INTEGER (C_INT), VALUE :: I 
         REAL*8, VALUE :: T
         REAL*8, DIMENSION(*) :: C, V 
       END SUBROUTINE coleso_pointvalue 

       END INTERFACE 
       END MODULE COLESO_FORTRAN
