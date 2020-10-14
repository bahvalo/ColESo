       include "coleso_fortran.h"

       subroutine solution_using_prescribed_parameters
       USE COLESO_FORTRAN
       USE, INTRINSIC :: ISO_C_BINDING

!      string buffers
       CHARACTER(LEN=100, KIND=C_CHAR) :: NAME, VALUE
!      handle for the solution
       integer ID
!      time and coordinates
       real*8 t
       real*8 c(3)
!      buffer for the solution
       real*8 v(5)

!      Creating an object for the solution
       NAME = '4peak' // C_NULL_CHAR
       call coleso_add_function(NAME, ID)

!      Setting solution parameters
!      Xmin and Xmax prescribes a linear coordinate mapping
       NAME = 'Xmin' // C_NULL_CHAR
       write(VALUE,*) 0.0, C_NULL_CHAR
       call coleso_set_parameter(NAME, VALUE)
       NAME = 'Xmax' // C_NULL_CHAR
       write(VALUE,*) 1.0, C_NULL_CHAR
       call coleso_set_parameter(NAME, VALUE)
!      Set additional space of size 0.5 between pulses
       NAME = 'Period' // C_NULL_CHAR
       write(VALUE,*) 1.5, C_NULL_CHAR
       call coleso_set_parameter(NAME, VALUE)

!      Applying parameters to the solution object
       call coleso_read_set(ID)

!      Initializing of the solution
!      (for this particular case this call has no effect)
       call coleso_init(ID)

!      Calculating the solution and printing it to test1.dat
       open(1, file='test1.dat')
       t = 0.4;
       c(1) = 0d0
       c(2) = 0d0
       c(3) = 0d0
       do i = 0, 384
           c(1) = i / 128d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) c(1), v(1)
       enddo
       close(1)
       end subroutine solution_using_prescribed_parameters


       subroutine solution_using_parameters_in_file
       USE COLESO_FORTRAN
       USE, INTRINSIC :: ISO_C_BINDING

!      string buffers
       CHARACTER(LEN=100, KIND=C_CHAR) :: NAME, VALUE
!      handle for the solution
       integer ID
!      time and coordinates
       real*8 t
       real*8 c(3)
!      buffer for the solution
       real*8 v(5)

!      First write the parameters to a file
       open(1, file = 'test.txt')
!      Xmin and Xmax prescribes a linear coordinate mapping
       write(1,*) 'Xmin 0.0'
       write(1,*) 'Xmax 1.0'
!      Set additional space of size 0.5 between pulses
       write(1,*) 'Period 1.5'
       close(1)

!      Creating an object for the solution
       NAME = '4peak' // C_NULL_CHAR
       call coleso_add_function(NAME, ID)

!      Reading solution parameters from a file
       NAME = 'test.txt' // C_NULL_CHAR
       call coleso_read_file(ID, NAME)

!      Initializing of the solution
!      (for this particular case this call has no effect)
       call coleso_init(ID)

!      Calculating the solution and printing it to test1.dat
       open(1, file='test2.dat')
       t = 0.4;
       c(1) = 0d0
       c(2) = 0d0
       c(3) = 0d0
       do i = 0, 384
           c(1) = i / 128d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) c(1), v(1)
       enddo
       close(1)
       end subroutine solution_using_parameters_in_file

       PROGRAM MAIN
       CALL solution_using_prescribed_parameters
       CALL solution_using_parameters_in_file
       END
