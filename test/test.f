       include "coleso_fortran.h"

       program main
       USE COLESO_FORTRAN
       USE, INTRINSIC :: ISO_C_BINDING
       CHARACTER(LEN=100, KIND=C_CHAR) :: NAME, VALUE
       
       integer ID, i, j
       real*8 t, x, y, z, r, phi, C(3), V(10), V1(10), V2(10), V3(10)
       real*8 PiNumber
       
       PiNumber = datan(1d0) * 4d0

       write(*,*) 'DATA1D/4peak.dat'
       call coleso_add_function('4peak' // C_NULL_CHAR, ID)
       call coleso_read_file(ID, 'PARAMS/es_4peak.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/4peak.dat')
       do i = -1000, 1000
           t = 0d0; C(1) = i / 1000d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1)
       enddo
       close(1)
       
       write(*,*) 'DATA1D/planargauss.dat'
       call coleso_add_function('planargauss' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_planargauss.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/planargauss.dat')
       do i = 4300, 5700
           t = 25d0; C(1) = i / 100d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1)
       enddo
       close(1)       

       write(*,*) 'DATA1D/planarsinus_1.dat'
       call coleso_add_function('planarsinus' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_planarsinus_1.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/planarsinus_1.dat')
       do i = 0, 1000
           t = 50d0; C(1) = i / 10d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1)
       enddo
       close(1)       
       
       write(*,*) 'DATA2D/planarsinus_2.dat'
       call coleso_add_function('planarsinus' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_planarsinus_2.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/planarsinus_2.dat')
       do i = 0, 50
       do j = 0, 200
           t = 10d0; C(1) = j * 0.02d0; C(2) = i * 0.02d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/entropyvortex_1.dat'
       call coleso_add_function('entropyvortex' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_entropyvortex.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/entropyvortex_1.dat')
       do i = 370, 640, 5
       do j = 410, 680, 5
           t = 5d0; C(1) = j * 0.1d0; C(2) = i * 0.1d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(2)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/entropyvortex_2.dat'
       call coleso_add_function('entropyvortex' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_entropyvortex.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/entropyvortex_2.dat')
       do i = 370, 640, 5
       do j = 410, 680, 5
           t = 5d0; C(1) = j * 0.1d0; C(2) = i * 0.1d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA1D/source1d_1.dat'
       call coleso_add_function('source1d' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_source1d_1.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/source1d_1.dat')
       do i = 0, 1000
           t = 20d0; C(1) = i / 10d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1), V(2)
       enddo
       close(1)       

       write(*,*) 'DATA1D/source1d_2.dat'
       call coleso_add_function('source1d' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_source1d_2.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/source1d_2.dat')
       do i = 0, 1000
           t = 20d0; C(1) = i / 10d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1), V(2)
       enddo
       close(1)       

       write(*,*) 'DATA1D/acousticshock.dat'
       call coleso_add_function('acousticshock' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_acousticshock.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/acousticshock.dat')
       do i = 0, 460
           C(1) = i; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, 0d0, c, v1)
           call coleso_pointvalue(ID, 35d0, c, v2)
           call coleso_pointvalue(ID, 105d0, c, v3)
           write(1,*) C(1), V1(1), V2(1), V3(1)
       enddo
       close(1)       

       write(*,*) 'DATA2D/gaussian2d.dat'
       call coleso_add_function('gaussian2d' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_gaussian2d.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/gaussian2d.dat')
       do i = 0, 340, 5
       do j = 0, 400, 5
           t = 30d0; C(1) = j * 0.1d0; C(2) = i * 0.1d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/source2d.dat'
       call coleso_add_function('source2d' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_source2d.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/source2d.dat')
       do i = 0, 100
       do j = 0, 100
           t = 12d0; C(1) = i * 1d0; C(2) = j; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/gaussian3d.dat'
       call coleso_add_function('gaussian3d' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_gaussian3d.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/gaussian3d.dat')
       do i = 0, 340, 5
       do j = 0, 400, 5
           t = 30d0; C(1) = j * 0.1d0; C(2) = i * 0.1d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/source3d.dat'
       call coleso_add_function('source3d' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_source3d.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/source3d.dat')
       do i = 0, 500, 15
       do j = -250, 1000, 15
           t = 25d0; C(1) = j * 0.1d0; C(2) = i * 0.1d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(2)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/pointsource.dat'
       call coleso_add_function('pointsource' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_pointsource.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/pointsource.dat')
       do i = 0, 150, 2
       do j = 0, 300, 2
           t = 50d0; C(1) = j; C(2) = i; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(2)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/rotatingdipole.dat'
       call coleso_add_function('rotatingdipole' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_rotatingdipole.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/rotatingdipole.dat')
       do i = -100, 100
       do j = -100, 100
           t = 0d0; C(1) = j * 0.01d0; C(2) = i * 0.01d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(5)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/coaxial_1.dat'
       call coleso_add_function('coaxial' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_coaxial_1.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/coaxial_1.dat')
       do i = 30, 100
       do j = -100, 100
           r = i * 0.01; phi = j * PiNumber / 100d0
           t = 1d0; C(1) = r * dcos(phi); C(2) = r * dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/coaxial_2.dat'
       call coleso_add_function('coaxial' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_coaxial_2.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/coaxial_2.dat')
       do i = 10, 110
       do j = 0, 50
           t = 17d0; C(1) = i; C(2) = 0d0; C(3) = j
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(3), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/corner_1.dat'
       call coleso_add_function('corner' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_corner_1.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/corner_1.dat')
       do i = -42, 66, 2
       do j = -85, 110, 2
           t = 69d0; C(1) = j; C(2) = i; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/corner_2.dat'
       call coleso_add_function('corner' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_corner_2.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/corner_2.dat')
       do i = 0, 100
       do j = 0, 100
           r = i * 0.01; phi = j * (PiNumber * 2d0 / 3d0) / 100d0
           t = 0.9d0; C(1) = r*dcos(phi); C(2) = r*dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/cornerplanar_1.dat'
       call coleso_add_function('cornerplanar' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_cornerplanar_1.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/cornerplanar_1.dat')
       do i = -130, 160, 5
       do j = -250, 200, 5
           t = 60d0; C(1) = j / 10d0; C(2) = i / 10d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/cornerplanar_2.dat'
       call coleso_add_function('cornerplanar' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_cornerplanar_2.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/cornerplanar_2.dat')
       do i = 0, 100
       do j = 0, 100
           r = i * 0.01; phi = j * (PiNumber * 2d0 / 3d0) / 100d0
           t = 0.4d0; C(1) = r*dcos(phi); C(2) = r*dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/cylinder.dat'
       call coleso_add_function('cylinder' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_cylinder.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/cylinder.dat')
       do i = 5, 100
       do j = 0, 200
           r = i * 0.1; phi = j * (PiNumber) / 200d0
           t = 8d0; C(1) = r*dcos(phi); C(2) = r*dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/sinusvisc.dat'
       call coleso_add_function('sinusvisc' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_sinusvisc.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/sinusvisc.dat')
       do i = 0, 100
       do j = 0, 100
           t = 0d0; C(1) = j / 100d0; C(2) = i / 100d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(1)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA1D/vortexincylinder.dat'
       call coleso_add_function('vortexincylinder' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_vortexincylinder.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/vortexincylinder.dat')
       do i = 0, 100
           C(1) = i / 100d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, 0d0, c, v1)
           call coleso_pointvalue(ID, 0.1d0, c, v2)
           write(1,*) C(1), V1(3), V2(3)
       enddo
       close(1)       

       write(*,*) 'DATA2D/waveinchannel_1.dat'
       call coleso_add_function('waveinchannel' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_waveinchannel_1.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/waveinchannel_1.dat')
       do i = -25, 25
       y = 0d0
       do j = 0, 10000
          if(y > 0.5d0) cycle
          t = 0.5d0; C(1) = i / 100d0; C(2) = y; C(3) = 0d0
          call coleso_pointvalue(ID, t, c, v)
          write(1,*) C(1), C(2), V(2)
          y = y + 0.02*(0.5-y)+0.001
       enddo
       write(1,*)
       enddo
       close(1)  
       
       write(*,*) 'DATA2D/waveinchannel_2.dat'
       call coleso_add_function('waveinchannel' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_waveinchannel_2.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/waveinchannel_2.dat')
       do i = 0, 50
       do j = -25, 25
           t = 1d0; C(1) = j / 100d0; C(2) = i / 100d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(2)
       enddo
       write(1,*)
       enddo
       close(1)       
        
       write(*,*) 'DATA1D/waveinchannel_3.dat'
       call coleso_add_function('waveinchannel' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_waveinchannel_3.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/waveinchannel_3.dat')
       do i = 0, 1000
           t = 0.1d0; C(1) = 0d0; C(2) = i / 1000d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(2), V(1), V(3), V(5)
       enddo
       close(1)
       
       write(*,*) 'DATA2D/waveinchannel_4.dat'
       call coleso_add_function('waveinchannel' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_waveinchannel_4.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/waveinchannel_4.dat')
       do i = 0, 100
       r = 0d0
       do j = 0, 10000
          if(r > 1d0) cycle
          phi = PiNumber / 4d0 * (i / 100d0)
          t = 0.0d0; C(1)=r*dcos(phi); C(2)=r*dsin(phi); C(3) = 0d0
          call coleso_pointvalue(ID, t, c, v)
          write(1,*) C(1), C(2), V(1)
          r = r + 0.01*(1.0-r)+1d-5
       enddo
       write(1,*)
       enddo
       close(1) 

       write(*,*) 'DATA1D/vortexes.dat'
       call coleso_add_function('finitevortex' // C_NULL_CHAR, ID1)
       call coleso_add_function('gaussianvortex' // C_NULL_CHAR, ID2)
       call coleso_add_function('rankinevortex' // C_NULL_CHAR, ID3)
       call coleso_read_file
     1 (ID1, 'PARAMS/es_finitevortex.txt' // C_NULL_CHAR)
       call coleso_read_file
     1 (ID2, 'PARAMS/es_gaussianvortex.txt' // C_NULL_CHAR)
       call coleso_read_file
     1 (ID3, 'PARAMS/es_rankinevortex.txt' // C_NULL_CHAR)
       call coleso_init(ID1)
       call coleso_init(ID2)
       call coleso_init(ID3)
       open(1, file='DATA1D/vortexes.dat')
       do i = 0, 80
           t = 0d0; C(1) = i / 10d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID1, t, c, v1)
           call coleso_pointvalue(ID2, t, c, v2)
           call coleso_pointvalue(ID3, t, c, v3)
           write(1,*) C(1), V1(3), V2(3), V3(3), V1(5), V2(5), V3(5)
       enddo
       close(1)       

       write(*,*) 'DATA1D/riemann.dat'
       call coleso_add_function('riemann' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_riemann.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/riemann.dat')
       do i = 175, 310
           t = 4d0; C(1) = i; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1)
       enddo
       close(1)  
       
       write(*,*) 'DATA1D/simplewave.dat'
       call coleso_add_function('simplewave' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_simplewave.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/simplewave.dat')
       do i = -70, 50
           C(1) = i / 100d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, 0d0, c, v1)
           call coleso_pointvalue(ID, 0.1d0, c, v2)
           write(1,*) C(1), V1(1), V1(5), V2(1), V2(5)
       enddo
       close(1)  
       
       write(*,*) 'DATA1D/viscshock_1.dat'
       call coleso_add_function('viscshock' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_viscshock_1.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/viscshock_1.dat')
       do i = -200, 200
           t = 0d0; C(1) = i / 10d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1), V(2), V(5)
       enddo
       close(1)  

       write(*,*) 'DATA1D/viscshock_2.dat'
       call coleso_add_function('viscshock' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_viscshock_2.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/viscshock_2.dat')
       do i = -200, 200
           t = 0d0; C(1) = i / 10d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1), V(2), V(5)
       enddo
       close(1)  

       write(*,*) 'DATA1D/viscshock_3.dat'
       call coleso_add_function('viscshock' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_viscshock_3.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA1D/viscshock_3.dat')
       do i = -200, 200
           t = 0d0; C(1) = i / 10d0; C(2) = 0d0; C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), V(1), V(2), V(5)
       enddo
       close(1)  
       
       write(*,*) 'DATA2D/conccyl.dat'
       call coleso_add_function('conccyl' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_conccyl.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/conccyl.dat')
       do i = 100, 200
       do j = 0, 50
           r = i * 0.01; phi = PiNumber * (0.25d0 + j*0.02d0/3d0)
           t = 0d0; C(1) = r * dcos(phi); C(2) = r * dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(5)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/curlfreecylinder.dat'
       call coleso_add_function('curlfreecylinder' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_curlfreecylinder.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/curlfreecylinder.dat')
       do i = 50, 500
       do j = 0, 50
           r = i * 0.01; phi = j * PiNumber / 50d0
           t = 0d0; C(1) = r * dcos(phi); C(2) = r * dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(5)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/potentialsphere.dat'
       call coleso_add_function('potentialsphere' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_potentialsphere.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/potentialsphere.dat')
       do i = 50, 500
       do j = 0, 50
           r = i * 0.025; phi = j * PiNumber / 50d0
           t = 0d0; C(1) = r * dcos(phi); C(2) = r * dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(5)
       enddo
       write(1,*)
       enddo
       close(1)       

       write(*,*) 'DATA2D/viscsphere.dat'
       call coleso_add_function('viscsphere' // C_NULL_CHAR, ID)
       call coleso_read_file
     1 (ID, 'PARAMS/es_viscsphere.txt' // C_NULL_CHAR)
       call coleso_init(ID)
       open(1, file='DATA2D/viscsphere.dat')
       do i = 50, 500
       do j = 0, 50
           r = i * 0.025; phi = j * PiNumber / 50d0
           t = 0d0; C(1) = r * dcos(phi); C(2) = r * dsin(phi); C(3) = 0d0
           call coleso_pointvalue(ID, t, c, v)
           write(1,*) C(1), C(2), V(5)
       enddo
       write(1,*)
       enddo
       close(1)       

       end
