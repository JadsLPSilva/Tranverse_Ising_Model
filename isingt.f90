!----------------------!
 module systemvariables
 implicit none

 integer :: ll                       ! system linear size
 integer :: nn                       ! number of spins (nn=ll*ll)
 character*80 outname		     ! name of the output file 
 integer, allocatable :: spin(:)     ! spins array

 integer :: bins
 integer(8) steps1,steps2	     ! # of bins, # of warm up steps, # of measurements steps 
 integer(8) nskips
 real*8 beta			     ! Temperature
 integer :: seed2		     ! seed

 integer(8) :: n_acep		     ! # of accepted flips, # of total measurements 

 real*8, allocatable :: energ(:) 
 real*8, allocatable :: energ2(:)   
 real*8, allocatable :: magnet(:),magnet2(:),magnet4(:) 
 real*8, allocatable :: ce(:), sus(:),U2(:)                                

 integer :: nb 
 
!--------------------transverse field-------------------------!
   integer :: trd                                             ! transverse dimension
   real*8 tn,jz,hx                                        ! tn : Trotter number; jz: campo em z; hx: campo em x. 
!-------------------------------------------------------------!
   
  end module systemvariables
!--------------------------!


!---------------!
 program isingt2d ! Main Program
!---------------!

 use systemvariables
 USE IFPORT
 implicit none

 integer(8) i
 integer :: j,loopf                !!!
 real*8 :: hmax,hmin,interv   !!! Loop
 character*80 magstring,mag2string,mag4string,susstring,binderstring
                                     !!!
 open(1,file='readisingt.in',status='old')
 read(1,2010) outname
 read(1,*) ll,trd,tn, jz
 read(1,*) hmax,hmin,interv


 read(1,*) steps1,steps2,bins
 read(1,*) nskips
 read(1,*) seed2
2010    format(A40)
 close(1)

! Your total number of steps should be a multiple of the number of bins
 if ( mod(steps2,bins) .ne. 0) then
 write(*,*) "BAD choice for bins and/or nskips!"
 stop
 end if

 !------------- transverse variables-----------------!
   !jz=1.d0
   beta=1.d0
   nn=ll*trd	!total number of spins > square lattice
   
   
   ! Writing partial output - part 1
   call output1()
   
  !------------------Open--------------!
  !------------------------------------!

    magstring = 'mag'//outname
    mag2string = 'mag2'//outname
    mag4string = 'mag4'//outname
    susstring = 'sus'//outname
    binderstring = 'binder'//outname

    open(unit=29,file=adjustl(adjustr(magstring)//'.dat'),status='replace')
    open(unit=30,file=adjustl(adjustr(mag2string)//'.dat'),status='replace')
    open(unit=31,file=adjustl(adjustr(mag4string)//'.dat'),status='replace')
    open(unit=32,file=adjustl(adjustr(susstring)//'.dat'),status='replace')
    open(unit=33,file=adjustl(adjustr(binderstring)//'.dat'),status='replace') 
   ! open(unit=34,file='ce.dat',status='replace')
  
 !-------------------------------------! 
    ! allocating arrays
    allocate (spin(0:nn-1))
    allocate(energ(0:bins-1)) 
    allocate(energ2(0:bins-1)) 
    allocate(magnet(0:bins-1))   
    allocate(magnet2(0:bins-1)) 
    allocate(magnet4(0:bins-1)) 
    allocate(sus(0:bins-1))     
    allocate(ce(0:bins-1))     
    allocate(U2(0:bins-1)) 
    
!--------------------------------------loop
   
   loopf=nint((hmax-hmin)/interv) + 1
   
   
   
   do j=1, loopf
   
   hx=hmin+(j-1)*interv
   
   
   write(*,*) " o programa está rodando p/ jz,hx,ll,trd,tn,beta=",  jz ,hx , ll, trd, tn, beta
   



! Reinitializating rand() and/or irand() functions
! call srand(seed2)

! random initial distribution for spins
 do i=0,nn-1   
  spin(i) = 2*mod(irand(),2) - 1 
 enddo

! Resetting the variables 
 call zeroas() 

!# loop de aquecimento
 do i=1,steps1
 call mc_steps()
 end do

! Writing partial output - part 2 
 write(11,*) "----------------------------"
 write(11,*) "Acceptance (warmup)=", dfloat(n_acep)/(dfloat(nn)*dfloat(steps1))
 write(11,*) "----------------------------"

! Resetting the variables (novamente para zerar o nÃºmero de aceitaÃ§Ã£o)
 call zeroas() 

! MC steps for measurements + subroutine for measurements
 do i=1,steps2

     	! Escrever na tela a porcentagem de passos jÃ¡ rodados
 	if (mod(i,steps2/100) .eq. 0)  write(*,*) "bin, percentage", nb, (i*100)/steps2, &
& dfloat(sum(spin))/dfloat(nn)  

	! MC steps
	call mc_steps()

	! Measurements (You will measure skipping 'nskips' steps)
	if ( mod(i,nskips) .eq. 0) then
	call measurements()
	end if

     	! Selecionando o bin		!  - new
	if( mod(i,(steps2/bins)) .eq. 0 ) then
	nb=nb+1
	end if

 enddo

!-------Susceptibility and Specific heat-------!


 sus(:)=(magnet2(:)-(magnet(:))**2)*(dfloat(nn)*beta)
 
 
 ce(:)=(energ2(:)-(energ(:))**2)*(dfloat(nn)*(beta)**2)

 U2(:)= 1.5d0-0.5d0*(magnet4(:)/(magnet2(:))**2)
!----------------------------------------------------------
 call output2()

 end do          !---------------------------------end loop

! deallocating arrays
 deallocate(spin)
 deallocate(energ)
 deallocate(energ2)
 deallocate(magnet)
 deallocate(magnet2)
 deallocate(magnet4)
 deallocate(sus)
 deallocate(ce)
 deallocate(U2)
 
! close
 
    close(11)
    close(29)
    close(30)
    close(31)
    close(32)
    close(33)
  !  close(34)
    
!-------------------!
 end program isingt2d ! Main Program
!-------------------!





!*************************************************************************
!****************************** SUBROUTINES ******************************
!*************************************************************************

!-------------------!
 subroutine output1()
!-------------------!
 use systemvariables
 implicit none

 character*80 acstring

 acstring = 'dados'//outname              
 open(unit=11,file=adjustl(adjustr(acstring)//'.dat'),position='append')

 write(11,*) "############################"
 write(11,*) "Ising 2D"
 write(11,*) "############################"
 write(11,*) " "
 write(11,*) "----------------------------"
 write(11,*) "Input data"
 write(11,*) "----------------------------"
 write(11,*) "Linear size=", ll , "por", trd
 write(11,*) "beta=", beta
 write(11,*) "Warming steps=", steps1
 write(11,*) "Measurements steps=", steps2
 write(11,*) "Number of bins=", bins
 write(11,*) "seed=", seed2

 return
 end subroutine output1
!---------------------!


!-------------------!
 subroutine output2()
!-------------------!
 use systemvariables
 implicit none

 real*8 aval, erval  


 write(11,*) "hx=", hx

 write(11,*) "----------------------------"
 write(11,*) "Acceptance (meas)=", dfloat(n_acep)/(dfloat(nn)*dfloat(steps2))
 write(11,*) "----------------------------"

 write(11,*) "----------------------------"

 call values1(energ,aval,erval,bins)
 write(11,*) "Energy_per_site=", aval, "+-", erval

!--------------------------------------------------------- 

!---------------------------------------------------------
 call values1(energ2,aval,erval,bins)
 write(11,*) "Energy2_per_site=", aval, "+-", erval 
!---------------------------------------------------------
 
!---------------------------------------------------------
 write(11,*) "----------------------------"
 write(11,*) "                            "
 write(11,*) "----------------------------"
 
 call values1(magnet,aval,erval,bins)
 write(11,*) "magnetization_per_site=", aval, "+-" , erval 
 write(29,*) hx , aval, erval
!---------------------------------------------------------------
 
!---------------------------------------------------------------
 
 call values1(magnet2,aval,erval,bins)
 write(11,*) "magnetization2_per_site=", aval, "+-" , erval 
 write(30,*) hx , aval, erval
!---------------------------------------------------------------- 
 
!----------------------------------------------------------------
 call values1(magnet4,aval,erval,bins)
 write(11,*) "magnetization4_per_site=", aval, "+-" , erval 
 write(31,*) hx , aval, erval
 write(11,*) "----------------------------"

!--------------------------------------------------------------- 

!----------------------------------------------------------------
 write(11,*) "----------------------------"
 write(11,*) "                            "
 write(11,*) "----------------------------"

 call values1(sus,aval,erval,bins)
 write(11,*) "susceptibility_per_site=", aval, "+-" , erval
 write(32,*) hx , aval, erval
!-----------------------------------------------------------------

!-----------------------------------------------------------------
 
!  write(11,*) "----------------------------"
!  
!  call values1(ce,aval,erval,bins)
!  write(11,*) "specific-heat_per_site=", aval, "+-" , erval 
!  write(34,*) hx , aval, erval
!-----------------------------------------------------------------

!-----------------------------------------------------------------

write(11,*) "----------------------------"
write(11,*) "----------------------------"

call values1(U2,aval,erval,bins)
 write(11,*) "binder-cumulante=", aval, "+-" , erval 
 write(33,*) hx , aval, erval
 write(11,*) "----------------------------"

 return
 end subroutine output2
!---------------------!


!-------------------!
 subroutine zeroas()  
!-------------------!
 use systemvariables
 implicit none

 n_acep=0
 nb=0		

 magnet=0.d0
 magnet2=0.d0
 magnet4=0.d0
 energ=0.d0
 energ2=0.d0
 sus=0.d0
 ce=0.d0
 U2=0.d0

 return
 end subroutine zeroas
!---------------------!




!---------------------! 
 subroutine mc_steps()
!--------------------!
   use systemvariables
   implicit none
   real*8 ran2
   integer :: ss, i
   integer :: linha,coluna
   real*8 :: s,sran
   integer :: s1,s2,s3,s4          
   real*8 deltU, heatb

   !##################
   !integer :: trd                       
   real*8 kx,ky
   !##################
   
   
  do i=1,nn !do_onestep

   112 continue 
   s= ran2(seed2)
   ss=nint(s*nn)
   if (ss==nn) go to 112

!----------------------------------------! definindo a variÃ¡vel linha e coluna:
   linha = mod(ss,ll)
                                         ! ok
   coluna= (ss - linha)/ll 

!-----------------------------------------------!

! Determinando os vizinhos do sitio escolhido : -seja s1,s2,s3,s4 estes vizinhos.
!! Periodic Boundary Conditions


!   call rede(linha,coluna,trd,s1,s2,s3,s4)
! 
     s1=spin(linha+mod(coluna+1,trd)*ll)     ! spin acima do sitio
     s2=spin(linha+mod(coluna-1+trd,trd)*ll)  ! spin abaixo do sitio
     s3=spin(mod(linha-1+ll,ll)+coluna*ll)  ! spin a direita do sitio
     s4=spin(mod(linha+1,ll)+coluna*ll)     ! spin a esquerda do sitio

!-----------------------------------------------!
    
      
    deltU = 2.d0*dfloat(spin(ss))*(tn*jz)*dfloat(s3+s4) 
    deltU=deltU + 2.d0*dfloat(spin(ss))*(-0.5d0*dlog(dtanh(tn*hx)))*dfloat(s1+s2)

   heatb= dexp(-beta*deltU) / ( 1.d0 + dexp(-beta*deltU) )
    
    sran=ran2(seed2)
 
   if ( sran .lt. heatb ) then 
   spin(ss)=-spin(ss)
   n_acep = n_acep + 1
   end if


  end do !do_onestep

  return
  end subroutine
!-----------------------------------------------!




!---------------------! 
 subroutine measurements()
!---------------------!
  use systemvariables
  implicit none
  real*8 :: E,E2,M,M2,M4
  integer :: H
  integer :: se,s1,s2,s3,s4
  integer(8) nmed
  integer :: linha,coluna
  integer :: i
  

!------bLoco de inicializacao das variaveis 
  nmed=steps2/nskips   ! numero de measurements
  E=0.d0                  ! E: energia do bloco de spins
  E2=0.d0                 !E2: quadrado da energia do bloco
  m=0.d0                  ! m: magnetizacao do bloco de spins
  m2=0.d0                 ! m2: quadrado da magnetizacao do bloco de spins
  m4=0.d0
 !------! bloco da energia !--------!

! ! !   do i=0,nn-1
! ! ! 
! ! !    se=spin(i)
! ! !                     !determinando os vizinhos
! ! !    linha = mod(i,ll)
! ! !    coluna= (i - linha)/ll
! ! ! 
! ! !     s1=spin(linha+mod(coluna+1,ll)*ll)     ! spin acima do sitio
! ! !     s2=spin(linha+mod(coluna-1+ll,ll)*ll)  ! spin abaixo do sitio
! ! !     s3=spin(mod(linha-1+ll,ll)+coluna*ll)  ! spin a direita do sitio
! ! !     s4=spin(mod(linha+1,ll)+coluna*ll)     ! spin a esquerda do sitio
! ! !                                          
! ! !   
! ! !    H = -se*(s1+s2+s3+s4)           ! Caso sem campo externo
! ! ! 
! ! !    E = E + 0.5d0*dfloat(H)/dfloat(nn) !0.5 pq a contagem acima dobrara a energia; dividir por nn Ã© pra termos densidade de energia 
! ! !   end do
! ! ! 
! ! ! ! Quadrado da energia:
! ! !    E2 = E**2
! ! ! 
! ! ! !------bloco calc. energia media
! ! !   
! ! !   energ(nb)=energ(nb) + (E/dfloat(nmed/bins)) !  > Notar que temos que dividir pelo numero de bins tambem
! ! !   energ2(nb)=energ2(nb) + (E2/dfloat(nmed/bins)) 

!-------! bloco da Magnetizacao !---------!

 m = dabs(dfloat(sum(spin)))/dfloat(nn)     
 m2= m**2
 m4= m**4

!-------! calculo da magnetizacao( m , m**2 , m**4 )
 magnet(nb)=magnet(nb) + m/dfloat(nmed/bins)
 magnet2(nb)=magnet2(nb) + m2/dfloat(nmed/bins)
 magnet4(nb)=magnet4(nb) + m4/dfloat(nmed/bins)


  return
  end subroutine
!-----------------------------------------------!


!-----------------------ERROR_BAR-----------------------------!
 subroutine values1(vm,aval,erval,bins)
!-------------------------------------------------------------!
 implicit none
 integer bins
 real*8 aval, erval
 real*8 vm(0:bins-1), vr(bins)

 aval= sum(vm(:))/dfloat(bins)		!mean value

 vr(:) = (vm(:) - aval)**2		!standard deviation

 erval = dsqrt( sum(vr)/dfloat(bins - 1) )	!error

 return
 end subroutine values1
!-------------------------------------------------------------!



!*****************************************************************
      DOUBLE PRECISION FUNCTION RAN2(IDUM)
      implicit none
      double precision rm
      integer ia,ic,idum,iff,ir,iy,j,m
      save 
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112d-6)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M) 
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
!     IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END
