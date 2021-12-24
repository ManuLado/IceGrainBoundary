program bordedegrano
  implicit none
  integer,allocatable, dimension(:,:)::M,P
  real(kind(1.d0)),allocatable,dimension(:)::vector_A,y,mmult,mmult2,y0
  real(kind(1.d0))::E,A, B, h, dt, xmax,tangbeta,x0,t,ymax,ymin, tt, depth
  integer::i, j,N, k, imax,superpaso,paso
  character(len=60)::archivo,tiempo,directorio,karakt,chain2,chain3
  character,allocatable,dimension(:)::chain
  B=1.52E-19 !coeficiente de Difusion (cm⁴/s)
  E=0*(-6.22E-6) !coef de Evaporacion libre(cm/s)
  A=6.59E-12 !coef.de Evaporacion-condensacion (cm/s)
  dt=0.001 !intervalo de tiempo (s)
  xmax=50.0E-4 !intervalo en x desde 0 (cm)
  tangbeta=0.254 !tang. angulo diahedral
  N=250 !Dimension total (define precision del algoritmo)
  superpaso=1000000 !Define cuantas iteraciones totales hay
  paso=superpaso/100 !Define cuantas curvas (y,x) va a graficar
  !------------------------------------------------------------------------------
  open(1, file="plot.plt", status="replace")!Comandos Gnuplot para graficar las
  curvas (y,x)
  write(1,*) "# plot.plt"
  write(1,*) 'set xlabel "x"'
  write(1,*) 'set yrange [-0.000025:0.000016]'
  write(1,*) 'set ylabel "y"'
  write(1,*) 'set xrange [0:0.001]'
  write(1,*) "set terminal png"
  write(1,*) "set nokey"
  write(1,*) "set grid"
  do k=1,superpaso
  if (mod(k,paso)==0) then
  write(karakt,'(I10)') k
  write(1,*) 'set title "perfil t=',trim(adjustl(karakt)),'.dt s"'
  write(chain2,*) trim('"perfil'//trim(adjustl(karakt))//'pasos' //
  '.dat"'//' with linespoints')
  write(1,*)'set output ' ,'"',trim(adjustl(karakt)),'.png"'
  write(1,*) "plot ",chain2
  end if
  end do
  !------------------------------------------------------------------------------
  allocate
  (vector_A(1:N),EE(1:N),M(1:N,1:N),P(1:N,1:N),Y(1:N),mmult(1:N),mmult2(1:N),y0(1
  :N))
  x0=0.d0
  h=(xmax-x0)/real(N,kind(1.d0))
  5 format(2(E19.13))
  6 format(3(1X,F45.13))
  7 format(100(I2))
  open(2, file="profundidad.dat", status="replace")
  open(3, file="semiancho.dat", status="replace")
  write(2,*) "m(Bt)1/4 profundidad "write(3,*) "(Bt)1/4 semiancho ymax "
  !----------------------- matriz M(i,j) (de dif.superficial) i=columnas, j=filas
  M=0
  M(1,1)=28
  M(1,2)=-36
  M(1,3)=12
  M(2,1)=-39
  M(2,2)=64
  M(2,3)=-40
  M(3,1)=12
  M(3,2)=-39
  M(3,3)=56
  M(4,1)=-1
  M(4,2)=12
  M(4,3)=-39
  M(5,1)=0
  M(5,2)=-1
  M(5,3)=12
  M(6,3)=-1
  do i=4,N-3
  do j=4,N-3
  if (i==j) then
  M(i,j)=56
  M(i+1,j)=-39
  M(i+2,j)=12
  M(i+3,j)=-1
  M(i-3,j)=-1
  M(i-2,j)=12
  M(i-1,j)=-39
  end if
  end do
  end do
  M(N,N)=28
  M(N,N-1)=-28
  M(N,N-2)=11
  M(N-1,N)=-39
  M(N-1,N-1)=56
  M(N-1,N-2)=-39
  M(N-2,N)=12
  M(N-2,N-1)=-39
  M(N-2,N-2)=56
  M(N-3,N)=-1
  M(N-3,N-1)=12
  M(N-3,N-2)=-39
  M(N-4,N)=0
  M(N-4,N-1)=-1
  M(N-4,N-2)=12
  M(N-5,N)=0
  M(N-5,N-1)=0
  M(N-5,N-2)=-1
  M(N-6,N)=0
  M(N-6,N-1)=0
  M(N-6,N-2)=0
  !**************************************Defino el vector vector_A
  vector_A=E
  vector_A(1)=E-3.d0*tangbeta*B/(h**3)-(26.d0/12.d0)*tangbeta*A/(h)
  vector_A(2)=E+3.d0*tangbeta*B/(h**3)+(1.d0/6.d0)*tangbeta*A/(h)
  vector_A(3)=E-(2.d0/6.d0)*tangbeta*B/(h**3)
  !************************************* DEFINO la matriz P (de evapo.equilibr)
  P=0P(1,1)=-27
  P(1,2)=16
  P(1,3)=-1
  P(2,1)=28
  P(2,2)=-31
  P(2,3)=16
  P(3,1)=-1
  P(3,2)=16
  P(3,3)=-30
  P(4,1)=0
  P(4,2)=-1
  P(4,3)=16
  P(5,1)=0
  P(5,2)=0
  P(5,3)=-1
  do i=4,N-2
  do j=4,N-2
  if (i==j) then
  P(i,j)=-30
  P(i+1,j)=16
  P(i+2,j)=-1
  P(i-2,j)=-1
  P(i-1,j)=16
  end if
  end do
  end do
  P(N,N)=-15
  P(N,N-1)=15
  P(N,N-2)=-1
  P(N-1,N)=16
  P(N-1,N-1)=-30
  P(N-1,N-2)=16
  P(N-2,N)=-1
  P(N-2,N-1)=16
  P(N-2,N-2)=-30
  P(N-3,N)=0
  P(N-3,N-1)=-1
  P(N-3,N-2)=16
  P(N-4,N)=0
  P(N-4,N-1)=0
  P(N-4,N-2)=-1
  !**************************************Defino el vector y
  y=0.d0
  t=0.
  ymax=0.d0
  ymin=0.d0
  imax=0.d0
  do k=0, superpaso !cuantas curvas quiero de t
  mmult=0.d0
  mmult2=0.d0
  y0=y
  x0=0.d0
  do i=1,N !producto matriz por vector
  do j=1,N
  mmult(i)=mmult(i)+real(M(i,j),kind(1.d0))*y0(j)
  mmult2(i)=mmult2(i)+real(P(i,j),kind(1.d0))*y0(j)
  end doy(i)=y0(i)+(A*mmult2(i)/(12.d0*(h**2))+vector_A(i)-
  B*mmult(i)/(6.d0*(h**4)))*dt
  end do
  ymax=maxval(y)
  ymin=minval(y)
  imax=maxloc(y,1) !calculo del ancho
  depth=abs(ymax-ymin) !calculo de la profundidad
  t=t+dt
  tt=(B*t)**0.25
  x0=0.d0
  !----------------------------------- Comandos para que escriba los datos
  ! en archivos con nombre "perfil-paso.dat"
  if (mod(k,paso)==0) then
  Write(tiempo,'(i10)') k
  write(*,*) k
  archivo=trim( 'perfil' //trim(adjustl(tiempo))// 'pasos' // '.dat')
  open(unit=30,file=archivo,status="replace")
  do i=1,N
  write(30,6)x0, y(i)
  x0=x0+h
  end do
  close(30)
  do i=1,N
  x0=x0+h
  if (i==imax) then
  write(3,*)tt,x0, y(imax)
  write(2,*)tt*tangbeta, depth
  end if
  end do
  end if
  end do
  close (2)
  close(3)
  !*************************************
  call system('gnuplot -p plot.plt')
end program bordedegrano
