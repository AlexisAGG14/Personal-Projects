
! ======================================================================
       module  parameteFuerza
          real (kind=8) :: cut2,Rcero,Kfe,kcon,Fw
          real (kind=8) :: Rcu,Altura_box
          Integer (kind=8) :: iprint
          Integer (kind=8) :: Tequili
        end module parameteFuerza


! ======================================================================
       module coor
          integer, parameter :: npmax=5000
          real (kind=8) :: x(npmax),y(npmax),z(npmax)
          real (kind=8) :: q,s,ps,AL
          Integer (kind=8):: npart,G

!     x(i),y(i),z(i)    : position particle i
!     npart             : actual number of particles
!     q                 : Nose Hoover Mass
!     g                 : number of degrees of freedom
!     s,ps              : Nose-Hoover coupling parameter
        end module coor


! ======================================================================
        module  veloc
          integer, parameter :: npmax=5000
          real (kind=8) :: vx(npmax),vy(npmax),vz(npmax)
          integer (kind=8) :: itt
!     vx(i),vy(i),vz(i)    : velocity particle i
        end module veloc


! ======================================================================
! ============ Program para polimero ===================================
! ======================================================================




      program Poly
        use parameteFuerza  ! parametros para el potencial.
        use coor            ! Coordenadas de posicion.
        use veloc           ! Coordenadas de velocidad.
        
        implicit none
        integer (kind=8) :: nuno, ndos, ntre, ncua, ncin !#####


!  --  Definicion de tipo de variable que se usan (real,entera, etc):
        INCLUDE 'ArchivoDefin.inc'
        real (kind=8) :: Lx,Ly,Lz,Lxcm,Lycm,Lzcm             !  - C. Momento angular -
        real (kind=8) :: sumrx,sumry,sumrz,sumvx,sumvy,sumvz !  - C. Momento angular -
        integer (kind=8) :: iTem
        real (kind=8) :: pp   ! contador par saber cuantos pasos se interan
       
!  --  Nombre de los archivos donde se colocan resultados:
        INCLUDE 'ArchivoOpen.inc'
        

!  -- VAlor de los parametros usados durante la simulacion:
        INCLUDE 'ArchivoValores.inc'

!  -- VAlor del parametro que se esta variando: (TEMPERATURA)
        !INCLUDE 'sig.inc'


        INCLUDE 'valor_corte.inc'  !#####


! Condicion para deterner el programa, por estructura del codigo.
         If (n.GT. npmaxx) Then
             Write (0,*) 'Numero de pariculas muy grandote '
             stop
         EndIf

! coordenadas iniciales del polimero- cof. incial lineal
        Call inicio(n)

! Coordenadas inciales de velocidad con momento lineal cero.
         Call vel_init_poli(n,sig,sumv2)

! -------- C. Momento angular: inicial --------
        sumrx = 0.d0
        sumry = 0.d0
        sumrz = 0.d0

        sumvx = 0.d0
        sumvy = 0.d0
        sumvz = 0.d0
        
        Lxcm = 0.d0
        Lycm = 0.d0
        Lzcm = 0.d0
 
       do i=1, n
        sumrx = sumrx + X(i)
        sumry = sumry + Y(i)
        sumrz = sumrz + Z(i)

        sumvx = sumvx + vX(i)
        sumvy = sumvy + vY(i)
        sumvz = sumvz + vZ(i)
       enddo

        sumrx = sumrx/n
        sumry = sumry/n
        sumrz = sumrz/n

        sumvx = sumvx/n
        sumvy = sumvy/n
        sumvz = sumvz/n
        
          Lxcm = Lxcm + sumry*n*sumvz-sumvy*n*sumrz
          Lycm = Lycm + sumrz*n*sumvx-sumvz*n*sumrx
          Lzcm = Lzcm + sumrx*n*sumvy-sumvx*n*sumry

        Lx = 0.d0
        Ly = 0.d0
        Lz = 0.d0
         do i=1,n
          Lx = Lx + y(i)*vz(i)-vy(i)*z(i)
          Ly = Ly + z(i)*vx(i)-vz(i)*x(i)
          Lz = Lz + x(i)*vy(i)-vx(i)*y(i)
        enddo

        Write(45,*)Lx,Ly,Lz,Lxcm,Lycm,Lzcm
! --------Fin: C. Momento angular --------

                   ! Imprime datos en archvo 0_Salida.dat
       Write (0,*) 
       Write (0,*)'Imprime datos cada,total datos:',iprint,npasos/iprint
       Write (0,*)
       Write (0,*)'Temperatura incial:Tem,sumv2',sumv2/(3.d0*n),sumv2
       Write (0,*)
       Write (0,*)'Tiempo de simulacion:',SNGL(npasos*dt)
       Write (0,*)


          npart = n
          G =  3*npart ! Valores iniciales de parametros de integracion -NH
          S = 0
          PS = (sumv2-G*sig)/Q


          ekin = sumv2*0.5d0

! Parametros de fluctuacion: Inicializa : muestreo en cero.
         i=1
         Call SAMPLE(npart,0,i,En,ekin,hami)


!----------- Fuerzas Totales: primer paso de interacion-----------
          Call Fuerzas(i,n,fx,fy,fz,epot)






!########################## termo n=0  #################################
!- - - - - - Inicia el proceso de interacion - - - - - - - - - - - - - -
       Do k = 1, 100000

          CALL SOLVE(1, fx, fy, fz, ekin, dt, 123456, sig, epot, hami)

           
          Call Fuerzas(k,n,fx,fy,fz,epot)
      

          CALL SOLVE(2, fx, fy, fz, ekin, dt, 123456, sig, epot, hami)

           
!.----------------------------------------------------------------------
!.----------------------------------------------------------------------
      if(mod(k,20000).EQ.0) then
         
         WRITE (43, *) k, SNGL(2*ekin/DBLE(3*N)),SNGL(Hami)

      endif
!.----------------------------------------------------------------------
!.----------------------------------------------------------------------


      if(mod(k,20000).EQ.0) then ! Inicio:para Pymol
          write(10,26)n      

              write(10,30)1,x(1),y(1),z(1),22,1+1
              write(10,30)2,x(2),y(2),z(2),22,2+1
              write(10,30)3,x(3),y(3),z(3),22,3+1
              write(10,30)4,x(4),y(4),z(4),22,4+1
              write(10,30)5,x(5),y(5),z(5),22,5+1
              write(10,30)6,x(6),y(6),z(6),22,6+1
              write(10,30)7,x(7),y(7),z(7),22,7+1
              write(10,30)8,x(8),y(8),z(8),22,8+1
              write(10,29)9,x(9),y(9),z(9),22,9+1

          do i=10,nuno
             if (i.EQ.n) then
              write(10,31)i,x(i),y(i),z(i),22,n
             else
              write(10,31)i,x(i),y(i),z(i),22,i+1
             endif
          enddo
       
       endif
       
       Enddo
      !cierra numero de interaciones para n=0.!#########################






      Do iTem=1,1
!########################## termo  #####################################
!- - - - - - Inicia el proceso de interacion - - - - - - - - - - - - - -
       Do k = 1, 500000


          CALL SOLVE(1, fx, fy, fz, ekin, dt, 123456, sig, epot, hami)


          Call Fuerzas(k,n,fx,fy,fz,epot)


          CALL SOLVE(2, fx, fy, fz, ekin, dt, 123456, sig, epot, hami)

!.------------------------------------------------------------------------------------
      if(mod(k,5000).EQ.0) then

         WRITE (43, *) k, SNGL(2*ekin/DBLE(3*N)),SNGL(Hami)

         Call RGyCM(k,n)      ! Radio de giro y centro de masa.

         write(3,*)k,SNGL(epot/(n)),SNGL(ekin/n),SNGL(epot/(n)+ekin/n)

      endif
!.------------------------------------------------------------------------------------

      if(mod(k,5000).EQ.0) then ! Inicio:para Pymol
          write(10,26)n

        !n1:
              write(10,30)1,x(1),y(1),z(1),22,1+1
              write(10,30)2,x(2),y(2),z(2),22,2+1
              write(10,30)3,x(3),y(3),z(3),22,3+1
              write(10,30)4,x(4),y(4),z(4),22,4+1
              write(10,30)5,x(5),y(5),z(5),22,5+1
        !n2:
              write(10,30)6,x(6),y(6),z(6),22,6+1
              write(10,30)7,x(7),y(7),z(7),22,7+1
              write(10,30)8,x(8),y(8),z(8),22,8+1
              write(10,29)9,x(9),y(9),z(9),22,9+1

          do i=10,nuno
             if (i.EQ.n) then
              write(10,31)i,x(i),y(i),z(i),22,n
             else
              write(10,31)i,x(i),y(i),z(i),22,i+1
             endif
          enddo

       endif

      End do

       !  if(mod(iTem,2).EQ.0) then !
       !     sig = 1.5d0
       !  else
            sig = 0.4d0
       !  endif
      end do
      !cierra numero de interaciones.!##################################
      
      




      
      
!########################## termo final  ###############################
      sig = 0.4d0
!- - - - - - Inicia el proceso de interacion - - - - - - - - - - - - - -
       Do k = 1, npasos


          CALL SOLVE(1, fx, fy, fz, ekin, dt, 123456, sig, epot, hami)


          Call Fuerzas(k,n,fx,fy,fz,epot)


          CALL SOLVE(2, fx, fy, fz, ekin, dt, 123456, sig, epot, hami)
      
! --------- ! imprime salida de coordenadas.. -------------
       if(mod(k,2000).EQ.0) then
        do i=1,n
         write(1,*)x(i),y(i),z(i)
        end do

      !  do i=1,nuno
      !   write(11,*)x(i),y(i),z(i)
      !  end do

      !  do i=nuno+1,ndos
      !   write(12,*)x(i),y(i),z(i)
      !  end do

      !  do i=ndos+1,ntre
      !   write(13,*)x(i),y(i),z(i)
      !  end do
       endif
       

!.-------calcula parametros---------------------------------------------
       If(k .GT. Tequili) then
!.----------------------------------------------------------------------
      if(mod(k,400).EQ.0) then

         WRITE (43, *) k, SNGL(2*ekin/DBLE(3*N)),SNGL(Hami)

         Call RGyCM(k,n)          ! Radio de giro y centro de masa.

         write(3,*)k,SNGL(epot/(n)),SNGL(ekin/n),SNGL(epot/(n)+ekin/n)

      endif
!.----------------------------------------------------------------------
      ENDIF
!.----------------------------------------------------------------------
      
      
      ! Inicio:para Pymol ----------------------------------------------
      if(mod(k,10000).EQ.0) then
          write(14,26)n

              write(14,30)1,x(1),y(1),z(1),22,1+1
              write(14,30)2,x(2),y(2),z(2),22,2+1
              write(14,30)3,x(3),y(3),z(3),22,3+1
              write(14,30)4,x(4),y(4),z(4),22,4+1
              write(14,30)5,x(5),y(5),z(5),22,5+1
              write(14,30)6,x(6),y(6),z(6),22,6+1
              write(14,30)7,x(7),y(7),z(7),22,7+1
              write(14,30)8,x(8),y(8),z(8),22,8+1
              write(14,29)9,x(9),y(9),z(9),22,9+1

          do i=10,nuno
             if (i.EQ.n) then
              write(14,31)i,x(i),y(i),z(i),22,n
             else
              write(14,31)i,x(i),y(i),z(i),22,i+1
             endif
          enddo
       endif
      ! Fin: para Pymol ------------------------------------------------
      
      
      End do
      !cierra termo final ##############################################
      !#################################################################
      
      
      

        Write(0,*) 'para. de hist.',cont5,((npasos-Tequili)/5000)
        write(0,*)
        write(0,*),'temperature fianl=',(2*ekin)/(3*n)
        write(0,*)
        write(0,*),'n,dt, npasos, Fw, kcon, sig Q'
        write(0,*) n, dt, npasos, Fw, kcon, sig, Q
        write(0,*)
        write(0,*) ,'AL, Altura_box,Rcero,Rcu,cut2'
        write(0,*) AL, Altura_box,Rcero,Rcu,cut2 



2      format(I10,4(4x,f20.16))
1      format(6(4x,F20.16))
14     format(I10,3(4x,F21.8),4(4x,F20.16))

      ! formato:pymol
30    format(4x,I3,2x,'O',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,3x,I3) !formato:n1
29    format(4x,I3,2x,'O',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,4x,I3) !formato:n1
31    format(4x,I3,2x,'O',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,4x,I3) !formato:n2

32    format(4x,I3,2x,'C',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,3x,I3)!formato:n3

35    format(4x,I3,2x,'O',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,3x,I3)!formato:n4
37    format(4x,I3,2x,'O',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,5x,I3)!formato:n4

34    format(4x,I3,2x,'C',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,5x,I3)!formato:n5


33    format(4x,I3,2x,'O',3x,f10.5,2x,f10.5,2x,f10.5,4x,I3,4x,I3)!formato:





26    format(2x,I4)
16    format(2x,I3)

          Close(0)
          Close(1)
          Close(3)
          Close(4)

          Close(43)
          Close(44)
          Close(45)
          Close(10)

          Close(11)
          Close(12)
          Close(13)
          Close(14)


!###        print*,'eso esooo esss ess totootot dodooo amigoos (enter)'
      end program Poly


!=======================================================================
!======================= Sub rutinas ===================================
!=======================================================================




!**********************************************************************
!        Se calcula el radio de giro y el centro de masa.
!**********************************************************************
       subroutine RGyCM(k,n)

        use coor

        implicit none

        integer (kind=4), INTENT(IN) :: n
        
        integer (kind=8), INTENT(IN) ::k
        
        integer (kind=4) :: i,j
        
        real (kind=8) ::rg(6),dpf(6)

        real (kind=8) ::srg
        
        real (kind=8) ::yy,zz

        integer (kind=8) :: nuno,ndos,ntre,ncua,ncin

        integer (kind=8) :: sumn1,sumn2,sumn3,sumn4

        INCLUDE 'valor_corte.inc'


            do j=1,6
             rg(j) = 0.d0
             dpf(j)= 0.d0
            enddo  


!Rg-Total******************************
          srg = 0.d0
          do i=1,n
            do j=1,n
                 yy=y(i)-y(j)
                 zz=z(i)-z(j)            
              srg = srg+(x(i)-x(j))**2.+(yy)**2+(zz)**2
            enddo
          enddo
              rg(1)=srg/(2.d0*(n**2.d0))
!Rg-1
      !    srg = 0.d0
      !    do i=1,nuno
      !      do j=1,nuno
      !           yy=y(i)-y(j)
      !           zz=z(i)-z(j)
      !        srg = srg+(x(i)-x(j))**2.+(yy)**2+(zz)**2
      !      enddo
      !    enddo
      !        rg(2)=srg/(2.d0*((nuno)**2.d0))




!*******************************************


!Distancia-Total              
                 yy=y(1)-y(n)
                 zz=z(1)-z(n)
              dpf(1) = (x(1)-x(n))**2.+(yy)**2+(zz)**2
              dpf(1) = SQRT(dpf(1))
!Distancia-1
      !           yy=y(1)-y(nuno)
      !           zz=z(1)-z(nuno)
      !        dpf(2) = (x(1)-x(nuno))**2.+(yy)**2+(zz)**2
      !        dpf(2) = SQRT(dpf(2))




       write(4,*)SNGL(rg(1)), SNGL(dpf(1))

   
4      format(I10,7(4x,F20.16))  !#####
       end subroutine
        


!***********************************************************************
!              Inicializacion de velocidades: momento cero
!***********************************************************************
      SUBROUTINE vel_init_poli(n,T,vm2)

       use veloc
               
      integer (kind=4), INTENT(IN) :: n
      Integer :: cont
      real (kind=8), INTENT(IN) :: T
      Real (kind=8) :: vmsum,vmsum2
      Real (kind=8) :: vmsumx,vmsumy,vmsumz
      Real (kind=8) :: vx0t,vy0t,vz0t
      Real (kind=8) :: vm2,vr,f


         vmsumx=0.0d0
         vmsumy=0.0d0
         vmsumz=0.0d0

         vm2 = 0.d0

          do cont=1,n
            vmsumx = vmsumx + vx(cont)
            vmsumy = vmsumy + vy(cont)
            vmsumz = vmsumz + vz(cont)
          enddo

         vmsumx = vmsumx/Real(n)
         vmsumy = vmsumy/Real(n)
         vmsumz = vmsumz/Real(n)

         vx0t = 0.d0
         vy0t = 0.d0
         vz0t = 0.d0

         vm2 = 0.d0

         do cont=1,n
            vx(cont) = (vx(cont)-vmsumx)
            vy(cont) = (vy(cont)-vmsumy)
            vz(cont) = (vz(cont)-vmsumz)

            vm2 = vm2+vx(cont)**2+vy(cont)**2+vz(cont)**2
         enddo
         
           f = sqrt(3.0*n*T/vm2)
           
           
         vm2 = 0.d0
         do cont=1,n
            vx(cont) = vx(cont)*f
            vy(cont) = vy(cont)*f
            vz(cont) = vz(cont)*f
            vm2 = vm2+vx(cont)**2+vy(cont)**2+vz(cont)**2

         end do

         do cont=1,n
            vx0t = vx0t + vx(cont)
            vy0t = vy0t + vy(cont)
            vz0t = vz0t + vz(cont)
         enddo

         vx0t = vx0t/Real(n)
         vy0t = vy0t/Real(n)
         vz0t = vz0t/Real(n)
         
      WRITE (0, 99001) vm2/(3.d0*n)
      WRITE (0,*)
      WRITE (0, 99002) Vx0t, Vy0t, Vz0t
      WRITE (0,*)

99001 FORMAT (' Initial temperature en sub. Velocidad : ', f5.3)
99002 FORMAT (' Velocity centre of mass : ', /, '          x = ', e8.2,
     &        /, '          y = ', e8.2, /, '          z = ', e8.2)

      END SUBROUTINE
      
      


!***********************************************************************
! ================ Suma de Fuerzas     =================================
!***********************************************************************

      Subroutine Fuerzas(k,n,fx,fy,fz,epot)

         use coor
         use parameteFuerza
         

         real (kind=8), INTENT(OUT) :: fx(n),fy(n),fz(n)
         
         real (kind=8), INTENT(OUT) :: epot
         
         integer (kind=4), INTENT(IN) :: n
         
         integer (kind=8), INTENT(IN) :: k
         
         real (kind=8) :: fxLJ(n),fyLJ(n),fzLJ(n)
         real (kind=8) :: nfxLJ(n),nfyLJ(n),nfzLJ(n)

         real (kind=8) :: fxf(n),fyf(n),fzf(n)
         real (kind=8) :: fbenx(n),fbeny(n),fbenz(n)
         real (kind=8) :: fxwall(n)
         
         real (kind=8) :: eneFle,eneFENE,eneLJ,enewall,neneLJ


         Do i = 1, n  !fuerzas en cero, inciales a cada paso.
            fx(i) = 0.d0
            fy(i) = 0.d0
            fz(i) = 0.d0

            fxwall(i)= 0.0

            fxf(i) = 0.d0
            fyf(i) = 0.d0
            fzf(i) = 0.d0

      !      fxLJ(i) = 0.d0
      !      fyLJ(i) = 0.d0
      !      fzLJ(i) = 0.d0

            nfxLJ(i) = 0.d0
            nfyLJ(i) = 0.d0
            nfzLJ(i) = 0.d0

            fbenx(i) = 0.d0
            fbeny(i) = 0.d0
            fbenz(i) = 0.d0

            
         end do

        !Fuerzas: que se toman en cuenta
           !Call fLJ(n,fxLJ,fyLJ,fzLJ,eneLJ)

           Call nfLJ(n,nfxLJ,nfyLJ,nfzLJ,neneLJ)

           Call fFENE(n,fxf,fyf,fzf,eneFENE)

           Call FBEND(n,fbenx,fbeny,fbenz,eneFle)


          do i = 1, n    !Suman fuerzas.
             fx(i) = fxf(i)+fxLJ(i)+fbenx(i)+nfxLJ(i)
             fy(i) = fyf(i)+fyLJ(i)+fbeny(i)+nfyLJ(i)
             fz(i) = fzf(i)+fzLJ(i)+fbenz(i)+nfzLJ(i)
          end do

           epot =  eneLJ + eneFENE + eneFle + neneLJ


      If(k .GT. 400) then ! guardar los datos de cada potencial.
      if(mod(k,iprint).EQ.0) then
      write(44,*)k,SNGL(eneFENE/n),SNGL(eneFle/n),
     & SNGL(eneLJ/n)+SNGL(eneFle/n)+SNGL(neneLJ/n)
      endif
      ENDIF


      END SUBROUTINE
      

!***********************************************************************
! ================ Fuerzas  de LJ  =====================================
!     Este es el potencial LJ de Model.
!***********************************************************************

      Subroutine fLJ(n,fxLJ,fyLJ,fzLJ,eneLJ)
      
        use coor
        use parameteFuerza
        
        
        implicit none
        
        real (kind=8), INTENT(OUT) :: fxLJ(n),fyLJ(n),fzLJ(n)

        real (kind=8), INTENT(OUT) :: eneLJ
          
        integer (kind=4), INTENT(IN) :: n

        Integer :: i,j
        real (kind=8) :: fxx,fyy,fzz
        real (kind=8) :: xx,yy,zz,r2,r1,r6,pot,rr
        
            eneLJ = 0.d0 ! monitoriar energia
            
           do i = 1, n-2   ! potencial de LJ. mono-mono - - - - 1 !
           do j = i+2, n
              xx = x(i)-x(j)
              yy = y(i)-y(j)
              zz = z(i)-z(j)

               r2 = xx**2+yy**2+zz**2

                  r1 = 1/r2
                  r6 = r1**3

                   pot=1.0d0*r6*(r6-1)
                   eneLJ = eneLJ + pot

                      rr =12.0d0*r6*r1*(r6-0.5d0)

                    fxx = rr*xx
                    fyy = rr*yy
                    fzz = rr*zz

                   fxLJ(i) = fxLJ(i)+fxx
                   fyLJ(i) = fyLJ(i)+fyy
                   fzLJ(i) = fzLJ(i)+fzz
                    fxLJ(j) = fxLJ(j)-fxx
                    fyLJ(j) = fyLJ(j)-fyy
                    fzLJ(j) = fzLJ(j)-fzz

           end do
           end do
           
      END SUBROUTINE
      
      
      
!***********************************************************************
! ================ Fuerzas  de FENE  ===================================
!     Potencial que forma el enlace: LJ(parete cortada)+ FENE
!***********************************************************************
      subroutine fFENE(n,fxf,fyf,fzf,eneFENE)

        use coor
        use parameteFuerza
        
        implicit none
        
        real (kind=8), INTENT(OUT) :: fxf(n),fyf(n),fzf(n)
        real (kind=8), INTENT(OUT) :: eneFENE

        integer (kind=4), INTENT(IN) :: n

        Integer :: i,j
        real (kind=8) :: fxx,fyy,fzz,nRcero
        real (kind=8) :: xx,yy,zz,r2,r1,r6,pot,rr
        
         eneFENE = 0.d0    ! monitoriar energia
         
         nRcero=1.d0/Rcero ! parametro para FENE: radio cero





!-------- potencial de LJ solo para enlace ----------------------------
           do i = 1, n-1   ! potencial de LJ. mono-mono - - - - 1 !
           j = i+1
              xx = x(i)-x(j)
              yy = y(i)-y(j)
              zz = z(i)-z(j)


               r2 = xx**2+yy**2+zz**2

                if (r2**0.5d0 .lt. Rcero) then

                  r1 = 1/r2
                  r6 = r1**3

                   pot=4.0d0*r6*(r6-1)

                   eneFENE = eneFENE + pot

                    rr = 48.d0*r6*r1*(r6-0.5d0)

                    fxx = rr*xx
                    fyy = rr*yy
                    fzz = rr*zz

                   fxf(i) = fxf(i)+fxx
                   fyf(i) = fyf(i)+fyy
                   fzf(i) = fzf(i)+fzz
                    fxf(j) = fxf(j)-fxx
                    fyf(j) = fyf(j)-fyy
                    fzf(j) = fzf(j)-fzz
                endif

                 if (r2**0.5d0 .lt. Rcero) then

         pot = -0.5d0*kfe*(Rcero**2.d0)*log(1.d0-(nRcero**2.d0)*r2)

                   eneFENE = eneFENE + pot


                   r1 = (-Kfe*Rcero**2.d0)/(Rcero**2.d0-r2)

                   fxx = r1*xx
                   fyy = r1*yy
                   fzz = r1*zz

                   fxf(i) = fxf(i)+fxx
                   fyf(i) = fyf(i)+fyy
                   fzf(i) = fzf(i)+fzz
                    fxf(j) = fxf(j)-fxx
                    fyf(j) = fyf(j)-fyy
                    fzf(j) = fzf(j)-fzz
                 end if
         end do




      end subroutine
      
      

      


!***********************************************************************
! ============= Solve: de ec. de mov  ==================================
!***********************************************************************
      SUBROUTINE SOLVE(Switch,Fx,Fy,Fz,Enkin,Delt,Iseed,Temp,Enpot,H)

!     using Velocity Verlet algorithm: Swope et al JCP 76 (1982) 637
!     switch = 1 : update position + first step velocity
!     switch = 2 : update second step velocity

        !use paramete
        use coor
        use veloc
        
        IMPLICIT NONE
        
       integer, parameter :: npmaxx=5000
        
      real (kind=8)  :: Fx(*), Fy(*), Fz(*), Enkin, Delt, delth, delt2,
     &                 Temp, vxn(npmaxx), vyn(npmaxx), vzn(npmaxx),
     &                 vxo(npmaxx), vyo(npmaxx), vzo(npmaxx), pso, psn,
     &                 err, sumv2, Enpot, H, ri, di, bx(npmaxx),
     &                 by(npmaxx), bz(npmaxx), delps
      INTEGER  i, Switch, Iseed, iter
      LOGICAL ready

      IF (Switch.EQ.1) THEN

!        --input: force and velocity at time t
!        --out: position time t+1 and part velocity

!        ===solve equations of motion

         delt2 = Delt*Delt/2
         delth = Delt/2
         sumv2 = 0

         DO i = 1, NPART
            X(i) = X(i) + Delt*VX(i) + delt2*(Fx(i)-PS*VX(i))
            Y(i) = Y(i) + Delt*VY(i) + delt2*(Fy(i)-PS*VY(i))
            Z(i) = Z(i) + Delt*VZ(i) + delt2*(Fz(i)-PS*VZ(i))
            sumv2 = sumv2 + VX(i)**2 + VY(i)**2 + VZ(i)**2
            VX(i) = VX(i) + delth*(Fx(i)-PS*VX(i))
            VY(i) = VY(i) + delth*(Fy(i)-PS*VY(i))
            VZ(i) = VZ(i) + delth*(Fz(i)-PS*VZ(i))
         END DO

!        S = S + PS*Delt + (sumv2-G*Temp)*delt2/Q !(esto se modifico para que Hami sea contante.)
         S = S +  (sumv2-G*Temp)*delt2/Q
         PS = PS + (sumv2-G*Temp)*delth/Q
         
      ELSE IF (Switch.EQ.2) THEN

!        ---final update velocity
         delth = Delt/2
         err = 1.D-10
         sumv2 = 0
         DO i = 1, NPART
            vxn(i) = VX(i)
            vyn(i) = VY(i)
            vzn(i) = VZ(i)
            sumv2 = sumv2 + vxn(i)*vxn(i) + vyn(i)*vyn(i) + vzn(i)
     &              *vzn(i)
         END DO
         psn = PS
         ready = .FALSE.
         iter = 0
         DO WHILE (.NOT.ready.AND.iter.LT.100)
            iter = iter + 1
            pso = psn
            delps = 0
            DO i = 1, NPART
               vxo(i) = vxn(i)
               vyo(i) = vyn(i)
               vzo(i) = vzn(i)
               bx(i) = -delth*(Fx(i)-pso*vxo(i)) - (VX(i)-vxo(i))
               ri = vxo(i)*Delt/Q
               delps = delps + ri*bx(i)
               by(i) = -delth*(Fy(i)-pso*vyo(i)) - (VY(i)-vyo(i))
               ri = vyo(i)*Delt/Q
               delps = delps + ri*by(i)
               bz(i) = -delth*(Fz(i)-pso*vzo(i)) - (VZ(i)-vzo(i))
               ri = vzo(i)*Delt/Q
               delps = delps + ri*bz(i)
            END DO
            di = -(pso*delth+1)
            delps = delps - di*((-sumv2+G*Temp)*delth/Q-(PS-pso))
            delps = delps/(-Delt*delth*sumv2/Q+di)
            sumv2 = 0
            DO i = 1, NPART
               vxn(i) = vxn(i) + (bx(i)+delth*vxo(i)*delps)/di
               vyn(i) = vyn(i) + (by(i)+delth*vyo(i)*delps)/di
               vzn(i) = vzn(i) + (bz(i)+delth*vzo(i)*delps)/di
               sumv2 = sumv2 + vxn(i)*vxn(i) + vyn(i)*vyn(i) + vzn(i)
     &                 *vzn(i)
            END DO
            psn = pso + delps

!           ---test for convergence
            ready = .TRUE.
            i = 0
            DO WHILE (i.LE.NPART.AND.ready)
               i = i + 1
               IF (i.LE.NPART) THEN
                  IF (ABS((vxn(i)-vxo(i))/vxn(i)).GT.err)
     &                ready = .FALSE.
                  IF (ABS((vyn(i)-vyo(i))/vyn(i)).GT.err)
     &                ready = .FALSE.
                  IF (ABS((vzn(i)-vzo(i))/vzn(i)).GT.err)
     &                ready = .FALSE.
               ELSE
                  IF (ABS((psn-pso)/psn).GT.err) ready = .FALSE.
               END IF
            END DO
         END DO
         DO i = 1, NPART
            VX(i) = vxn(i)
            VY(i) = vyn(i)
            VZ(i) = vzn(i)
         END DO

         PS = psn
         Enkin = sumv2/2.0
         H = Enkin + Enpot + (PS**2*Q)/2 + G*Temp*S

      ELSE
         STOP 'error switch'

      END IF
      
      
      end subroutine
      
      
      
!***********************************************************************
! ============= Calculo de fluctuaciones temp   ========================
!***********************************************************************
      
      SUBROUTINE SAMPLE(n,Switch,Is,En,Enk,H)

        use veloc

        IMPLICIT NONE
        
        integer, parameter :: NHIsmax=250

       integer :: Switch
       integer (kind=8) :: is
       integer (kind=8) :: i,itemp,ivel
       integer (kind=8) :: bintemp(NHIsmax),binvel(NHIsmax)
       integer (kind=8) :: binvel2(NHIsmax),n

       real (kind=8) :: v0,factemp,factvel,tempavv,tempav2,v2av,v4av,ill
       real (kind=8) ::  sum,temp,v2,H, En, Enk,factvel2


      IF(switch.Eq.1) then

      
        IF (N.NE.0) THEN

            temp = 2*Enk/DBLE(3*N)

            itemp = INT(temp*factemp)
            bintemp(itemp) = bintemp(itemp) + 1

           do i=1,n
             ivel = INT((VX(i)+v0)*factvel)
               binvel(ivel) = binvel(ivel) + 1
             ivel = INT((VY(i)+v0)*factvel)
               binvel(ivel) = binvel(ivel) + 1
             ivel = INT((VZ(i)+v0)*factvel)
               binvel(ivel) = binvel(ivel) + 1

               v2 = VX(i)**2 + VY(i)**2 + VZ(i)**2
               v2av = v2av + v2
               v4av = v4av + v2*v2
               ivel = INT(v2*factvel2) + 1
               binvel2(ivel) = binvel2(ivel) + 1
            END DO
            itt = itt + 1


            tempavv = tempavv + temp
            tempav2 = tempav2 + temp*temp
        ELSE

           stop  'cero particulas'

        END IF
            WRITE (43, *) Is, SNGL(temp),H!,SNGL(enp)



      Else if (switch.EQ.0) then
       ! valores inciales.
         itt = 0
         v0 = 10.D0
         factemp = NHIsmax/5.D0
         factvel = NHIsmax/(6*v0)
         factvel2 = NHIsmax/v0**2
         tempavv = 0
         tempav2 = 0
         v2av = 0
         v4av = 0
         ILL = 0
           

      Else if (switch.EQ.2) then
         sum =0
         IF(itt.NE.0) then
            write(0,*) '------ En Sub. SAMPLE s=2 ----'
            write(0,*) ' itt samples : ', itt
            write(0,*)
            
            sum =0
           do i=1, NHIsmax
             sum= sum + (binvel(i)/DBLE(N*itt*3))/factvel
           enddo
             sum = 1/sum

          do i=1, NHIsmax
             IF (bintemp(i).NE.0) THEN
                WRITE (0, '(2(2x,f8.3))') (i+0.5)/factemp,
     &                   (bintemp(i)/DBLE(itt))/factemp
             endif

             IF (binvel(i).NE.0) THEN
                WRITE (0, '(2(2x,f8.3))') (i+0.5)/factvel-v0,
     &                 sum*(binvel(i)/DBLE(N*itt*3))
             END IF

             IF (binvel2(i).NE.0) THEN
                WRITE (0, '(2(2x,f10.3))') (i-0.5)/factvel2,
     &                   (binvel2(i)/DBLE(N*itt))/factvel2
             END IF
          ENDDO

            tempavv = tempavv/itt
            write(0,*) ' temperature-2 ', SNGL(tempavv)
            write(0,*)
            
            
            tempav2 = tempav2/itt
            write(0,*) ' sigma(T) ',
     &            SNGL(SQRT((tempav2-tempavv**2)/tempavv**2)),
     &            SNGL(SQRT(N*(tempav2-tempavv**2)/tempavv**2))
            write(0,*)


            v2av = v2av/(itt*N)
            v4av = v4av/(itt*N)

            write(0,*) ' Velocity v^2 av', SNGL(v2av)
            write(0,*)
            write(0,*) ' sigma(v) ', SNGL(SQRT((v4av-v2av**2)/v2av**2))
            write(0,*)
            write(0,*) ' sum ', sum
            write(0,*) ' -----------------------'


         END IF

      Else
        Stop 'Error (sample.f) switch'

      EndIf
      

      
      end subroutine
      
      
! **********************************************************************
!***** Genera cordenadas inciales **************************************
!***********************************************************************

      subroutine inicio(n)

       use coor
       use veloc
       use parameteFuerza

       IMPLICIT NONE

       integer (kind=4), INTENT(IN) :: n
       Integer (kind=8) :: i,j,k
       real (kind=8) :: ram,Tfx,Tfy,Tfz,Ftotal,epot,ram1
       real (kind=8) :: fx(n),fy(n),fz(n)


        ! cordenadas: a mano.
          TFy = (1.15d0)*(0.130526d0) ! cos(a)
          TFz = (1.15d0)*(0.8091444d0) ! sin(a)

          x(1) = Altura_box/2.0d0+0.046
          y(1) = 0.d0
          z(1) = 0.d0

          do i=2,n,2         
                   Call random_number(ram)
            x(i) = Altura_box/2.0d0-0.05*ram
            y(i) = TFy
            z(i) = TFz*(i-1.d0)
          enddo      

          do i=3,n,2
                   Call random_number(ram)
            x(i) = Altura_box/2.0d0-0.051*ram
            y(i) = 0.0d0
            z(i) = TFz*(i-1.d0)
          enddo      





      
      ! se escala todo al cetro de masas
            Tfy = 0.d0
            Tfz = 0.d0
        do i=1, n
           TFy = Tfy + y(i)
           TFz = Tfz + z(i)
        enddo
           TFy = Tfy/n
           TFz = Tfz/n

        do i=1, n
           y(i) = y(i)- Tfy
           z(i) = z(i)- Tfz
        enddo



         k=1
         Call Fuerzas(k,n,fx,fy,fz,epot)

            Tfx = 0.d0
            Tfy = 0.d0
            Tfz = 0.d0

         do i=1, n
            Tfx = Tfx+fx(i)
            Tfy = Tfy+fy(i)
            Tfz = Tfz+fz(i)
            
         enddo
           

         Ftotal= Sqrt(Tfx**2.d0+Tfy**2.d0+Tfz**2.d0)

         do i=1, n
              Call random_number(ram)
            vx(i) = -ram*Tfx/Ftotal
              Call random_number(ram)
            vy(i) = -ram*Tfy/Ftotal
              Call random_number(ram)
            vz(i) = -ram*Tfz/Ftotal
         enddo

      end subroutine

!***********************************************************************
! ============= Fuerza para flex. interacion angular ===================
!     Calculo de fuerza flex, por rutina de articulo de Tamara y new,
!     con potencial de articulo de Model
!***********************************************************************

      SUBROUTINE FBEND(n,fbenx,fbeny,fbenz,sumcos)

        use coor
        use parameteFuerza

        implicit none

        real (kind=8), INTENT(OUT) :: fbenx(n),fbeny(n),fbenz(n)
        integer (kind=4), INTENT(IN) :: n
        Integer (kind=8) :: i
       ! real (kind=8) :: dfcos(n), sumcos
       real (kind=8) ::  cs(n),dcosx(3),dcosy(3),dcosz(3),csi,isen4,Pi
       real (kind=8) ::  sumcos,pot
       real (kind=8) ::  Kcon1,Kcon2,Kcon3,Kcon4,Kcon5
       Integer (kind=8) :: nuno,ndos,ntre,ncua,ncin


       INCLUDE 'valor_Kcon.inc'  !#####
       INCLUDE 'valor_corte.inc'

        sumcos =0.d0  ! monitorea energia.

        csi = -0.965925d0 !-0.642787d0

        do i=1,n
           fbenx(i)=0.0d0
           fbeny(i)=0.0d0
           fbenz(i)=0.0d0
        end do


! Segmento 1
        do i=2,nuno-1                                                !DERCOS(n,Ka,Le,Me,Ne,dcosx,dcosy,dcosz,coseno)

           Call dercos(n,i-1,i,i+1,i,dcosx,dcosy,dcosz,cs)        !cs es abreviacion de coseno


           pot = (1.d0 + cs(i))*kcon1 ! potencial de model

           sumcos = sumcos + pot


      fbenx(i-1)=fbenx(i-1)+(-1.d0*Kcon1)*dcosx(1)   !Todos los calculos deben usar cs(i). Ver dibujo.
                                                    !Se hace el calculo en una sola l/A1nea para evitar problemas al sumar los valores...
      fbeny(i-1)=fbeny(i-1)+(-1.d0*Kcon1)*dcosy(1)   !... de las fuerzas durante las iteraciones.

      fbenz(i-1)=fbenz(i-1)+(-1.d0*Kcon1)*dcosz(1)


      fbenx(i)=fbenx(i)+(-1.d0*Kcon1)*dcosx(2)

      fbeny(i)=fbeny(i)+(-1.d0*Kcon1)*dcosy(2)

      fbenz(i)=fbenz(i)+(-1.d0*Kcon1)*dcosz(2)


      fbenx(i+1)=fbenx(i+1)+(-1.d0*Kcon1)*dcosx(3)

      fbeny(i+1)=fbeny(i+1)+(-1.d0*Kcon1)*dcosy(3)

      fbenz(i+1)=fbenz(i+1)+(-1.d0*Kcon1)*dcosz(3)

        end do











       END SUBROUTINE

      
!***********************************************************************
! ============= Calculo de derivada del coseno de fi  ==================
!***********************************************************************

      SUBROUTINE DERCOS(n,Ka,Le,Me,Ne,dcosx,dcosy,dcosz,coseno)

        use veloc
        use coor

        IMPLICIT NONE

        integer, parameter :: NHIsmax=250

        integer (kind=4), INTENT(IN) :: n

       !real (kind=8) ::  num,den,den2,dnum(2),dden(2),dcosa(2)    !Algunas de las variables originales
       real (kind=8) ::  coseno(n),cs
       real (kind=8) ::  dnumx,dnumy,dnumz,dden1x,dden1y,dden1z
       real (kind=8) ::  dden2x,dden2y,dden2z,ddenx,ddeny,ddenz
       real (kind=8) ::  xx1,yy1,zz1,xx2,yy2,zz2
       real (kind=8) ::  num,den1,den2,den,dcosx(3),dcosy(3),dcosz(3)
       integer (kind=8) :: Ka,Le,Me,Ne,P,i,j

       double precision der1

       ! Me=Le    !Para mantener la notacion.

       do i=1,n
       coseno(i)=0.0d0
       end do

        xx1 = x(Ka)-x(Le)
        yy1 = y(Ka)-y(Le)
        zz1 = z(Ka)-z(Le)

        xx2 = x(Me)-x(Ne)    !Checa porque lo corren as/A1
        yy2 = y(Me)-y(Ne)
        zz2 = z(Me)-z(Ne)


       num=xx1*xx2+yy1*yy2+zz1*zz2  ! /A5
       den1=xx1*xx1+yy1*yy1+zz1*zz1
       den2=xx2*xx2+yy2*yy2+zz2*zz2
       den=(den1*den2)**0.5d0    ! ~D

       coseno(Le)=num/den
       cs=coseno(Le)

       j=1

       do P=Le-1,Le+1

       dnumx=0.0d0
       dnumy=0.0d0
       dnumz=0.0d0

       dden1x=0.0d0
       dden1y=0.0d0
       dden1z=0.0d0

       dden2x=0.0d0
       dden2y=0.0d0
       dden2z=0.0d0

       dnumx=DER1(Me,Ne,P)*xx1+DER1(Ka,Le,P)*xx2   !Derivada del numerador
       dnumy=DER1(Me,Ne,P)*yy1+DER1(Ka,Le,P)*yy2
       dnumz=DER1(Me,Ne,P)*zz1+DER1(Ka,Le,P)*zz2

      ! end do

       dden1x=2.0d0*(DER1(Ka,Le,P)*xx1)
       dden1y=2.0d0*(DER1(Ka,Le,P)*yy1)
       dden1z=2.0d0*(DER1(Ka,Le,P)*zz1)

       dden2x=2.0d0*(DER1(Me,Ne,P)*xx2)
       dden2y=2.0d0*(DER1(Me,Ne,P)*yy2)
       dden2z=2.0d0*(DER1(Me,Ne,P)*zz2)

       ddenx=(dden1x*den2+dden2x*den1)/(2.0d0*den)  !Derivada del denominador
       ddeny=(dden1y*den2+dden2y*den1)/(2.0d0*den)
       ddenz=(dden1z*den2+dden2z*den1)/(2.0d0*den)

       dcosx(j)=(dnumx-ddenx*cs)/den           !Derivada del Coseno del angulo
       dcosy(j)=(dnumy-ddeny*cs)/den
       dcosz(j)=(dnumz-ddenz*cs)/den

       j = j +1

       end do

      END SUBROUTINE
!***********************************************************************
! ============= Definicion: para calculo de deltas de dirac=============
!     Se ocupa para caluculo de fuerzas en Flex.
!***********************************************************************
       double precision function DER1(ii,jj,kk)

       integer (kind=8) :: ii,jj,kk

      der1=0.0d0
      if (kk .eq. ii) der1=1.0d0
      if (kk .eq. jj) der1=-1.0d0

      return
      end

!***********************************************************************
! ============= Definicion: LJ segmentos                   =============
!
!***********************************************************************


      Subroutine nfLJ(n,nfxLJ,nfyLJ,nfzLJ,neneLJ)
      
        use coor
        use parameteFuerza
        
        
        implicit none
        
        real (kind=8), INTENT(OUT) :: nfxLJ(n),nfyLJ(n),nfzLJ(n)

        real (kind=8), INTENT(OUT) :: neneLJ
          
        integer (kind=4), INTENT(IN) :: n

        Integer :: i,j
        real (kind=8) :: fxx,fyy,fzz
        real (kind=8) :: xx,yy,zz,r2,r1,r6,pot,rr,r5
        real (kind=8) :: pLJ1,pLJ2,pLJ3,pLJ4,pLJ5 !#####
        real (kind=8) :: fLJ1,fLJ2,fLJ3,fLJ4,fLJ5

        integer (kind=8) :: nuno, ndos, ntre, ncua, ncin

        INCLUDE 'valor_LJ.inc'
        INCLUDE 'valor_corte.inc'
        
            neneLJ = 0.d0 ! monitoriar energia

! Segmento 1
            
           do i = 1, nuno-2   ! potencial de LJ. mono-mono - - - - 1 !
           do j = i+2, nuno

              xx = x(i)-x(j)
              yy = y(i)-y(j)
              zz = z(i)-z(j)

               r2 = xx**2+yy**2+zz**2
                if (r2**0.5d0 .LT. 2.5d0) then
                  r1 = 1/r2
                  r5 = r1**5

                   pot=pLJ1*r5*(5.d0*r1-6.d0)+pLJ1*0.00054526d0
                   neneLJ = neneLJ + pot

                      rr = 60.d0*fLJ1*r5*r1*(r1-1.0d0)

                    fxx = rr*xx
                    fyy = rr*yy
                    fzz = rr*zz

                   nfxLJ(i) = nfxLJ(i)+fxx
                   nfyLJ(i) = nfyLJ(i)+fyy
                   nfzLJ(i) = nfzLJ(i)+fzz
                    nfxLJ(j) = nfxLJ(j)-fxx
                    nfyLJ(j) = nfyLJ(j)-fyy
                    nfzLJ(j) = nfzLJ(j)-fzz
               endif

           end do
           end do









           
      END SUBROUTINE
      

