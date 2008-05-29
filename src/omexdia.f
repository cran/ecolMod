!==========================================================================
!==========================================================================
! initialisation subroutine: 
! initialises the common block with 25 parameter values, 
! followed by thicknesses, porosities, bioturbation values
!==========================================================================
!==========================================================================

       subroutine initomexdia (steadyparms)
       external steadyparms
       integer,parameter :: N=100
       integer,parameter :: nc = 2*N + 3*(N+1) + 25

       double precision parms(nc)
       common /myparms/parms

       call steadyparms(nc, parms)
       
       return
       end


!==========================================================================
!==========================================================================
! subroutine calculating the rate of change
! the omexdia model
!==========================================================================
!==========================================================================

      subroutine omexdiamod (neq, t, Conc, dConc, yout, ip)
      
      implicit none


!......................... declaration section.............................
      integer           :: neq, ip, i
      integer,parameter :: N=100         
      double precision  :: t, Conc(*), dConc(*), yout(*), Flux(N+1)
      double precision :: por(N),intpor(N+1),Db(N+1),dx(N),dxInt(N+1)
      double precision :: Fdet(N),Sdet(N),O2(N),NO3(N),NH3(N),ODU(N)
      double precision :: dFdet(N),dSdet(N),dO2(N),dNO3(N),               &
     &                     dNH3(N),dODU(N)

      double precision :: cflux,Cprod(N),Nprod(N),Rescale(N)
      double precision :: Oxicminlim(N),Denitrilim(N),Anoxiclim(N) 
      double precision :: Oxicmin(N),Denitrific(N),anoxicmin(N)
      double precision :: nitri(N),oduox(N),odudepo(N)
      double precision :: pdepo

       
      double precision :: MeanFlux,rFast,rSlow ,pFast,w,NCrFdet,          & 
     &  NCrSdet,bwO2,bwNO3,bwNH3,bwODU,NH3Ads,rnit,ksO2nitri,             &
     &  rODUox,ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,       &
     &  kinO2anox,DispO2,DispNO3,DispNH3,DispODU   

      common /myparms    /MeanFlux,rFast,rSlow ,pFast,w,NCrFdet,          & 
     &  NCrSdet,bwO2,bwNO3,bwNH3,bwODU,NH3Ads,rnit,ksO2nitri,             &
     &  rODUox,ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,       &
     &  kinO2anox,DispO2,DispNO3,DispNH3,DispODU,dx,dxint,                &
     &  por,intpor,Db   

      double precision :: O2flux, O2deepflux, NO3flux, NO3deepflux
      double precision :: NH3flux,NH3deepflux,ODUflux,ODUdeepflux
      common /myout    /O2flux, O2deepflux, NO3flux, NO3deepflux,          &
     &                  NH3flux,NH3deepflux,ODUflux,ODUdeepflux

      character(len=80) msg
!............................ statements ..................................
!     check memory for output variables
!      if (ip(1) < 8)  call rexit("nout should be at least 8") 

! from Conc to fdet, sdet, o2,...
      CALL extract(N,6,Conc,Fdet,Sdet,O2,NO3,NH3,ODU)

      CFlux  = MeanFlux * (1+sin(2.*3.14158*t/365.))
       
! Rate of change due to transport: solid substances
      CALL tran1Dsol(FDET,cFlux*pFast,Db,w,por,intpor,dx,dxint,            &
     &               Flux,dFDET)

      CALL tran1Dsol(SDET,cFlux*(1.-pFast),Db,w,por,intpor,dx,dxint,       &
     &               Flux,dSDET)

! Rate of change due to transport: solute substances
      CALL tran1dliq(O2 ,bwO2 ,dispO2 ,w,por,intpor,dx,dxint,Flux,dO2)     
      O2flux = Flux(1)
      O2deepflux = Flux(N)
                        
      CALL tran1dliq(NO3,bwNO3,dispNO3,w,por,intpor,dx,dxint,Flux,dNO3)                        
      NO3flux = Flux(1)
      NO3deepflux = Flux(N)

      CALL tran1dliq(NH3,bwNH3,dispNH3/(1+NH3Ads),w,por,intpor,dx,dxint,   &
     &               Flux,dNH3)                        
      NH3flux = Flux(1)
      NH3deepflux = Flux(N)

      CALL tran1dliq(ODU,bwODU,dispODU,w,por,intpor,dx,dxint,Flux,dODU)                        
      ODUflux = Flux(1)
      ODUdeepflux = Flux(N)


! Production of DIC and DIN, expressed per cm3 LIQUID/day
      Cprod= (rFast*FDET        +rSlow*SDET        )*(1.-por)/por
      Nprod= (rFast*FDET*NCrFdet+rSlow*SDET*NCrSdet)*(1.-por)/por

! Oxic mineralisation, denitrification, anoxic mineralisation
! first the limitation terms
      Oxicminlim = O2/(O2+ksO2oxic)                ! limitation terms
      Denitrilim = (1.-O2/(O2+kinO2denit))*NO3/(NO3+ksNO3denit)
      Anoxiclim  = (1.-O2/(O2+kinO2anox))*(1-NO3/(NO3+kinNO3anox))
      Rescale    = 1./(Oxicminlim+Denitrilim+Anoxiclim)

! then the mineralisation rates
      OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
      Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
      AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
      Nitri      = rnit  *NH3*O2/(O2+ksO2nitri)
      OduOx      = rODUox*ODU*O2/(O2+ksO2oduox)

      pDepo      = min(1.d0,0.233*(w*365)**0.336 )
      OduDepo    = AnoxicMin*pDepo 

! Update the rate of change with biogeochemical processes
      dFDET = dFDET  - rFast*FDET
      dSDET = dSDET  - rSlow*SDET
      dO2   = dO2    - OxicMin          -2* Nitri -      OduOx
      dNH3  = dNH3   + (Nprod             - Nitri) / (1.+NH3Ads)
      dNO3  = dNO3   - 0.8*Denitrific     + Nitri 
      dODU  = dODU   + AnoxicMin                  - OduOx - OduDepo

! dfdet, dsdet, do2,... to dconc

      CALL reorder(N,6,dConc,dFdet,dSdet,dO2,dNO3,dNH3,dODU)
      CALL getout(yout)
      return
      end

!==========================================================================
! put output variables in one vector
!==========================================================================

      subroutine getout(yout)
      double precision :: yout(*), out(8)
      integer :: i

      common /myout    /out
      do i = 1, 8
       yout(i) = out (i)
      enddo       
 
      end
       
!==========================================================================
! put parameter values in a vector - NOT USED
!==========================================================================

      subroutine getparms(Np,is,yout)
       integer, parameter :: N=100
       integer,parameter ::nc = 2*N + 3*(N+1) + 25
       integer np,is,i,j
       character(len=80) msg
       
       double precision parms(nc),yout(*)
       common /myparms/parms
       j=1
       do i = is,np+is-1
         yout(j) = parms(i)
         j=j+1
       enddo
       end
      
!==========================================================================
! from ordered per species to ordered per slice
!==========================================================================

       subroutine extract(N,nspec,y,A,B,C,D,E,F)
       integer, intent(in) :: N, nspec
       integer :: i,j,ioff,ntot
       double precision, intent(in):: y(N*nspec)
       double precision, intent(out):: A(N),B(N),C(N),D(N),E(N),F(N)
       
       ntot = N*nspec
       ioff = 1
       j    = 1
       do i = ioff,ntot,nspec 
          A(j) = y(i)
          B(j) = y(i+1)          
          C(j) = y(i+2)          
          D(j) = y(i+3)                              
          E(j) = y(i+4)          
          F(j) = y(i+5)          
          j = j+1
       enddo
       end subroutine

!==========================================================================       
! from ordered per slice to ordered per species
!==========================================================================

       subroutine reorder(N,nspec,y,A,B,C,D,E,F)
       integer, intent(in) :: N, nspec
       integer :: i,j,ioff,ntot
       double precision, intent(out):: y(N*nspec)
       double precision, intent(in) :: A(N),B(N),C(N),D(N),E(N),F(N)
       
       ntot = N*nspec
       ioff = 1
       j    = 1
       do i = ioff,ntot,nspec 
          y(i)   = A(j) 
          y(i+1) = B(j) 
          y(i+2) = C(j) 
          y(i+3) = D(j) 
          y(i+4) = E(j) 
          y(i+5) = F(j) 
          j = j+1
       enddo
       end subroutine       
       
!==========================================================================
! transport of liquids
!==========================================================================

       subroutine tran1dliq(y,yup,disp,v,por,porint,dx,dxint,Flux,dy)

       integer,parameter :: N=100         
       double precision, intent(in):: y(N),disp,v,por(N),porint(N+1)
       double precision, intent(in):: dx(N),dxint(N+1),yup 
       double precision, intent(out):: flux(N+1),dy(N)
       integer :: i
       double precision :: ydown,porfac(N+1)

! boundaries
       ydown = y(N)       

! advection, corrected for porosity
         do i = 1,N+1
           porfac(i) = porint(N+1)/porint(i)
         ENDDO

! dispersion flux  
         do i =1,N-1
           Flux(i+1) = (y(i)-y(i+1))/dxint(i+1)
         enddo
         Flux(1)     = (yup-y(1))/dxint(1)  
         Flux(N+1)   = (y(N)-ydown)/dxint(N+1)
         do i =1,N+1
          Flux(i) = Disp*Flux(i)
         enddo 

! advection flux in direction of axis - centered differences
         do i =2,N+1
           Flux(i) = Flux(i) + 0.5*v*porfac(i)*y(i-1)
         enddo
           Flux(1) = Flux(1) + 0.5*v*porfac(1)*yup
         do i =1,N
           Flux(i) = Flux(i)     + 0.5*v*porfac(i)*y(i)
         enddo
           Flux(N+1) = Flux(N+1) + 0.5*v*porfac(N+1)*ydown
              
! rate of change = flux gradient, corrected for surface (porosity) changes 
       do i =1,N
         dy(i) = (Flux(i)*porint(i)-Flux(i+1)*porint(i+1))                &
     &             /por(i)/dx(i)  
       enddo
       
       end subroutine

!==========================================================================
! transport of solids
!==========================================================================

       subroutine tran1dsol(y,fluxup,disp,v,por,porint,dx,dxint,Flux,dy)

       integer,parameter :: N=100         
       double precision, intent(in):: y(N),fluxup,disp(N+1),v,por(N)
       double precision, intent(in):: dx(N),dxint(N+1),porint(N+1)
       double precision, intent(out):: flux(N+1),dy(N)
       integer :: i
       double precision :: yup,ydown,porfac(N+1)
       character(len=80) msg

! boundaries
       yup   = y(1)       
       ydown = y(N)       

! dispersion flux  
         do i =1,N-1
           Flux(i+1) = (y(i)-y(i+1))/dxint(i+1)
         enddo
         Flux(1)     = (yup - y(1) )/dxint(1)  
         Flux(N+1)   = (y(N)- ydown)/dxint(N+1)
         do i =1,N+1
          Flux(i) = Disp(i)*Flux(i)
         enddo 
! advection, corrected for porosity
         do i = 1,N+1
           porfac(i) = (1.-porint(N+1))/(1.-porint(i))
         ENDDO
! advection flux in direction of axis
         do i =2,N+1
           Flux(i) = Flux(i) + 0.5*v*porfac(i)*y(i-1)
         enddo
           Flux(1) = Flux(1) + 0.5*v*porfac(1)*yup
         do i =1,N
           Flux(i) = Flux(i)     + (1.-0.5)*v*porfac(i)*y(i)
         enddo
           Flux(N+1) = Flux(N+1) + (1.-0.5)*v*porfac(N+1)*ydown

! fluxes per bulk
       Flux = Flux * (1-porint)

! boundary fluxes
       Flux(1) = fluxup
              
! rate of change = flux gradient, corrected for surface (porosity) changes 
       do i =1,N
         dy(i) = (Flux(i)-Flux(i+1)) /(1.-por(i))/dx(i)         
       enddo
       
       end subroutine
       