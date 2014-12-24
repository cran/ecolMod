
C***********************************************************
C     PERFORMS nstep COMPETITION STEPS IN LATTICE GRID
C***********************************************************

      SUBROUTINE lattice(nspec,ncell,nstep,cell,replace,sumdens,seed)
       IMPLICIT NONE
       INTEGER nspec,ncell,nstep,seed,ii
       INTEGER cell(ncell,ncell),newcell(ncell,ncell)
       INTEGER sumdens(nstep,nspec)
       DOUBLE PRECISION replace(nspec,nspec)

       INTEGER I,J,K,L
       INTEGER neighbour(4)
       DOUBLE PRECISION  replacement(4),rep
       DOUBLE PRECISION  rnd
       DOUBLE PRECISION normrnd 
       
C given species nr per cell and replacement probability, performs nstep steps and
C returns the updated species composition at each grid point for the final time step
C also returns the summed densities at each time step

C       call sRand(seed)    ! provokes an error on SUN - so it was removed
       CALL rndstart()
       DO I = 1,nstep
         DO J = 1, ncell
           DO K = 1, ncell

             rnd = normrnd()

             ii = cell(J,K)
             
             IF(J .GT. 1) THEN
               neighbour(1)= cell(J-1,K)
             ELSE
               neighbour(1)= cell(ncell,K)
             ENDIF
             
             IF(J .LT. ncell) THEN
               neighbour(2)= cell(J+1,K)
             ELSE
               neighbour(2)= cell(1,K)
             ENDIF
             
             IF(K .GT. 1) THEN
               neighbour(3)= cell(J,K-1)
             ELSE
               neighbour(3)= cell(J,ncell)
             ENDIF
             
             IF(K .LT. ncell) THEN
               neighbour(4)= cell(J,K+1)
             ELSE
               neighbour(4)= cell(J,1)
             ENDIF
             
             Rep          = 0
             newcell(J,K) = ii

C cumulative probabilities
             DO L = 1, 4
               replacement(L) = rep+replace(neighbour(L),ii)/4

               IF (rnd .GE. rep .AND. rnd.LT.replacement(L)) THEN
                 newcell(J,K) = neighbour(L)
                 exit
               ENDIF
               rep = replacement(L)

             ENDDO
           ENDDO 
         ENDDO

         DO J = 1, nspec
           sumdens(I,J) = 0
         ENDDO
         DO J = 1, ncell
           DO K = 1, ncell
             ii        = newcell(J,K) 
             cell(J,K) = ii 
             sumdens(i,ii) = sumdens(i,ii)+1
           ENDDO
         ENDDO

       ENDDO
       CALL rndend()       
       END SUBROUTINE
       
       
       
