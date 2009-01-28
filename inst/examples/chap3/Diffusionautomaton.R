###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 3.                                ##
## Spatial components and transport          ##
## Chapter 3.6.3.                            ##
## Cellular automaton model of diffusion.    ##
###############################################

#=====================================================

Diffuse <- function (Particles)
#-----------------------------------------------------
# Performs one diffusion step
#-----------------------------------------------------

  {
  
  # select the positions where there are Particles 
  x      <- which(Particles>0)
  
  # if random number > 0.5, particle moves to right (+1) if < 0.5 to left (-1)
  rnd            <- runif(length(x))   
  move           <- c(x[rnd<0.5]-1,x[rnd>=0.5]+1)

  # wrap around edges if necessary 
  move[move<1]     <- move[move<1]    +ncell
  move[move>ncell] <- move[move>ncell]-ncell  

  # Mark the new positions that are not empty              
  notfree <- which (Particles[move]>0)

  # movement allowed only when there is just one element moves to new position             
  duplo  <- which(move%in%move[duplicated(move)])
  free   <- move[-c(duplo,notfree)]
  source <- c(x[rnd<0.5],x[rnd>=0.5])[-c(duplo,notfree)]

  # source particles are emptied, new positions become 1
  Particles[source] <- 0
  Particles[free]   <- 1  
  return(Particles) 

} # END SUBROUTINE Diffuse

#=====================================================



##############################
# Initialising               #
##############################
  
ncell      <- 200
Particles  <- rep (0,ncell) 
Particles [c(10:30,50:70,90:110,130:150,170:190)] <- 1
nsteps     <- 1000

##############################
# RUNNING the model:         #
##############################

# perform nsteps diffusion steps 

Grid <- matrix(ncol=nsteps, nrow=ncell)

for (j in 1:nsteps) {
  Particles <- Diffuse(Particles)
  Grid[,j]  <- Particles
         }   

##############################
# PLOTTING model output:     #
##############################

windows(width=4,height=8)   ## A4
par(oma = c(0,0,3,0))
image (y=1:nsteps,z=Grid, col = c(0,1), axes=FALSE, 
       ylim=c(nsteps,1),xlab="",ylab="")
box()
mtext(outer=TRUE,side=3,"Diffusion in cellular automaton",cex=1.5)
   