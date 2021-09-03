####function for finding initial conditions on variables that don't 
#make it to equilibrium before 200 yrs

#usage: i0=current initial condition, i1=condition at last time step, returns the 
#value that should be inserted into the initial conditions spot
findInitials<-function(i0, i1){
  inew<-ifelse((i1 < i0+i0*0.001) && (i1 > i0-i0*0.001), as.numeric(i0), as.numeric(i1))
  return(inew)}

