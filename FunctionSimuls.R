InitialisePopulation <- function(nind=100, VA=1, VE=1, year=1980)
{
  pop <- data.frame(id=1:nind, sex=c(rep("F", times=round(nind/2)),
                              rep("M", times=nind - round(nind/2))),
                    birthyear = year,
                    deathyear = NA,
                    alive = 1,
                    dam = NA,
                    sire = NA,
                    a = rnorm(nind, 0, sd= sqrt(VA)),
                    e = rnorm(nind, 0, sd= sqrt(VE)))
  pop$z <- pop$a+pop$e
  
  return(pop)
}#end InitialisePopulation()

survfuction <- function(z=0, theta=0, omega=1) {
  exp( (-(theta - z)^2) / (2*(omega^2) ) )
}

RunPopulation <- function(start=1980, end=2000, maturity=1,
                          baserepro=2, reprovarf=1, reprovarm=1, K = 200, VE=1, VA=1,
                          maxsurvp = 0.8, theta=0, omega=1,
                          pop) 
{
  if(length(theta)==1)
  {
    theta <- rep(theta, times = end-start + 1)
  }else{if(length(theta)!= (end-start + 1))
    stop("Mismatch between theta and simulation duration")
  }
   
  for (year in start:end)
  {
    adultfemales <- which(pop$alive==1 & pop$sex=="F" & (year - pop$birthyear)>=maturity)
    adultmales <- which(pop$alive==1 & pop$sex=="M" & (year - pop$birthyear)>=maturity)
    
    offsp <- rpois(length(adultfemales),
                   lambda = exp(log(baserepro) + rnorm(length(adultfemales), 0, sqrt(reprovarf))) )
    # regulation by reproduction limitation, with random componant
    if(sum(offsp)+sum(pop$alive) > K)
    {
      p <- K/(sum(offsp)+sum(pop$alive))
      offsp <- sapply(offsp, function(x) rbinom(n = 1, size = x, prob = p))
    }
    
    if(any(is.na(offsp))){stop("here")}
    malereprofitness <- cumsum(exp(rnorm(length(adultmales), 0, sqrt(reprovarm))))
    
    if( sum(offsp)>0) #if any reproduction at all
    {
      newborn <- data.frame(id=(max(pop$id)+1):(max(pop$id)+sum(offsp)), 
                            sex=sample(c("F", "M"), size = sum(offsp), replace = TRUE),
                            birthyear = year,
                            deathyear = NA,
                            alive = 1,
                            dam = NA,
                            sire = NA,
                            a = 0,
                            e = rnorm(sum(offsp), 0, sd= sqrt(VE)))
      
      for(o in 1:sum(offsp))
      {
        currentoffsp <- min(which(cumsum(offsp) >= runif(1, min=0, max=sum(offsp))))
        newborn$dam[o] <- adultfemales[currentoffsp]
        offsp[currentoffsp] <- offsp[currentoffsp] - 1
        newborn$sire[o] <-adultmales[min(which(malereprofitness >= runif(1, min=0, max=max(malereprofitness))))]
        
        newborn$a[o] <- rnorm(n=1,
                              mean = (pop$a[pop$id==newborn$dam[o]] + pop$a[pop$id==newborn$sire[o]]) /2, 
                              sd = sqrt(VA/2))
        
      }
      newborn$z <- newborn$a + newborn$e
      pop <- rbind(pop, newborn)
    }#end if( sum(offsp)>0) 

    # viability selection
    survivors <- rbinom(n = sum(pop$alive), 
           size = 1,
           prob = maxsurvp*survfuction(pop$z[pop$alive==1],
                                       theta = theta[year - start + 1],
                                       omega = omega))
    
    sexofsurvivors <- pop$sex[pop$alive==1][which(survivors==1)]
    
    if((sum(sexofsurvivors=="F")*sum(sexofsurvivors=="M"))==0){
      warning(paste("Population went exctinct on year", year))
      pop$deathyear[pop$alive==1][which(survivors==0)] <- year
      pop$alive[pop$alive==1] <- survivors
      return(pop)
    }
    pop$deathyear[pop$alive==1][which(survivors==0)] <- year
    pop$alive[pop$alive==1] <- survivors
    
  }#end  for (year in start:end)
  
  return(pop)
}#end PreEventYears()


main_simul <- function(start=1980, end=2000, maturity=1,
                       baserepro=2, reprovarf=1, reprovarm=1, K = 200, VE=1, VA=1,
                       maxsurvp = 0.8, theta=0, omega=1){
  
  pop <- InitialisePopulation(nind = 100, VA = VA, VE = VE, year = start)
  pop <- RunPopulation(pop = pop, start=start, end=end, maturity = maturity, 
                       baserepro = baserepro, reprovarf = reprovarf, reprovarm = reprovarm,
                       K = K, VE=VE, VA=VA,
                       maxsurvp = maxsurvp, theta=theta, omega=omega)
  return(pop)
}
  
produce_pop_data <- function(pop)
{
  years <- vector()
  pop$deathyear[is.na(pop$deathyear)] <- max(pop$birthyear)
  for(i in 1:nrow(pop))
  {
    years <- c(years, pop$birthyear[i]:pop$deathyear[i])
  }
  return(as.data.frame(table(years)))
}
