source("FunctionSimuls.R")
set.seed(123)

pop <- InitialisePopulation()
pop <- RunPopulation(pop = pop, start=1981, end=2000)


pop$birthyear:pop$deathyear



plot(table(pop$birthyear))
table(pop$deathyear)

plot(pop$birthyear, pop$a)
plot(tapply(pop$a, pop$birthyear, mean))

summary(lm(a ~ birthyear, data=pop))


survfuction <- function(z=0, theta=0, omega=1) {
  exp( (-(theta - z)^2) / (2*(omega^2) ) )
}

curve(expr = survfuction(z=x), from = -5, to = 5, n = 100)
