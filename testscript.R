source("FunctionSimuls.R")
set.seed(123)

pop <- InitialisePopulation()
pop2 <- PreEventYears(pop = pop, start=1981, end=2100, survp = 0.8)
table(pop2$birthyear)
table(pop2$deathyear)

plot(pop2$birthyear, pop2$a)
plot(tapply(pop2$a, pop2$birthyear, mean))

summary(lm(a ~ birthyear, data=pop2))

