# Discrete Version of Kumaraswamy Log-Logistic distribution

# Discrete Kumaraswamy log-logistic probability distribution
fklld <- function(t,gamma,alpha,a,b){ 
  prob <- (1-(((t)^gamma)/((t)^gamma + alpha^gamma))^a)^b - (1-(((t+1)^gamma)/((t+1)^gamma + alpha^gamma))^a)^b
  return(prob)
}
# Discrete Kumaraswamy log-logistic survival function 
sklld <- function(t,gamma,alpha,a,b){ 
  kum <- 1 - (1 - (((t+1)^gamma)/((t+1)^gamma + alpha^gamma))^a)^b
  sob <- 1- kum
  return(sob)
}

# Discrete Kumaraswamy log-logistic Log-Likelihood function
lvero_klld <- function(par){
  p1 <- sum(status*log(fklld(t,par[1],par[2],par[3],par[4])) + (1-status)*log(sklld(t,par[1],par[2],par[3],par[4])))
  return(p1)
}



# Example with kidney data

library(survival)
data("kidney") #data form the survival package

status <- kidney$status #getting the status to fit the function
t <- kidney$time #getting the time to fit the function

t <- findInterval(t, seq(0,450,7)) # Putting in 7 days intervals as this is a Discrete function


# Optimization to find the estimators

# c(.01,.01,2,.1)
# c(1,1,1,1)
ML_est <- optim(c(.01,.01,2,.1) , lvero_klld,
                 control=list(fnscale=-1,maxit=1000), hessian=T)
# fnscale: . If negative, turns the problem into a maximization problem
#convergence: An integer code. 0 indicates successful completion
ML_est$convergence==0

invR <- solve((-1)*ML_est$hessian) #Hessian Matrix
any(diag(invR)<0)# Check if any Var is negative
var_betas <- diag(invR)[1:4]
std <- sqrt(var_betas)  # getting the standard deviation


gamma <- ML_est$par[1]
alpha <- ML_est$par[2]
a <- ML_est$par[3]
b <- ML_est$par[4]


km0 <- survfit(Surv(t,status)~1,conf.int=T)


kumaj <- data.frame(tempo=km0$time,sob=km0$surv,kuma=sklld(km0$time,gamma,alpha,a,b))


library(ggplot2)
ajuste <- ggplot(kumaj) + geom_step(aes(tempo,sob,
                                        color ="Kaplan-Meier",linetype="Kaplan-Meier")) + 
  geom_step(aes(tempo,kuma,color ="Discrete K-LogLog",linetype="Discrete K-LogLog")) +
  labs(y="S(t)", x="Time") +
  scale_colour_manual("", breaks = c("Kaplan-Meier","Discrete K-LogLog"),
                      values = c("deeppink3",1))+
   scale_linetype_manual("",breaks = c("Kaplan-Meier","Discrete K-LogLog"),
                        values = c("Kaplan-Meier"=1,"Discrete K-LogLog"=2)) +
  scale_y_continuous(limits=c(0,1)) + theme_light() +
  theme(legend.position =c(0.75, 0.75),legend.text=element_text(size=rel(.9)),legend.title = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"))  

ajuste



#------------------------- Regression Model

x0 <- rep(1,length(t))
x1 <- kidney$frail

log_reg <- function(par){
  #beta0, beta1, gamma are the parameters of parametrization
  beta0 <- par[1]
  beta1 <- par[2]
  gamma <- par[3]
  a <- par[4]
  b <- par[5]
  alpha <- exp(beta0*x0+beta1*x1)
  p1 <- sum(status*log(fklld(t,gamma,alpha,a,b)) + (1-status)*log(sklld(t,gamma,alpha,a,b)))
  return(p1)
}
ML_est_reg <- optim(c(1,.1,.1,1,1), log_reg,
                    control=list(fnscale=-1,maxit=10000),hessian=T)

invR <- solve((-1)*ML_est_reg$hessian)
any(diag(invR==0))
var_betas <- diag(invR)[1:2]
std <- sqrt(var_betas)
betas <- ML_est_reg$par[1:2]
tbetas <- betas/std
pvalue <- round(2*(1-pnorm(abs(tbetas))),6)
Est_betas <- cbind(betas,std, pvalue)
rownames(Est_betas) <- paste("beta",c(0,1),sep="")
Est_betas


