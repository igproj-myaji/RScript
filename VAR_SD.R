VAR <- function(x){VAR<-var(x)*(length(x)-1)/length(x);VAR}
SD <- function(x){SD<-sqrt(VAR(x));SD}
