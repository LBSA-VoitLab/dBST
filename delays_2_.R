# Depays and DDEs # Figure 1



rm(list = ls())



if(!require('deSolve',character.only = TRUE)) {install.packages('deSolve');library(deSolve)}		# install and call deSolve package

t = seq(0,20,.1)		# set time for simulation

Pars = c(			# Parameters		
	tau = 2,
	Bolus  = 0		
      )

######################################################################################################################################
##### system of delayed dif. eq. X

InitStateX=c(		# initial state of the variables
		X1 = 1,
		X2 = 1	
		)

EquationsX <- function(t, y, pars, tau)  		# DDE model 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(y)), 
            {	
		tcount = t				# time sanity check
		tlag = t - tau			# time lag
		if (tlag <= 0) {X1_lag = 1} else {X1_lag = lagvalue(t=tlag,nr=1)}		# making sore that the lag does not get negative time values in the initial time points
		# VERY IMPORTANT : "lagvalue" is the function that retrives values from past times. If you are using more than one variable, you have to specify that variable are you retriving in "nr"
		# If  you want the lag value of the derivative, use "lagderiv". It works in a similar way to "lagvalue".

		eq1 = 1 - (X1_lag)^0.5
		eq2 = (X1_lag)^0.5 - X2

		dX1 = eq1
		dX2 = eq2

            return(list(dy = c(dX1,dX2),count=c(eq1,eq2),tcount=tcount,tlag=tlag,X1_lag=X1_lag))
            })
        }

perturbX = data.frame(var = 'X1', time = 3, value = 3, method = "add") #  method ("replace", "add", "multiply")		# perturbation done without rootfun

rootfun <- function(t, y, p) {return(t - 3)}; eventfun <- function(t, y, p) {return (c(y[1] +3, y[2]))}			# perturbation with rootfun triggering the events

outX = dede(			# run DEDE solver
		times = t,
		y = InitStateX,
		parms = Pars,
		func = EquationsX,
		tau = Pars[1],
		events = list(data = perturbX)) # 

		#rootfun = rootfun,
		#events = list(func = eventfun, root = TRUE)


par(mfrow=c(1,2))
plot(outX[,1],outX[,2],type='l',col='red',lwd=3,ylim=c(0,6),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
legend(10,5,legend=c('X1','Y1','Z1'),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','darkgreen','blue'),bty="n")
plot(outX[,1],outX[,3],type='l',col='red',lwd=3,ylim=c(0,2),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
legend(10,2,legend=c('X2','Y2','Z6'),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','darkgreen','blue'),bty="n")

##### system of delayed dif. eq. X
######################################################################################################################################



######################################################################################################################################
##### system of dif. eq. Y

InitStateY=c(		# initial state of the variables
		Y1 = 1,
		Y2 = 1 	
		)

EquationsY <- function(t, x, pars)  		# ODE model 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {	

		dY1 = 2 + Bolus - tau * Y1^0.5
		dY2 = tau * Y1^0.5 - tau * Y2

            return(list(c(dY1,dY2),count=c(),signal = c()))
            })
        }

perturbY = data.frame(var = 'Y1', time = 3, value = 3, method = "add") #  method ("replace", "add", "multiply")		# perturbation

outY = ode(			# run ODE solver
		times = t,
		y = InitStateY,
		parms = Pars,
		func = EquationsY,
		events = list(data = perturbY)) #  

par(mfrow=c(1,2))
plot(outY[,1],outY[,2],type='l',col='darkgreen',lwd=3,ylim=c(0,6),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
legend(10,5,legend=c('X1','Y1','Z1'),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','darkgreen','blue'),bty="n")
plot(outY[,1],outY[,3],type='l',col='darkgreen',lwd=3,ylim=c(0,2),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
legend(10,2,legend=c('X2','Y2','Z6'),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','darkgreen','blue'),bty="n")

##### system of dif. eq. Y
######################################################################################################################################



######################################################################################################################################
##### system of dif. eq. Z

InitStateZ=c(		# initial state of the variables
		Z1 = 1,
		Z2 = 1, 	
		Z3 = 1,
		Z4 = 1, 	
		Z5 = 1,
		Z6 = 1 	
		)

EquationsZ <- function(t, x, pars)   		# ODE model 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {	

		dZ1 = 2 + Bolus - tau * Z1^0.5
		dZ2 = tau * Z1^0.5 - tau * Z2
		dZ3 = tau * Z2 - tau * Z3
		dZ4 = tau * Z3 - tau * Z4
		dZ5 = tau * Z4 - tau * Z5
		dZ6 = tau * Z5 - 2 * Z6

            return(list(c(dZ1,dZ2,dZ3,dZ4,dZ5,dZ6),count=c(),signal = c()))
            })
        }

perturbZ = data.frame(var = 'Z1', time = 3, value = 3, method = "add") #  method ("replace", "add", "multiply")		# perturbation

outZ = ode(			# run ODE solver
		times = t,
		y = InitStateZ,
		parms = Pars,
		func = EquationsZ,
		events = list(data = perturbZ)) #  

par(mfrow=c(1,2))
plot(outZ[,1],outZ[,2],type='l',lty=2,col='blue',lwd=3,ylim=c(0,6),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
legend(10,5,legend=c('X1','Y1','Z1'),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','darkgreen','blue'),bty="n")
plot(outZ[,1],outZ[,7],type='l',lty=2,col='blue',lwd=3,ylim=c(0,2),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
legend(10,2,legend=c('X2','Y2','Z6'),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','darkgreen','blue'),bty="n")

##### system of dif. eq. Y
######################################################################################################################################



# tiff("Fig_delays1.tiff", height = 10, width = 20, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(1,2))

plot(outX[,1],outX[,2],type='l',col='red',lwd=3,ylim=c(0,5),xlab='t',ylab='',cex.lab=1,cex.axis=1)
lines(outY[,1],outY[,2],col='blue',lwd=3)
lines(outZ[,1],outZ[,2],lty=2,col='green',lwd=3)
legend(10,5,legend=c(expression('X'[1]),expression('Y'[1]),expression('Z'[1])),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','blue','green'),bty="n")
text(2,4.5,labels = "D",cex=2)

plot(outX[,1],outX[,3],type='l',col='red',lwd=3,ylim=c(0,5),xlab='t',ylab='',cex.lab=1,cex.axis=1)
lines(outY[,1],outY[,3],type='l',col='blue',lwd=3)
lines(outZ[,1],outZ[,7],lty=2,col='green',lwd=3)
legend(10,5,legend=c(expression('X'[2]),expression('Y'[2]),expression('Z'[6])),lty=c(1,1,2),lwd=c(3,3,3),col=c('red','blue','green'),bty="n")
text(2,4.5,labels = "E",cex=2)

# dev.off()
