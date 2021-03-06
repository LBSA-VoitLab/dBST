# Limit cycles # Figure 4



rm(list = ls())



######################################################################################################################################
##### time, initState and truePars 

t = seq(0,10000,.1)		# set time for simulation

state21 = c(		# initial state that produces damped oscilations
		X1 = 1,		
		X2 = 1.3	
		) 

state22 = c(		# initial state that produces increasing oscilations
		X1 = 1.1,		
		X2 = 1.8	
		) 

state23 = c(		# stable orbit
		X1 = 0.7249193,
		X2 = 0.8685822
		) 

Pars2 = c(			# Parameters	
		a1 = 0.011,
		a2 = 0.01
      	)

# equations
equa1 = expression(a1 * (X1^2 * X2^3 - X1 * X2))
equa2 = expression(a2 * (X1^(-1) * X2^4 - X1^3 * X2^5))


##### time, initState and truePars
######################################################################################################################################



######################################################################################################################################
##### system of dif. eq.

if(!require('deSolve',character.only = TRUE)) {install.packages('deSolve');library(deSolve)}

Equations2 <- function(t, x, pars)  		# ODE model
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {
	
		eq1 = eval(equa1)
		eq2 = eval(equa2)

		dX1 = eq1
		dX2 = eq2

            return(list(c(dX1,dX2),count=c(eq1,eq2)))
            })
        }

##### system of dif. eq.
######################################################################################################################################



######################################################################################################################################
##### ode 

ode_output21 = ode_output2 = ode(times = t, y = state21, parms = Pars2, func = Equations2)			# run ODE solver
ode_output22 = ode_output2 = ode(times = t, y = state22, parms = Pars2, func = Equations2)			# run ODE solver
ode_output23 = ode_output2 = ode(times = t, y = state23, parms = Pars2, func = Equations2)			# run ODE solver
# tail(out) # edit(edit) # View(out)

par(mfcol=c(3,2))
for (i in 1:2)
	{
	plot(ode_output21[,1],ode_output21[,i+1],pch=20,col='grey',lwd=3,xlab='t',ylab=colnames(ode_output2)[i+1],xlim=c(0,1000))
	plot(ode_output22[,1],ode_output22[,i+1],pch=20,col='grey',lwd=3,xlab='t',ylab=colnames(ode_output2)[i+1],xlim=c(0,1000))
	plot(ode_output23[,1],ode_output23[,i+1],pch=20,col='grey',lwd=3,xlab='t',ylab=colnames(ode_output2)[i+1],xlim=c(0,1000))
	}
# legend(60,1.0,legend = c('true model'),pch = c(20),col = c('grey'),bty = "n")	# ,lty = c(0,1)

##### ode
######################################################################################################################################



######################################################################################################################################
##### discrete 

DiscreteSolver = function (stateD,ode_out,letter=NULL)		# discrete solver
	{ 
	d_output2 = c(t=min(t),stateD)					# initiating variable to store the solution of the discrete solver
	for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )		# unnecessarily extremely complicated way to cycle through the values of time
		{
		if ( is.null(dim(d_output2))==TRUE )			# if d_output is a vector, i.e., if i=1
			{
			X1 = d_output2[2]
			X2 = d_output2[3]
			} else 							# if it is a matrix
			{
			X1 = d_output2[dim(d_output2)[1],2]
			X2 = d_output2[dim(d_output2)[1],3]
			}
	
		Xs = with(as.list(c(c(X1,X2), Pars2)), 
	            {
			
			X1_temp = eval(equa1) * (t[2]-t[1]) + X1
			X2_temp = eval(equa2) * (t[2]-t[1]) + X2
	
	            return(list(c(X1_temp,X2_temp)))
	            })
	
		d_output2 = rbind(d_output2,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2]))		# discrete result
		}

	# plots
	plot(d_output2[,1],d_output2[,2],pch=20,col = "blue",xlab='t',ylim=c(0,2.5),ylab = '',cex.lab=1.5,cex.axis=1.5,xlim=c(0,1000))
	#lines(ode_out[,1],ode_out[,2],col='grey')
	points(d_output2[,1],d_output2[,3],pch=20,col = "darkcyan")
	#lines(ode_out[,1],ode_out[,3],col='grey')
	legend(600,.5,legend=c(expression('X'[1]),expression('X'[2])),pch=c(20,20),pt.cex = c(2,2),lty=c(NA,NA),lwd=c(NA,NA),col=c('blue','darkcyan'),bty="n",cex=.9)
	text(10,2.3,labels = letter,cex=2)

	return(d_output2)
	}

# tiff("Fig_limit_cycles.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)

par(mfrow=c(2,2))
d_output21 = DiscreteSolver(stateD = state21, ode_out = ode_output21,letter="A")
d_output22 = DiscreteSolver(stateD = state22, ode_out = ode_output22,letter="B")
d_output23 = DiscreteSolver(stateD = state23, ode_out = ode_output23,letter="C")

plot(d_output21[,2],d_output21[,3],type='l',col='darkgreen',lwd=2,xlim=c(0.6,1.4),ylim=c(0,2.5),xlab=expression('X'[1]),ylab=expression('X'[2]),cex.lab=1.2,cex.axis=1.5)
points(d_output21[1,2],d_output21[1,3])
lines(d_output22[,2],d_output22[,3],col='green',lwd=2)
points(d_output22[1,2],d_output22[1,3])
lines(d_output23[,2],d_output23[,3],col='cyan',lwd=2)
points(d_output23[1,2],d_output23[1,3])
legend(.6,.7,
		legend=c('Initial Cond.',expression(paste('X'[1])~' = 1, '~'X'[2]~' = 1.3'),
						expression(paste('X'[1])~' = 1.1, '~'X'[2]~' = 1.8'),	
						expression(paste('X'[1])~' = 0.7249193, '~'X'[2]~' = 0.8685822')),	
				lty = c(NA,1,1,1),lwd= c(NA,2,2,2),col=c(NA,'darkgreen','green','cyan'),bty="n",cex=.8)
text(0.65,2.3,labels = "D",cex=2)
# dev.off()

##### discrete 
######################################################################################################################################
legend(.6,.7,legend=c('Initial Cond.',expression(paste('X'[1])~' = 1')))


# ODEs phase plane
plot(ode_output21[,2],ode_output21[,3],type='l',col='darkgreen',lwd=2,xlim=c(0.6,1.4),ylim=c(0,2.5),xlab='x1',ylab='x2',cex.lab=1.5,cex.axis=1.5)
points(ode_output21[1,2],ode_output21[1,3])
lines(ode_output22[,2],ode_output22[,3],col='green',lwd=2)
points(ode_output22[1,2],ode_output22[1,3])
lines(ode_output23[,2],ode_output23[,3],col='cyan',lwd=2)
points(ode_output23[1,2],ode_output23[1,3])
legend(.6,.7,legend=c('Initial Cond.','x1 = 1, x2 = 1.3','x1 = 1.1, x2 = 1.8','x1 = 0.7249193, x2 = 0.8685822'),lty = c(NA,1,1,1),lwd= c(NA,2,2,2),col=c(NA,'darkgreen','green','cyan'),bty="n",cex=.8)





