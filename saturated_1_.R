# adatpation to perturbation # Figure 3 and 6



rm(list = ls())



######################################################################################################################################
##### system of dif. eq.

if(!require('deSolve',character.only = TRUE)) {install.packages('deSolve');library(deSolve)}		# install and call deSolve package

t = seq(0,10,.5)		# set time for simulation

InitState=c(		# initial state of the variables
		x1 = 1.181653,
		x2 = 0.64
		)

Pars = c(			# Parameters
	gamma.0_1 = 2, x0 = 1, f.0_1.0 = 1,			
	r = 1.75, f.1_2.1 = .8,
	gamma.2_ = 2.5, f.2_.2 = .5				
      )

x0 = function(t)
	{	
	# Create a function to alter independent variables #
	if (t<5) {return(1)}
	if (t>=5) {return(1.2)}	
	}
x0(t=3)
x0(t=5)
x0(t=7)

Equations <- function(t, x, pars)  		# ODE model 
        { 
        ### returns rate of change
        # t = model's time structure
        # x = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {
		
		X0 = x0(t)		# evaluate x0 against time

		flux.0_1 = gamma.0_1 * X0^f.0_1.0
		flux.1_2 = r * x1^f.1_2.1
		flux.2_ = gamma.2_ * x2^f.2_.2		

		dx1 = flux.0_1 - flux.1_2
		dx2 = flux.1_2 - flux.2_

            return(list(c(dx1,dx2),count=c(flux.0_1, flux.1_2, flux.2_),signal = X0))
            })
        }


##### system of dif. eq.
######################################################################################################################################

# perturb = data.frame(var = 'x1', time = 6, value = 2, method = "replace") #  method ("replace", "add", "multiply")

out1 = ode(			# run ODE solver
		times = t,
		y = InitState,
		parms = Pars,
		func = Equations) # events = list(data = perturb) 

#windows()
plot(out1[,1],out1[,2],type='l',col='blue',lwd=3,ylim=c(0,2),xlab='t',ylab='')
points(out1[,1],out1[,3],type='l',col='lightblue',lwd=3)
legend(15,.3,legend=c('x1','x2'),lty=c(1,1),lwd=c(3,3),col=c('blue','lightblue'),bty="n")





######################################################################################################################################
##### discrete 

d_output1 = c(t=min(t),InitState)					# initiating variable to store the solution of the discrete solver
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )		# unnecessarily extremely complicated way to cycle through the values of time
	{
	if ( is.null(dim(d_output1))==TRUE )			# if d_output is a vector, i.e., if i=1
		{
		x1 = d_output1[2]
		x2 = d_output1[3]
		} else 							# if it is a matrix 
		{
		x1 = d_output1[dim(d_output1)[1],2]
		x2 = d_output1[dim(d_output1)[1],3]
		}

	Xs = with(as.list(c(c(x1,x2), Pars)), 
            {
		
		if (i>5) {x0 = 1.2}

		x1_temp = ( gamma.0_1 * x0^f.0_1.0 - r * x1^f.1_2.1 ) * (t[2]-t[1]) + x1
		x2_temp = ( r * x1^f.1_2.1 - gamma.2_ * x2^f.2_.2 ) * (t[2]-t[1]) + x2

            return(list(c(x1_temp,x2_temp)))
            })

	d_output1 = rbind(d_output1,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2]))		# discrete result
	}

# tiff("Fig_saturated.tiff", height = 15, width = 15, units = 'cm',compression = "lzw", res = 300)
plot(d_output1[,1],d_output1[,2],pch=20,col='blue',cex = 2,xlab='t',ylab=colnames(out1)[i+1],ylim=c(0,2),cex.lab=2,cex.axis=2)
lines(out1[,1],out1[,2],col = "darkgrey",lwd=2)
points(d_output1[,1],d_output1[,3],pch=20,cex = 2,col='darkcyan',xlab='t',ylab=colnames(out1)[i+1])
lines(out1[,1],out1[,3],col = "darkgrey",lwd=2)
legend(7,.5,legend=c(expression('X'[1]),expression('X'[2]),'ODEs'),pch=c(20,20,NA),pt.cex = c(2,2,NA),lty=c(NA,NA,1),lwd=c(NA,NA,2),col=c('blue','darkcyan','darkgrey'),cex=1.5,bty="n")
# dev.off()

##### discrete 
######################################################################################################################################





######################################################################################################################################
##### discrete with noise

t = seq(0,20,.5)			# set time for simulation

out1 = ode(				# run ODE solver
		times = t,
		y = InitState,
		parms = Pars,
		func = Equations) # events = list(data = perturb) 

d_output2 = c(t=min(t),InitState)					# initiating variable to store the solution of the discrete solver
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )		# unnecessarily extremely complicated way to cycle through the values of time
	{
	if ( is.null(dim(d_output2))==TRUE )			# if d_output is a vector, i.e., if i=1
		{
		x1 = d_output2[2]
		x2 = d_output2[3]
		} else 							# if it is a matrix
 
		{
		x1 = d_output2[dim(d_output2)[1],2]
		x2 = d_output2[dim(d_output2)[1],3]
		}

	Xs = with(as.list(c(c(x1,x2), Pars)), 
            {
		# discrete equations
		if (i>5) {x0 = 1.2}				# perturbation at time t=5

		r_rand = rnorm(1,mean = r,sd = .1)		# random variable

		x1_temp = ( gamma.0_1 * x0^f.0_1.0 - r_rand * x1^f.1_2.1 ) * (t[2]-t[1]) + x1
		x2_temp = ( r_rand * x1^f.1_2.1 - gamma.2_ * x2^f.2_.2 ) * (t[2]-t[1]) + x2

            return(list(c(x1_temp,x2_temp)))
            })

	d_output2 = rbind(d_output2,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2]))		# discrete result
	}



d_output3 = c(t=min(t),InitState)					# initiating variable to store the solution of the discrete solver
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )		# unnecessarily extremely complicated way to cycle through the values of time
	{
	if ( is.null(dim(d_output3))==TRUE )			# if d_output is a vector, i.e., if i=1
		{
		x1 = d_output3[2]
		x2 = d_output3[3]
		} else 							# if it is a matrix 
		{
		x1 = d_output3[dim(d_output3)[1],2]
		x2 = d_output3[dim(d_output3)[1],3]
		}

	Xs = with(as.list(c(c(x1,x2), Pars)), 
            {
		# discrete equations		
		if (i>5) {x0 = 1.2}					# perturbation at time t=5

		r_rand = rnorm(1,mean = r,sd = .1)			# random variable

		x1_temp = ( gamma.0_1 * x0^f.0_1.0 - r_rand * x1^f.1_2.1 ) * (t[2]-t[1]) + x1
		x2_temp = ( r_rand * x1^f.1_2.1 - gamma.2_ * x2^f.2_.2 ) * (t[2]-t[1]) + x2

            return(list(c(x1_temp,x2_temp)))
            })

	d_output3 = rbind(d_output3,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2]))		# discrete result
	}

# tiff("Fig_saturated_rand.tiff", height = 15, width = 30, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,2))
plot(d_output2[,1],d_output2[,2],pch=20,col='blue',cex = 2,xlab='t',ylab=colnames(out1)[i+1],ylim=c(0,2),cex.lab=2,cex.axis=2)
lines(out1[,1],out1[,2],col = "darkgrey",lwd=2)
points(d_output2[,1],d_output2[,3],pch=20,cex = 2,col='darkcyan',xlab='t',ylab=colnames(out1)[i+1])
lines(out1[,1],out1[,3],col = "darkgrey",lwd=2)
legend(13,.5,legend=c(expression('X'[1]),expression('X'[2]),'ODEs'),pch=c(20,20,NA),pt.cex = c(2,2,NA),lty=c(NA,NA,1),lwd=c(NA,NA,2),col=c('blue','darkcyan','darkgrey'),cex=1.5,bty="n")
text(0,2,labels = "A",cex=2)

plot(d_output3[,1],d_output3[,2],pch=20,col='blue',cex = 2,xlab='t',ylab=colnames(out1)[i+1],ylim=c(0,2),cex.lab=2,cex.axis=2)
lines(out1[,1],out1[,2],col = "darkgrey",lwd=2)
points(d_output3[,1],d_output3[,3],pch=20,cex = 2,col='darkcyan',xlab='t',ylab=colnames(out1)[i+1])
lines(out1[,1],out1[,3],col = "darkgrey",lwd=2)
legend(13,.5,legend=c(expression('X'[1]),expression('X'[2]),'ODEs'),pch=c(20,20,NA),pt.cex = c(2,2,NA),lty=c(NA,NA,1),lwd=c(NA,NA,2),col=c('blue','darkcyan','darkgrey'),cex=1.5,bty="n")
text(0,2,labels = "B",cex=2)
# dev.off()

##### discrete with noise
######################################################################################################################################








