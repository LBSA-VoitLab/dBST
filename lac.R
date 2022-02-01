


rm(list = ls())



######################################################################################################################################
##### system of dif. eq.

#install.packages("deSolve")
library(deSolve)

t = seq(0,400,1)

InitState4=c(
		x1 = 0.43,
		x2 = 0.22,
		x3 = 0.46
		)

Pars4 = c(		
	g_34 = 1,
	x4 = 0.1	
      )


x4_data = cbind(c(0,10,60,110,160,210,260,310),c(0.1,1,0.5,0.1,1.2,0.4,1,0.1))
x4 = approxfun(x4_data,method='constant',rule=2)
plot(0:400,x4(0:400))



Equations4 <- function(t, x, pars) 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {
		
		X4 = x4(t)	

		dx1 = 0.2 * x3^2- 0.1*x1
		dx2 = 0.5 * x1 - x2
		dx3 = 0.1 * X4^(g_34) - 0.1 * x2 * x3

            return(list(c(dx1,dx2,dx3),count=c(),signal = X4))
            })
        }


##### system of dif. eq.
######################################################################################################################################

# perturb = data.frame(var = 'x1', time = 6, value = 2, method = "replace") #  method ("replace", "add", "multiply")

out4 = ode(
		times = t,
		y = InitState4,
		parms = Pars4,
		func = Equations4) # events = list(data = perturb) 

#windows()
plot(out4[,1],out4[,2],type='l',col='red',lwd=3,ylim=c(0,3),xlab='t',ylab='')
lines(out4[,1],out4[,3],type='l',col='blue',lwd=3)
lines(out4[,1],out4[,4],type='l',col='green',lwd=3)
lines(out4[,1],out4[,5],type='l',col='purple',lwd=3)
legend(15,.3,legend=c('x1','x2'),lty=c(1,1),lwd=c(3,3),col=c('blue','lightblue'),bty="n")





######################################################################################################################################
##### discrete 1

d_output4 = c(t=min(t),InitState4)
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )
	{
	if ( is.null(dim(d_output4))==TRUE )
		{
		x1 = d_output4[2]
		x2 = d_output4[3]
		x3 = d_output4[4]
		} else 
		{
		x1 = d_output4[dim(d_output4)[1],2]
		x2 = d_output4[dim(d_output4)[1],3]
		x3 = d_output4[dim(d_output4)[1],4]
		}

	Xs = with(as.list(c(c(x1,x2,x3))), 
            {

		x1_temp = ( 0.2 * x3^2- 0.1*x1 ) * (t[2]-t[1]) + x1
		x2_temp = ( 0.5 * x1 - x2 ) * (t[2]-t[1]) + x2
		x3_temp = ( 0.1 * x4(i) - 0.1 * x2 * x3 ) * (t[2]-t[1]) + x3

            return(list(c(x1_temp,x2_temp,x3_temp)))
            })

	d_output4 = rbind(d_output4,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2],x3=Xs[[1]][3]))
	}


# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_lac1.tiff", height = 10, width = 20, units = 'cm',compression = "lzw", res = 300)
plot(d_output4[,1],d_output4[,2],pch=20,col='red',lwd=3,xlab='t',ylab='',ylim=c(0,3),cex.lab=1.5,cex.axis=1.5)
lines(out4[,1],out4[,2],col = "darkgrey",lwd=2)
points(d_output4[,1],d_output4[,3],pch=20,col='blue',lwd=3,xlab='t',ylab='')
lines(out4[,1],out4[,3],col = "darkgrey",lwd=2)
points(d_output4[,1],d_output4[,4],pch=20,col='green',lwd=3,xlab='t',ylab='')
lines(out4[,1],out4[,4],col = "darkgrey",lwd=2)
lines(1:400,x4(1:400),pch=20,col='purple',lwd=3,xlab='t',ylab='')
legend(320,3,legend=c(expression('X'[1]),expression('X'[2]),expression('X'[3]),expression('X'[4]),'ODEs'),pch=c(20,20,20,NA,NA),pt.cex = c(2,2,2,NA,NA),col=c('red','blue','green','purple','darkgrey'),lty=c(NA,NA,NA,1,1),lwd=c(NA,NA,NA,3,2),bty="n")
# dev.off()

##### discrete 1
######################################################################################################################################



######################################################################################################################################
##### discrete 2

x5_data = 0
while (sum(x5_data)<400)
	{
	x5_data = c(x5_data,rexp(1,rate = 1/60))
	}
x5_data = x5_data[1:(length(x5_data)-1)]
sum(x5_data)
x5_data = cumsum(x5_data)

x5_data = cbind(x5_data,rnorm(n=length(x5_data),mean=.5,sd=.25))

x5 = approxfun(x5_data,method='constant',rule=2)
plot(0:400,x5(0:400))

d_output5 = c(t=min(t),InitState4)
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )
	{
	if ( is.null(dim(d_output5))==TRUE )
		{
		x1 = d_output5[2]
		x2 = d_output5[3]
		x3 = d_output5[4]
		} else 
		{
		x1 = d_output5[dim(d_output5)[1],2]
		x2 = d_output5[dim(d_output5)[1],3]
		x3 = d_output5[dim(d_output5)[1],4]
		}

	Xs = with(as.list(c(c(x1,x2,x3))), 
            {

		x1_temp = ( 0.2 * x3^2- 0.1*x1 ) * (t[2]-t[1]) + x1
		x2_temp = ( 0.5 * x1 - x2 ) * (t[2]-t[1]) + x2
		x3_temp = ( 0.1 * x5(i) - 0.1 * x2 * x3 ) * (t[2]-t[1]) + x3

            return(list(c(x1_temp,x2_temp,x3_temp)))
            })

	d_output5 = rbind(d_output5,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2],x3=Xs[[1]][3]))
	}


# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_lac2.tiff", height = 10, width = 20, units = 'cm',compression = "lzw", res = 300)
plot(d_output5[,1],d_output5[,2],pch=20,col='red',lwd=3,xlab='t',ylab='',ylim=c(0,3),cex.lab=1.5,cex.axis=1.5)
points(d_output4[,1],d_output4[,3],pch=20,col='blue',lwd=3,xlab='t',ylab='')
points(d_output4[,1],d_output4[,4],pch=20,col='green',lwd=3,xlab='t',ylab='')
lines(1:400,x5(1:400),pch=20,col='purple',lwd=3,xlab='t',ylab='')
legend(350,3,legend=c(expression('X'[1]),expression('X'[2]),expression('X'[3]),expression('X'[4])),pch=c(20,20,20,NA),pt.cex = 2,col=c('red','blue','green','purple'),lty=c(NA,NA,NA,1),lwd=c(NA,NA,NA,3),bty="n")
# dev.off()

##### discrete 2
######################################################################################################################################



######################################################################################################################################
##### discrete 3

x5_data = 0
while (sum(x5_data)<400)
	{
	x5_data = c(x5_data,rexp(1,rate = 1/60))
	}
x5_data = x5_data[1:(length(x5_data)-1)]
sum(x5_data)
x5_data = cumsum(x5_data)

x5_data = cbind(x5_data,rnorm(n=length(x5_data),mean=.5,sd=.25))

x5 = approxfun(x5_data,method='constant',rule=2)
plot(0:400,x5(0:400))

d_output5 = c(t=min(t),InitState4)
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )
	{
	if ( is.null(dim(d_output5))==TRUE )
		{
		x1 = d_output5[2]
		x2 = d_output5[3]
		x3 = d_output5[4]
		} else 
		{
		x1 = d_output5[dim(d_output5)[1],2]
		x2 = d_output5[dim(d_output5)[1],3]
		x3 = d_output5[dim(d_output5)[1],4]
		}

	Xs = with(as.list(c(c(x1,x2,x3))), 
            {
		if (x3<0.5) { x1_temp = ( 0.2 - 0.1*x1 ) * (t[2]-t[1]) + x1 }
		if (x3>=0.5 & x3<=0.8) { x1_temp = ( 0.5 * x3^2 - 1*x1 ) * (t[2]-t[1]) + x1 }
		if (x3>0.8) { x1_temp = ( 0.1 - 0.1*x1 ) * (t[2]-t[1]) + x1 }
		
		x2_temp = ( 0.5 * x1 - x2 ) * (t[2]-t[1]) + x2
		x3_temp = ( 0.1 * x5(i) - 0.1 * x2 * x3 ) * (t[2]-t[1]) + x3

            return(list(c(x1_temp,x2_temp,x3_temp)))
            })

	d_output5 = rbind(d_output5,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2],x3=Xs[[1]][3]))
	}


# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_lac3.tiff", height = 10, width = 20, units = 'cm',compression = "lzw", res = 300)
plot(d_output5[,1],d_output5[,2],pch=20,col='red',lwd=3,xlab='t',ylab='',ylim=c(0,3),cex.lab=1.5,cex.axis=1.5)
points(d_output4[,1],d_output4[,3],pch=20,col='blue',lwd=3,xlab='t',ylab='')
points(d_output4[,1],d_output4[,4],pch=20,col='green',lwd=3,xlab='t',ylab='')
lines(1:400,x5(1:400),pch=20,col='purple',lwd=3,xlab='t',ylab='')
legend(350,3,legend=c(expression('X'[1]),expression('X'[2]),expression('X'[3]),expression('X'[4])),pch=c(20,20,20,NA),pt.cex = 2,col=c('red','blue','green','purple'),lty=c(NA,NA,NA,1),lwd=c(NA,NA,NA,3),bty="n")
# dev.off()

##### discrete 3
######################################################################################################################################


