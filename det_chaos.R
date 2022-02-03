

rm(list = ls())


######################################################################################################################################
##### time, initState and truePars 

t = seq(0,2000,1)

state3 = c(
		X1 = 10,		
		X2 = 10,
		X3 = 30	
		) 

Pars3 = c()

# equations
equa1 = expression(.2 * (X2 - X1))
equa2 = expression(.6 * X1 - 0.02 * X2 - 0.02 * X1 * X3)
equa3 = expression(.02 * X1 * X2 - 0.05 * X3)

##### time, initState and truePars
######################################################################################################################################



######################################################################################################################################
##### system of dif. eq.

if(!require('deSolve',character.only = TRUE)) {install.packages('deSolve');library(deSolve)}

Equations2 <- function(t, x, pars) 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {
		# saturation should be like this but for each flux
		eq1 = eval(equa1)
		eq2 = eval(equa2)
		eq3 = eval(equa3)

		dX1 = eq1
		dX2 = eq2
		dX3 = eq3

            return(list(c(dX1,dX2,dX3),count=c(eq1,eq2,eq3)))
            })
        }

##### system of dif. eq.
######################################################################################################################################



######################################################################################################################################
##### ode 

ode_output3 = ode(times = t, y = state3, parms = Pars3, func = Equations2)
# tail(out) # edit(edit) # View(out)

par(mfrow=c(2,2))
for (i in 1:3)
	{
	plot(ode_output3[,1],ode_output3[,i+1],pch=20,col='grey',lwd=3,xlab='t',ylab=colnames(ode_output3)[i+1])
	}
# legend(60,1.0,legend = c('true model'),pch = c(20),col = c('grey'),bty = "n")	# ,lty = c(0,1)

##### ode
######################################################################################################################################



######################################################################################################################################
##### discrete 

d_output3 = c(t=min(t),state3)
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )
	{
	if ( is.null(dim(d_output3))==TRUE )
		{
		X1 = d_output3[2]
		X2 = d_output3[3]
		X3 = d_output3[4]
		} else 
		{
		X1 = d_output3[dim(d_output3)[1],2]
		X2 = d_output3[dim(d_output3)[1],3]
		X3 = d_output3[dim(d_output3)[1],4]
		}

	Xs = with(as.list(c(c(X1,X2,X3), Pars3)), 
            {
		
		X1_temp = eval(equa1) * (t[2]-t[1]) + X1
		X2_temp = eval(equa2) * (t[2]-t[1]) + X2
		X3_temp = eval(equa3) * (t[2]-t[1]) + X3

            return(list(c(X1_temp,X2_temp,X3_temp)))
            })

	d_output3 = rbind(d_output3,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2],X3=Xs[[1]][3]))
	}

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_chaos.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
layout(matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE))
plot(d_output3[,1],d_output3[,2],pch=20,col = "blue",xlab='t',ylab='',xlim=c(0,1000),ylim=c(-70,70),cex.lab=1.5,cex.axis=1.5)
lines(ode_output3[,1],ode_output3[,2],col='darkgrey',lwd=2)
points(d_output3[,1],d_output3[,3],pch=20,col = "green")
lines(ode_output3[,1],ode_output3[,3],col='darkgrey',lwd=2)
points(d_output3[,1],d_output3[,4],pch=20,col = "red")
lines(ode_output3[,1],ode_output3[,4],col='darkgrey',lwd=2)
legend(800,-15,legend=c(expression('X'[1]),expression('X'[2]),expression('X'[3]),'ODEs'),pch=c(20,20,20,NA),pt.cex = c(2,2,2,NA),lty=c(NA,NA,NA,1),lwd=c(NA,NA,NA,2),col=c('blue','green','red','darkgrey'),bty="n",cex=.9)
text(0,60,labels = "A",cex=2)

plot(d_output3[,3],d_output3[,4],type='l',col='darkcyan',lwd=1.5,xlim=c(-40,40),ylim=c(0,70),main = "Discrete phase plane",xlab = expression('X'[2]),ylab = expression('X'[3]),cex.lab=1.2,cex.axis=1.5)
points(d_output3[1,3],d_output3[1,4])
text(-35,65,labels = "B",cex=2)

plot(ode_output3[,3],ode_output3[,4],type='l',col='darkgrey',lwd=1.5,xlim=c(-40,40),ylim=c(0,70),main = "ODE phase plane",xlab = expression('X'[2]),ylab = expression('X'[3]),cex.lab=1.2,cex.axis=1.5)
points(d_output3[1,3],d_output3[1,4])
text(-35,65,labels = "C",cex=2)
# dev.off()

##### discrete 
######################################################################################################################################



