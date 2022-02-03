


rm(list = ls())



######################################################################################################################################
##### system of dif. eq.

if(!require('deSolve',character.only = TRUE)) {install.packages('deSolve');library(deSolve)}

t = seq(0,60,1)

InitState5 = c(
		x1 = 1,
		x2 = 4,
		x3 = .5,
		x4 = 2.5
		)

Pars5 = c(	
	x0 = 1,

	a1 = 0.25,  
	b1 = 0.08, 
	b2 = 0.2,
	b3 = 0.3, 

	c1 = 0.15,
	c2 = 0.25, 
	d1 = 0.1, 

	g1 = -2,
	h11 = 1,
	h12 = 2,
	h13 = 0.5,
	h2 = 0.4,
	h3 = 0.6,
	h4 = 0.8
      )



Equations5 <- function(t, x, pars) 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {

		dx1 = (a1 * x0 * x3^g1 - b1 * x1^h11 - c1 * x1^h12 * x3^h13) 
		dx2 = (b1 * x1^h11 - b2 * x2^h2) 
		dx3 = (b2 * x2^h2 - b3 * x3^h3) 
		dx4 = (c1 * x1^h12 * x3^h13 - c2 * x4^h4) 

            return(list(c(dx1,dx2,dx3,dx4),count=c(),signal = c()))
            })
        }


##### system of dif. eq.
######################################################################################################################################

out5 = ode(
		times = t,
		y = InitState5,
		parms = Pars5,
		func = Equations5) # events = list(data = perturb) 

#windows()
plot(out5[,1],out5[,2],type='l',col='red',lwd=3,ylim=c(0,4),xlab='t',ylab='')
lines(out5[,1],out5[,3],type='l',col='blue',lwd=3)
lines(out5[,1],out5[,4],type='l',col='green',lwd=3)
lines(out5[,1],out5[,5],type='l',col='darkcyan',lwd=3)
legend(15,4,legend=c('x1','x2','x3','x4'),lty=c(1,1,1,1),lwd=c(3,3,3,3),col=c('red','blue','green','darkcyan'),bty="n")

######################################################################################################################################
##### discrete 1

d_output5 = c(t=min(t),InitState5)
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )
	{
	if ( is.null(dim(d_output5))==TRUE )
		{
		x1 = d_output5[2]
		x2 = d_output5[3]
		x3 = d_output5[4]
		x4 = d_output5[5]
		} else 
		{
		x1 = d_output5[dim(d_output5)[1],2]
		x2 = d_output5[dim(d_output5)[1],3]
		x3 = d_output5[dim(d_output5)[1],4]
		x4 = d_output5[dim(d_output5)[1],5]
		}

	Xs = with(as.list(c(c(x1,x2,x3,x4), Pars5)), 
            {

		x1_temp = ( a1 * x0 * x3^g1 - b1 * x1^h11 - c1 * x1^h12 * x3^h13 ) * (t[2]-t[1]) + x1
		x2_temp = ( b1 * x1^h11 - b2 * x2^h2 ) * (t[2]-t[1]) + x2
		x3_temp = ( b2 * x2^h2 - b3 * x3^h3 ) * (t[2]-t[1]) + x3
		x4_temp = ( c1 * x1^h12 * x3^h13 - c2 * x4^h4 ) * (t[2]-t[1]) + x4

            return(list(c(x1_temp,x2_temp,x3_temp,x4_temp)))
            })

	d_output5 = rbind(d_output5,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2],x3=Xs[[1]][3],x3=Xs[[1]][4]))
	}


# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_lac1.tiff", height = 10, width = 20, units = 'cm',compression = "lzw", res = 300)
plot(d_output5[,1],d_output5[,2],pch=20,col='red',lwd=3,xlab='t',ylab='',ylim=c(0,4),cex.lab=1.5,cex.axis=1.5)
lines(out5[,1],out5[,2],col = "red",lwd=2)
points(d_output5[,1],d_output5[,3],pch=20,col='blue',lwd=3,xlab='t',ylab='')
lines(out5[,1],out5[,3],col = "blue",lwd=2)
points(d_output5[,1],d_output5[,4],pch=20,col='green',lwd=3,xlab='t',ylab='')
lines(out5[,1],out5[,4],col = "green",lwd=2)
points(d_output5[,1],d_output5[,5],pch=20,col='darkcyan',lwd=3,xlab='t',ylab='')
lines(out5[,1],out5[,5],col = "darkcyan",lwd=2)
legend(10,4,legend=c('x1','x2','x3','x4','ODEs'),pch=c(20,20,20,20,NA),pt.cex = c(2,2,2,2,NA),col=c('red','blue','green','darkcyan','black'),lty=c(NA,NA,NA,NA,1),lwd=c(NA,NA,NA,NA,2),bty="n")
# dev.off()

##### discrete 1
######################################################################################################################################



######################################################################################################################################
##### optimization, original data

par_est_data = read.csv('par_est_data.csv')

varcolors = c('red','blue','green','darkcyan')
plot(0,0,pch=20,cex=2,col='white',xlim=c(0,60),ylim=c(0,4),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
for(i in 2:5)
	{
	lines(out5[,1],out5[,i],col=varcolors[i-1],lwd=3)
	points(par_est_data[,1],par_est_data[,i],pch=6,cex=1,col=varcolors[i-1])
	}
legend(15,4,legend=c('x1','x2','x3','x4'),lty=c(1,1,1,1),lwd=c(3,3,3,3),col=c('red','blue','green','darkcyan'),bty="n")



t = seq(0,60,.5)

OBJ_FUN = function (initPars)
	{
	
	# calculate time course with new parameters  
	d_output_est = c(t=min(t),InitState5)
	for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )
		{
		if ( is.null(dim(d_output_est))==TRUE )
			{
			x1 = d_output_est[2]
			x2 = d_output_est[3]
			x3 = d_output_est[4]
			x4 = d_output_est[5]
			} else 
			{
			x1 = d_output_est[dim(d_output_est)[1],2]
			x2 = d_output_est[dim(d_output_est)[1],3]
			x3 = d_output_est[dim(d_output_est)[1],4]
			x4 = d_output_est[dim(d_output_est)[1],5]
			}	

		Xs = with(as.list(c(c(x1,x2,x3,x4), initPars)), 
	            {

			x1_temp = ( a1 * x0 * x3^g1 - b1 * x1^h11 - c1 * x1^h12 * x3^h13 ) * (t[2]-t[1]) + x1
			x2_temp = ( b1 * x1^h11 - b2 * x2^h2 ) * (t[2]-t[1]) + x2
			x3_temp = ( b2 * x2^h2 - b3 * x3^h3 ) * (t[2]-t[1]) + x3
			x4_temp = ( c1 * x1^h12 * x3^h13 - c2 * x4^h4 ) * (t[2]-t[1]) + x4

	            return(list(c(x1_temp,x2_temp,x3_temp,x4_temp)))
      	      })

		d_output_est = rbind(d_output_est,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2],x3=Xs[[1]][3],x3=Xs[[1]][4]))
		}

	error = 100000 * runif(1,.95,1.05)																# if out_est is the same length as the data (important to calculate errors)	
	if ( sum(is.nan(d_output_est))==0 )																				# if there is no NAs
		{
		data_est = par_est_data
		varcolors = c('red','blue','green','darkcyan')
		plot(0,0,pch=20,cex=2,col='white',xlim=c(0,60),ylim=c(0,4),xlab='t',ylab='',cex.lab=2,cex.axis=2)

		for (i in 2:(dim(par_est_data)[2]))
			{
			# do a spline of the estimates
			smoothSpline = smooth.spline(d_output_est[,1],d_output_est[,i],df=dim(d_output_est)[1])
			# get the values of the spline to the data timepoints
     			smooth = predict(object = smoothSpline, x = par_est_data[,1], deriv = 0)                                 					# get the points of the fit of that linear model
			data_est[,i] = smooth$y
			#lines(out5[,1],out5[,i],col=varcolors[i-1],lwd=1)
			points(par_est_data[,1],par_est_data[,i],pch=6,cex=1,col=varcolors[i-1])
			lines(d_output_est[,1],d_output_est[,i],lty=1,lwd=2,col=varcolors[i-1])
			#points(smooth$x,smooth$y,pch=20,col=varcolors[i-1])
			}
		error = sum( (par_est_data - data_est)^2 )
		}
	return(error)
	}
OBJ_FUN(initPars = Pars5)

set.seed(3000)
start_pars = start_pars1 = Pars5 + rnorm(length(Pars5),mean=0,sd=.18)
names(start_pars) = names(Pars5)
OBJ_FUN(initPars = start_pars)

check = TRUE
previous_error = 900000
while (check)
	{
	opt = optim(
		par = start_pars, 
		fn = OBJ_FUN,
		method = "Nelder-Mead"		# "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
	      #lower = par_lower_bound, upper = par_upper_bound
	      )
	opt
	start_pars = opt$par
	if ( abs(opt$value - previous_error) < 10^(-3) ) { check = FALSE }
	previous_error = opt$value
	flush.console()
	}


opt = optim(
	par = start_pars, 
	fn = OBJ_FUN,
	method = "Nelder-Mead"		# "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
      #lower = par_lower_bound, upper = par_upper_bound
      )
opt
# summary(opt)
# names(opt)
# start_pars = opt$par
# 
legend(10,4,legend=c('x1','x2','x3','x4','ODEs','data','disc','discPnts'),lty=c(1,1,1,1,1,NA,2,NA),lwd=c(3,3,3,3,1,NA,2,NA),pch=c(NA,NA,NA,NA,NA,6,NA,20),col=c('red','blue','green','darkcyan','black','black','black','black'),bty="n")
cbind("Initial pars" = start_pars1,"True pars" = Pars5, "Estimates" = opt$par, dif = Pars5 - opt$par)

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_par_est_orig_data_halp_step.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
OBJ_FUN(initPars = opt$par)
legend(10,4,legend=c(expression('X'[1]),expression('X'[2]),expression('X'[3]),expression('X'[4]),'data'),lty=c(1,1,1,1,NA),lwd=c(3,3,3,3,NA),pch=c(NA,NA,NA,NA,6),col=c('red','blue','green','darkcyan','black'),cex=1.5,bty="n")
text(5,4,labels = "A",cex=2)
# dev.off()

##### optimization, original data
######################################################################################################################################



######################################################################################################################################
##### optimization, reduced sparse data

par_est_data = read.csv('par_est_data.csv')
par_est_data3 = par_est_data[seq(0,60,3),]

varcolors = c('red','blue','green','darkcyan')
plot(0,0,pch=20,cex=2,col='white',xlim=c(0,60),ylim=c(0,4),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
for(i in 2:5)
	{
	lines(out5[,1],out5[,i],col=varcolors[i-1],lwd=3)
	points(par_est_data[,1],par_est_data[,i],pch=6,cex=1,col=varcolors[i-1])
	points(par_est_data3[,1],par_est_data3[,i],pch=2,cex=1,col=varcolors[i-1])
	}
legend(15,4,legend=c('x1','x2','x3','x4'),lty=c(1,1,1,1),lwd=c(3,3,3,3),col=c('red','blue','green','darkcyan'),bty="n")

t = seq(0,60,1)

OBJ_FUN = function (initPars)
	{
	
	# calculate time course with new parameters  
	d_output_est = c(t=min(t),InitState5)
	for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )
		{
		if ( is.null(dim(d_output_est))==TRUE )
			{
			x1 = d_output_est[2]
			x2 = d_output_est[3]
			x3 = d_output_est[4]
			x4 = d_output_est[5]
			} else 
			{
			x1 = d_output_est[dim(d_output_est)[1],2]
			x2 = d_output_est[dim(d_output_est)[1],3]
			x3 = d_output_est[dim(d_output_est)[1],4]
			x4 = d_output_est[dim(d_output_est)[1],5]
			}	

		Xs = with(as.list(c(c(x1,x2,x3,x4), initPars)), 
	            {

			x1_temp = ( a1 * x0 * x3^g1 - b1 * x1^h11 - c1 * x1^h12 * x3^h13 ) * (t[2]-t[1]) + x1
			x2_temp = ( b1 * x1^h11 - b2 * x2^h2 ) * (t[2]-t[1]) + x2
			x3_temp = ( b2 * x2^h2 - b3 * x3^h3 ) * (t[2]-t[1]) + x3
			x4_temp = ( c1 * x1^h12 * x3^h13 - c2 * x4^h4 ) * (t[2]-t[1]) + x4

	            return(list(c(x1_temp,x2_temp,x3_temp,x4_temp)))
      	      })

		d_output_est = rbind(d_output_est,c(t=i,x1=Xs[[1]][1],x2=Xs[[1]][2],x3=Xs[[1]][3],x3=Xs[[1]][4]))
		}

	error = 100000 * runif(1,.95,1.05)																# if out_est is the same length as the data (important to calculate errors)	
	if ( sum(is.nan(d_output_est))==0 )																				# if there is no NAs
		{
		data_est = par_est_data3
		varcolors = c('red','blue','green','darkcyan')
		plot(0,0,pch=20,cex=2,col='white',xlim=c(0,60),ylim=c(0,4),xlab='t',ylab='',cex.lab=2,cex.axis=2)

		for (i in 2:(dim(par_est_data3)[2]))
			{
			# do a spline of the estimates
			smoothSpline = smooth.spline(d_output_est[,1],d_output_est[,i],df=dim(d_output_est)[1])
			# get the values of the spline to the data timepoints
     			smooth = predict(object = smoothSpline, x = par_est_data3[,1], deriv = 0)                                 					# get the points of the fit of that linear model
			data_est[,i] = smooth$y
			#lines(out5[,1],out5[,i],col=varcolors[i-1],lwd=1)
			points(par_est_data3[,1],par_est_data3[,i],pch=6,cex=1,col=varcolors[i-1])
			lines(d_output_est[,1],d_output_est[,i],lty=1,lwd=2,col=varcolors[i-1])
			#points(smooth$x,smooth$y,pch=20,col=varcolors[i-1])
			}
		error = sum( (par_est_data3 - data_est)^2 )
		}
	return(error)
	}
OBJ_FUN(initPars = Pars5)

set.seed(3000)
start_pars = start_pars1 = Pars5 + rnorm(length(Pars5),mean=0,sd=.18)
names(start_pars) = names(Pars5)
OBJ_FUN(initPars = start_pars)

check = TRUE
previous_error = 900000
while (check)
	{
	opt = optim(
		par = start_pars, 
		fn = OBJ_FUN,
		method = "Nelder-Mead"		# "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
	      #lower = par_lower_bound, upper = par_upper_bound
	      )
	opt
	start_pars = opt$par
	if ( abs(opt$value - previous_error) < 10^(-3) ) { check = FALSE }
	previous_error = opt$value
	flush.console()
	}


opt = optim(
	par = start_pars, 
	fn = OBJ_FUN,
	method = "Nelder-Mead"		# "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
      #lower = par_lower_bound, upper = par_upper_bound
      )
opt
# summary(opt)
# names(opt)
# start_pars = opt$par
# 
legend(10,4,legend=c('x1','x2','x3','x4','ODEs','data','disc','discPnts'),lty=c(1,1,1,1,1,NA,2,NA),lwd=c(3,3,3,3,1,NA,2,NA),pch=c(NA,NA,NA,NA,NA,6,NA,20),col=c('red','blue','green','darkcyan','black','black','black','black'),bty="n")
cbind("Initial pars" = start_pars1,"True pars" = Pars5, "Estimates" = opt$par, dif = Pars5 - opt$par)

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_par_est_sparse_data_one_step.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
OBJ_FUN(initPars = opt$par)
legend(10,4,legend=c(expression('X'[1]),expression('X'[2]),expression('X'[3]),expression('X'[4]),'data'),lty=c(1,1,1,1,NA),lwd=c(3,3,3,3,NA),pch=c(NA,NA,NA,NA,6),col=c('red','blue','green','darkcyan','black'),cex=1.5,bty="n")
text(5,4,labels = "B",cex=2)
# dev.off()

##### optimization, v data
######################################################################################################################################





######################################################################################################################################
##### optimization, ODEs original data

par_est_data = read.csv('par_est_data.csv')

varcolors = c('red','blue','green','darkcyan')
plot(0,0,pch=20,cex=2,col='white',xlim=c(0,60),ylim=c(0,4),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
for(i in 2:5)
	{
	lines(out5[,1],out5[,i],col=varcolors[i-1],lwd=3)
	points(par_est_data[,1],par_est_data[,i],pch=6,cex=1,col=varcolors[i-1])
	}
legend(15,4,legend=c('x1','x2','x3','x4'),lty=c(1,1,1,1),lwd=c(3,3,3,3),col=c('red','blue','green','darkcyan'),bty="n")


OBJ_FUN = function (initPars)
	{
	
	out_est = try(ode(
		times = par_est_data[,1],
		y = InitState5,
		parms = initPars,
		func = Equations5),TRUE) # events = list(data = perturb) 

	error = 100000 * runif(1,.95,1.05)
	if (class(out_est)!="try-error")
		{ 																					# if try does not detect an error (may create warnnings)	
																				# if out_est is the same length as the data (important to calculate errors)	
		if ( sum(is.nan(out_est))==0 )																				# if there is no NAs
			{
			varcolors = c('red','blue','green','darkcyan')
			plot(0,0,pch=20,cex=2,col='white',xlim=c(0,60),ylim=c(0,4),xlab='t',ylab='',cex.lab=2,cex.axis=2)
	
			for (i in 2:(dim(par_est_data3)[2]))
				{
				#lines(out5[,1],out5[,i],col=varcolors[i-1],lwd=1)
				points(par_est_data[,1],par_est_data[,i],pch=6,cex=1,col=varcolors[i-1])
				lines(out_est[,1],out_est[,i],lty=1,lwd=2,col=varcolors[i-1])
				}
			error = sum( (par_est_data - out_est)^2 )
			}
		}
	return(error)
	}
OBJ_FUN(initPars = Pars5)

set.seed(3000)
start_pars = start_pars1 = Pars5 + rnorm(length(Pars5),mean=0,sd=.18)
names(start_pars) = names(Pars5)
OBJ_FUN(initPars = start_pars)

check = TRUE
previous_error = 900000
while (check)
	{
	opt = optim(
		par = start_pars, 
		fn = OBJ_FUN,
		method = "Nelder-Mead"		# "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
	      #lower = par_lower_bound, upper = par_upper_bound
	      )
	opt
	start_pars = opt$par
	if ( abs(opt$value - previous_error) < 10^(-3) ) { check = FALSE }
	previous_error = opt$value
	flush.console()
	}


opt = optim(
	par = start_pars, 
	fn = OBJ_FUN,
	method = "Nelder-Mead"		# "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
      #lower = par_lower_bound, upper = par_upper_bound
      )
opt
# summary(opt)
# names(opt)
# start_pars = opt$par
# 
legend(10,4,legend=c('x1','x2','x3','x4','ODEs','data','disc','discPnts'),lty=c(1,1,1,1,1,NA,2,NA),lwd=c(3,3,3,3,1,NA,2,NA),pch=c(NA,NA,NA,NA,NA,6,NA,20),col=c('red','blue','green','darkcyan','black','black','black','black'),bty="n")
cbind("Initial pars" = start_pars1,"True pars" = Pars5, "Estimates" = opt$par, dif = Pars5 - opt$par)

# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_par_est_orig_data_ODE.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
OBJ_FUN(initPars = opt$par)
legend(10,4,legend=c(expression('X'[1]),expression('X'[2]),expression('X'[3]),expression('X'[4]),'data'),lty=c(1,1,1,1,NA),lwd=c(3,3,3,3,NA),pch=c(NA,NA,NA,NA,6),col=c('red','blue','green','darkcyan','black'),cex=1.5,bty="n")
text(5,4,labels = "C",cex=2)
# dev.off()

##### optimization, ODEs original data
######################################################################################################################################




