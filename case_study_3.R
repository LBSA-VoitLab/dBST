# case study # AhR, AhRNT, AhRR


rm(list = ls())


if(!require('deSolve',character.only = TRUE)) {install.packages('deSolve');library(deSolve)}		# install and call deSolve package


######################################################################################################################################
##### system of dif. eq. 

t = seq(0,20,.1)		# set time for simulation

InitState=c(		# initial state of the variables
		X1 = 10,	# free AhR
		X2 = 0,	# AhR + ligand
		X3 = 10,	# free AhRNT
		X4 = 0,	# free AhRR
		X5 = 0,	# AhR + ligand + AhRNt
		X6 = 0,	# AhRR + AhRNT
		X7 = 0,	# mRNA_TP
		X8 = 0,	# mRNA_AhRR
		X9 = 0	# TP
		)

Pars = c(
		k1 = 1,  
		k2 = 2, 
		k3 = 1,
		k4 = 1, 
		k5 = 1,
		k6 = 2, 
		k7 = 2, 
		k8 = 2, 
		k9= 10, 
		k10 = 10,  
		k11 = 10,  
		k12 = 10, 
		k13 = 6, 
		g = -4,
		L = 0		# ligand
		)

Lstart = 2
Lend = 12
parPerturb = function(t)
	{	
	### Create a function to alter independent variables # 
	# In this case, the presence of the ligand
	if (t<Lstart) {return(0)}
	if (t>=Lstart & t<=Lend) {return(2)}
	if (t>Lend) {return(0)}		
	}
parPerturb(t=1)
parPerturb(t=2)
parPerturb(t=3)
parPerturb(t=13)

Equations_ode <- function(t, x, pars)  		# ODE model 
        { 
        ### returns rate of change
        # t = model's time structure
        # x = model initial state
        # pars = model parameters 
	
        with(as.list(c(x, pars)), 
            {	
		L = parPerturb(t)
		
		flux1 = k1*L*X1 
		flux2 = k2*X2
		flux3 = k3*X2*X3
		flux4 = k4*X5
		flux5 = k5*X3*X4
		flux6 = k6*X6
		flux7 = k7*X8
		flux8 = k8*X4
		flux9 = k9*X5 *(X6+1)^g
		flux10 = k10*X7
		flux11 = k11*X8
		flux12 = k12*X7
		flux13 = k13*X9
		
		dX1 = -flux1 + flux2
		dX2 = flux1 - flux2 - flux3 + flux4
		dX3 = -flux3 + flux4 - flux5 + flux6 
		dX4 = -flux5 + flux6 + flux7 - flux8
		dX5 = flux3 - flux4
		dX6 = flux5 - flux6
		dX7 = flux9 - flux10
		dX8 = flux9 - flux11
		dX9 = flux12 - flux13
		
            return(list(c(dX1,dX2,dX3,dX4,dX5,dX6,dX7,dX8,dX9),count=c(flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8,flux9,flux10,flux11,flux12,flux13),signal = c()))
            })
        }

#perturb_ode = data.frame(var = 'L', time = 2, value = 2, method = "add") #  method ("replace", "add", "multiply")		# perturbation

out_ode = ode(			# run ODE solver
		times = t,
		y = InitState,
		parms = Pars,
		func = Equations_ode) # ,events = list(data = perturb_ode)
	
# visualization	 
#par(mfrow=c(1,2))
line_tickness = c(0,2,2,2,4,2,2,4,4,4)
var_color = c('white','pink','purple','brown','red','grey','cyan','green','darkgreen','blue')
plot(out_ode[,1],out_ode[,2],type='l',col='darkgreen',lwd=1,ylim=c(-1,15),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
segments(x0 = 0, y0 = -.5, x1 = 1.9, y1 = -.5, lwd = 3)
segments(x0 = 2.1, y0 = -.5, x1 = 11.9, y1 = -.5, lwd = 3)
segments(x0 = 12.1, y0 = -.5, x1 = 20, y1 = -.5, lwd = 3)
text(x = c(1,7,16), y = c(-1,-1,-1), labels = c('L=0','L=2','L=0'),cex=1)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(out_ode[,1],out_ode[,i],col=var_color[i],lwd=line_tickness[i])
	}
#plot(0,0,type='l',col='white',)
legend(3,15,
	legend = c('X1 - free AhR',
		'X2 - AhR + ligand',
		'X3 - free ARNT',
		'X4 - free AhRR',
		'X5 - AhR + ligand + ARNT',
		'X6 - AhRR + ARNT',
		'X7 - mRNA_TP',
		'X8 - mRNA_AhRR',
		'X9 - TP'
		),
	lty = c(1,1,1,1,1,1,1,1,1),
	lwd = c(2,2,2,4,2,2,4,4,4),
	col = var_color[2:10],
	cex = 1,
	bty="n")

##### system of dif. eq. Y
######################################################################################################################################



######################################################################################################################################
##### system of delayed dif. eq. X

Equations_dde <- function(t, y, pars, tau)  		# DDE model 
        { 
        ### returns rate of change
        # t = model's time structure
        # y = model initial state
        # pars = model parameters 
	  # tau = delays

        with(as.list(c(y, pars)), 
            {	
		tcount = t					# time sanity check
		tlag1 = t - tau[1]			# time lag1
		tlag2 = t - tau[2]			# time lag2
		if (tlag1 <= 0) {X5_lag = 0} else {X5_lag = lagvalue(t=tlag1,nr=5)}
		if (tlag1 <= 0) {X6_lag = 0} else {X6_lag = lagvalue(t=tlag1,nr=6)}		# making sore that the lag does not get negative time values in the initial time points
		if (tlag2 <= 0) {X7_lag = 0} else {X7_lag = lagvalue(t=tlag2,nr=7)}			
		if (tlag2 <= 0) {X8_lag = 0} else {X8_lag = lagvalue(t=tlag2,nr=8)}	

		# VERY IMPORTANT : "lagvalue" is the function that retrives values from past times. If you are using more than one variable, you have to specify that variable are you retriving in "nr"
		# If  you want the lag value of the derivative, use "lagderiv". It works in a similar way to "lagvalue".

		L = parPerturb(t)

		flux1 = k1*L*X1 
		flux2 = k2*X2
		flux3 = k3*X2*X3
		flux4 = k4*X5
		flux5 = k5*X3*X4
		flux6 = k6*X6		# replace X6_lag by X6
		flux7 = k7*X8_lag
		flux8 = k8*X4
		flux9 = k9*X5_lag *(X6+1)^g 
		flux10 = k10*X7
		flux11 = k11*X8
		flux12 = k12*X7_lag
		flux13 = k13*X9

		dX1 = -flux1 + flux2
		dX2 = flux1 - flux2 - flux3 + flux4
		dX3 = -flux3 + flux4 - flux5 + flux6 
		dX4 = -flux5 + flux6 + flux7 - flux8
		dX5 = flux3 - flux4
		dX6 = flux5 - flux6
		dX7 = flux9 - flux10
		dX8 = flux9 - flux11
		dX9 = flux12 - flux13

            return(list(dy = c(dX1,dX2,dX3,dX4,dX5,dX6,dX7,dX8,dX9),count=c(flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8,flux9,flux10,flux11,flux12,flux13),tcount=tcount,tlag=tlag1,X6_lag=X6_lag))
            })
        }

# perturb_dde = data.frame(var = 'X1', time = 3, value = 3, method = "add") #  method ("replace", "add", "multiply")		# perturbation done without rootfun

# rootfun <- function(t, y, p) {return(t - 3)}; eventfun <- function(t, y, p) {return (c(y[1] +3, y[2]))}			# perturbation with rootfun triggering the events

out_dde = dede(			# run DEDE solver
		times = t,
		y = InitState,
		parms = Pars,
		func = Equations_dde,
		tau = c(3,4)) #,events = list(data = perturb_dde)	# tau = 6
		 

		#rootfun = rootfun,
		#events = list(func = eventfun, root = TRUE)

# visualization
#par(mfrow=c(1,2))
line_tickness = c(0,2,2,2,4,2,2,4,4,4)
var_color = c('white','pink','purple','brown','red','grey','cyan','green','darkgreen','blue')
plot(out_dde[,1],out_dde[,2],type='l',col='darkgreen',lwd=1,ylim=c(-1,15),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
segments(x0 = 0, y0 = -.5, x1 = 1.9, y1 = -.5, lwd = 3)
segments(x0 = 2.1, y0 = -.5, x1 = 11.9, y1 = -.5, lwd = 3)
segments(x0 = 12.1, y0 = -.5, x1 = 20, y1 = -.5, lwd = 3)
text(x = c(1,7,16), y = c(-1,-1,-1), labels = c('L=0','L=2','L=0'),cex=1)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(out_dde[,1],out_dde[,i],col=var_color[i],lwd=line_tickness[i])
	}
#plot(0,0,type='l',col='white',)
legend(3,15,
	legend = c('X1 - free AhR',
		'X2 - AhR + ligand',
		'X3 - free ARNT',
		'X4 - free AhRR',
		'X5 - AhR + ligand + ARNT',
		'X6 - AhRR + ARNT',
		'X7 - mRNA_TP',
		'X8 - mRNA_AhRR',
		'X9 - TP'
		),
	lty = c(1,1,1,1,1,1,1,1,1),
	lwd = c(2,2,2,4,2,2,4,4,4),
	col = var_color[2:10],
	cex = 1,
	bty="n")

##### system of delayed dif. eq. X
######################################################################################################################################





######################################################################################################################################
##### discrete without anything

d_output1 = c(t=min(t),InitState[1:9]) # initial state and starting responce matrix
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )	# cycle to all timepoints
	{
	# retriving the variable values from the responce matrix
	if ( is.null(dim(d_output1))==TRUE )	# If is the first step
		{
		X1 = d_output1[2]
		X2 = d_output1[3]
		X3 = d_output1[4]
		X4 = d_output1[5]
		X5 = d_output1[6]
		X6 = d_output1[7]
		X7 = d_output1[8]
		X8 = d_output1[9]
		X9 = d_output1[10]
		} else 					# if not first step
		{
		X1 = d_output1[dim(d_output1)[1],2]
		X2 = d_output1[dim(d_output1)[1],3]
		X3 = d_output1[dim(d_output1)[1],4]
		X4 = d_output1[dim(d_output1)[1],5]
		X5 = d_output1[dim(d_output1)[1],6]
		X6 = d_output1[dim(d_output1)[1],7]
		X7 = d_output1[dim(d_output1)[1],8]
		X8 = d_output1[dim(d_output1)[1],9]
		X9 = d_output1[dim(d_output1)[1],10]
		}

	Xs = with(as.list(c(c(X1,X2,X3,X4,X5,X6,X7,X8,X9), Pars)), 
            {
		# setting the ligand (L value)
		if (i>Lstart) {L = 2} else {L = 0}
		if (i>Lend) {L = 0}

		# fluxes
		flux1 = k1*L*X1 
		flux2 = k2*X2
		flux3 = k3*X2*X3
		flux4 = k4*X5
		flux5 = k5*X3*X4
		flux6 = k6*X6
		flux7 = k7*X8
		flux8 = k8*X4
		flux9 = k9 * X5 * (X6+1)^g 	
		flux10 = k10*X7
		flux11 = k11*X8
		flux12 = k12*X7
		flux13 = k13*X9

		# discrete equations
		X1_temp = ( -flux1 + flux2 ) * (t[2]-t[1]) + X1
		X2_temp = ( flux1 - flux2 - flux3 + flux4 ) * (t[2]-t[1]) + X2
		X3_temp = ( -flux3 + flux4 - flux5 + flux6 ) * (t[2]-t[1]) + X3
		X4_temp = ( -flux5 + flux6 + flux7 - flux8 ) * (t[2]-t[1]) + X4
		X5_temp = ( flux3 - flux4 ) * (t[2]-t[1]) + X5
		X6_temp = ( flux5 - flux6 ) * (t[2]-t[1]) + X6
		X7_temp = ( flux9 - flux10 ) * (t[2]-t[1]) + X7
		X8_temp = ( flux9 - flux11 ) * (t[2]-t[1]) + X8
		X9_temp = ( flux12 - flux13 ) * (t[2]-t[1]) + X9

            return(list(c(X1_temp,X2_temp,X3_temp,X4_temp,X5_temp,X6_temp,X7_temp,X8_temp,X9_temp)))
            })
	# updating the responce matrix
	d_output1 = rbind(d_output1,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2],X3=Xs[[1]][3],X4=Xs[[1]][4],X5=Xs[[1]][5],X6=Xs[[1]][6],X7=Xs[[1]][7],X8=Xs[[1]][8],X9=Xs[[1]][9]))
	}

# visualization
# setwd('C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_caseStudy1.tiff", height = 15, width = 17, units = 'cm',compression = "lzw", res = 300)
#par(mfrow=c(1,2))
line_tickness = c(0,2,2,2,4,2,2,4,4,4)
var_color = c('white','pink','purple','brown','red','grey','cyan','green','darkgreen','blue')
plot(d_output1[,1],d_output1[,2],type='l',col='darkgreen',lwd=1,ylim=c(-1,15),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
segments(x0 = 0, y0 = -.5, x1 = 1.9, y1 = -.5, lwd = 3)
segments(x0 = 2.1, y0 = -.5, x1 = 11.9, y1 = -.5, lwd = 3)
segments(x0 = 12.1, y0 = -.5, x1 = 20, y1 = -.5, lwd = 3)
text(x = c(1,7,16), y = c(-1,-1,-1), labels = c('L=0','L=2','L=0'),cex=1)
text(x = 1, y = 14, labels = 'A',cex = 2)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(d_output1[,1],d_output1[,i],col=var_color[i],lwd=line_tickness[i])
	}
#plot(0,0,type='l',col='white',)
legend(3,15,
	legend = c('X1 - free AhR',
		'X2 - AhR + ligand',
		'X3 - free ARNT',
		'X4 - free AhRR',
		'X5 - AhR + ligand + ARNT',
		'X6 - AhRR + ARNT',
		'X7 - mRNA_TP',
		'X8 - mRNA_AhRR',
		'X9 - TP'
		),
	lty = c(1,1,1,1,1,1,1,1,1),
	lwd = c(2,2,2,4,2,2,4,4,4),
	col = var_color[2:10],
	cex = 1,
	bty="n")
# dev.off()

##### discrete without anything
######################################################################################################################################



######################################################################################################################################
##### discrete with deplay

tau1 = 3	# RNA transcription delay
tau2 = 4	# RNA translation delay
X5_store=vector()	# these are not important, just to verify is the delays were working
X7_store=vector()
X8_store=vector()
d_output2 = c(t=min(t),InitState[1:9])	# initial state and starting responce matrix
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )	# cycle to all timepoints
	{
	# retriving the variable values from the responce matrix
	if ( is.null(dim(d_output2))==TRUE )	# If is the first step
		{
		X1 = d_output2[2]
		X2 = d_output2[3]
		X3 = d_output2[4]
		X4 = d_output2[5]
		X5 = d_output2[6]
		X6 = d_output2[7]
		X7 = d_output2[8]
		X8 = d_output2[9]
		X9 = d_output2[10]
		} else 					# if not first step
		{
		X1 = d_output2[dim(d_output2)[1],2]
		X2 = d_output2[dim(d_output2)[1],3]
		X3 = d_output2[dim(d_output2)[1],4]
		X4 = d_output2[dim(d_output2)[1],5]
		X5 = d_output2[dim(d_output2)[1],6]
		X6 = d_output2[dim(d_output2)[1],7]
		X7 = d_output2[dim(d_output2)[1],8]
		X8 = d_output2[dim(d_output2)[1],9]
		X9 = d_output2[dim(d_output2)[1],10]
		}

	Xs = with(as.list(c(c(X1,X2,X3,X4,X5,X6,X7,X8,X9), Pars)), 
            {
		# setting the ligand (L value)
		if (i>Lstart) {L = 2} else {L = 0}
		if (i>Lend) {L = 0}
		# setting delay times
		temp1 = round((i - tau1 )/(t[2]-t[1]),0)
		temp2 = round((i - tau2 )/(t[2]-t[1]),0)
		# getting the delay variables values
		if (temp1 <= 0) { X5_lag = InitState[5] } else { X5_lag = d_output2[ temp1,6] } 
		X5_store <<- c(X5_store,X5_lag)
		if (temp2 <= 0) { X7_lag = InitState[7] } else { X7_lag = d_output2[ temp2,8] }
		X7_store <<- c(X7_store,X7_lag)
		if (temp2 <= 0) { X8_lag = InitState[8] } else { X8_lag = d_output2[ temp2,9] }
		X8_store <<- c(X8_store,X8_lag)

		# fluxes
		flux1 = k1*L*X1 
		flux2 = k2*X2
		flux3 = k3*X2*X3
		flux4 = k4*X5
		flux5 = k5*X3*X4
		flux6 = k6*X6
		flux7 = k7*X8_lag
		flux8 = k8*X4
		flux9 = k9 * X5_lag * (X6+1)^g 	
		flux10 = k10*X7
		flux11 = k11*X8
		flux12 = k12*X7_lag
		flux13 = k13*X9

		# discrete equations
		X1_temp = ( -flux1 + flux2 ) * (t[2]-t[1]) + X1
		X2_temp = ( flux1 - flux2 - flux3 + flux4 ) * (t[2]-t[1]) + X2
		X3_temp = ( -flux3 + flux4 - flux5 + flux6 ) * (t[2]-t[1]) + X3
		X4_temp = ( -flux5 + flux6 + flux7 - flux8 ) * (t[2]-t[1]) + X4
		X5_temp = ( flux3 - flux4 ) * (t[2]-t[1]) + X5
		X6_temp = ( flux5 - flux6 ) * (t[2]-t[1]) + X6
		X7_temp = ( flux9 - flux10 ) * (t[2]-t[1]) + X7
		X8_temp = ( flux9 - flux11 ) * (t[2]-t[1]) + X8
		X9_temp = ( flux12 - flux13 ) * (t[2]-t[1]) + X9

            return(list(c(X1_temp,X2_temp,X3_temp,X4_temp,X5_temp,X6_temp,X7_temp,X8_temp,X9_temp)))
            })
	# updating the responce matrix
	d_output2 = rbind(d_output2,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2],X3=Xs[[1]][3],X4=Xs[[1]][4],X5=Xs[[1]][5],X6=Xs[[1]][6],X7=Xs[[1]][7],X8=Xs[[1]][8],X9=Xs[[1]][9]))
	}

# visualization
# setwd('C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_caseStudy2.tiff", height = 15, width = 17, units = 'cm',compression = "lzw", res = 300)
line_tickness = c(0,2,2,2,4,2,2,4,4,4)
var_color = c('white','pink','purple','brown','red','grey','cyan','green','darkgreen','blue')
plot(d_output2[,1],d_output2[,2],type='l',col='darkgreen',lwd=1,ylim=c(-1,15),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
segments(x0 = 0, y0 = -.5, x1 = 1.9, y1 = -.5, lwd = 3)
segments(x0 = 2.1, y0 = -.5, x1 = 11.9, y1 = -.5, lwd = 3)
segments(x0 = 12.1, y0 = -.5, x1 = 20, y1 = -.5, lwd = 3)
text(x = c(1,7,16), y = c(-1,-1,-1), labels = c('L=0','L=2','L=0'),cex=1)
text(x = 1, y = 14, labels = 'B',cex = 2)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(d_output2[,1],d_output2[,i],col=var_color[i],lwd=line_tickness[i])
	}
#plot(0,0,type='l',col='white',)
legend(3,15,
	legend = c('X1 - free AhR',
		'X2 - AhR + ligand',
		'X3 - free ARNT',
		'X4 - free AhRR',
		'X5 - AhR + ligand + ARNT',
		'X6 - AhRR + ARNT',
		'X7 - mRNA_TP',
		'X8 - mRNA_AhRR',
		'X9 - TP'
		),
	lty = c(1,1,1,1,1,1,1,1,1),
	lwd = c(2,2,2,4,2,2,4,4,4),
	col = var_color[2:10],
	cex = 1,
	bty="n")
# dev.off()

##### discrete with delay
######################################################################################################################################



######################################################################################################################################
##### discrete with deplay and noise

tau1 = 3	# RNA transcription delay
tau2 = 4	# RNA translation delay
var1 = .1
X5_store=vector()	# these are not important, just to verify is the delays were working
X7_store=vector()
X8_store=vector()
d_output3 = c(t=min(t),InitState[1:9])	# initial state and starting responce matrix
for ( i in seq(min(t)+(t[2]-t[1]),max(t),t[2]-t[1]) )	# cycle to all timepoints
	{
	# retriving the variable values from the responce matrix
	if ( is.null(dim(d_output3))==TRUE )	# If is the first step
		{
		X1 = d_output3[2]
		X2 = d_output3[3]
		X3 = d_output3[4]
		X4 = d_output3[5]
		X5 = d_output3[6]
		X6 = d_output3[7]
		X7 = d_output3[8]
		X8 = d_output3[9]
		X9 = d_output3[10]
		} else  					# if not first step
		{
		X1 = d_output3[dim(d_output3)[1],2]
		X2 = d_output3[dim(d_output3)[1],3]
		X3 = d_output3[dim(d_output3)[1],4]
		X4 = d_output3[dim(d_output3)[1],5]
		X5 = d_output3[dim(d_output3)[1],6]
		X6 = d_output3[dim(d_output3)[1],7]
		X7 = d_output3[dim(d_output3)[1],8]
		X8 = d_output3[dim(d_output3)[1],9]
		X9 = d_output3[dim(d_output3)[1],10]
		}

	Xs = with(as.list(c(c(X1,X2,X3,X4,X5,X6,X7,X8,X9), Pars)), 
            {
		# setting the ligand (L value)		
		if (i>Lstart) {L = 2 * rnorm(1,1,var1)} else {L = 0}
		if (i>Lend) {L = 0}
		# setting delay times
		temp1 = round((i - tau1 )/(t[2]-t[1]),0)
		temp2 = round((i - tau2 )/(t[2]-t[1]),0)
		# getting the delay variables values
		if (temp1 <= 0) { X5_lag = InitState[5] } else { X5_lag = d_output3[ temp1,6] } 
		X5_store <<- c(X5_store,X5_lag)
		if (temp2 <= 0) { X7_lag = InitState[7] } else { X7_lag = d_output3[ temp2,8] }
		X7_store <<- c(X7_store,X7_lag)
		if (temp2 <= 0) { X8_lag = InitState[8] } else { X8_lag = d_output3[ temp2,9] }
		X8_store <<- c(X8_store,X8_lag)

		# fluxes
		flux1 = k1*L*X1 
		flux2 = k2*X2
		flux3 = k3*X2*X3
		flux4 = k4*X5
		flux5 = k5*X3*X4
		flux6 = k6*X6
		flux7 = k7*X8_lag * rnorm(1,1,var1)
		flux8 = k8*X4
		flux9 = k9 * X5_lag * (X6+1)^g 
		flux10 = k10*X7
		flux11 = k11*X8
		flux12 = k12*X7_lag * rnorm(1,1,var1)
		flux13 = k13*X9

		# discrete equations
		X1_temp = ( -flux1 + flux2 ) * (t[2]-t[1]) + X1
		X2_temp = ( flux1 - flux2 - flux3 + flux4 ) * (t[2]-t[1]) + X2
		X3_temp = ( -flux3 + flux4 - flux5 + flux6 ) * (t[2]-t[1]) + X3
		X4_temp = ( -flux5 + flux6 + flux7 - flux8 ) * (t[2]-t[1]) + X4
		X5_temp = ( flux3 - flux4 ) * (t[2]-t[1]) + X5
		X6_temp = ( flux5 - flux6 ) * (t[2]-t[1]) + X6
		X7_temp = ( flux9 * rnorm(1,1,var1) - flux10 ) * (t[2]-t[1]) + X7
		X8_temp = ( flux9 * rnorm(1,1,var1) - flux11 ) * (t[2]-t[1]) + X8
		X9_temp = ( flux12 - flux13 ) * (t[2]-t[1]) + X9

            return(list(c(X1_temp,X2_temp,X3_temp,X4_temp,X5_temp,X6_temp,X7_temp,X8_temp,X9_temp)))
            })
	# updating the responce matrix
	d_output3 = rbind(d_output3,c(t=i,X1=Xs[[1]][1],X2=Xs[[1]][2],X3=Xs[[1]][3],X4=Xs[[1]][4],X5=Xs[[1]][5],X6=Xs[[1]][6],X7=Xs[[1]][7],X8=Xs[[1]][8],X9=Xs[[1]][9]))
	}

# visualization
# setwd('C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_caseStudy3.tiff", height = 15, width = 17, units = 'cm',compression = "lzw", res = 300)
line_tickness = c(0,2,2,2,4,2,2,4,4,4)
var_color = c('white','pink','purple','brown','red','grey','cyan','green','darkgreen','blue')
plot(d_output3[,1],d_output3[,2],type='l',col='darkgreen',lwd=1,ylim=c(-1,15),xlab='t',ylab='',cex.lab=1.5,cex.axis=1.5)
segments(x0 = 0, y0 = -.5, x1 = 1.9, y1 = -.5, lwd = 3)
segments(x0 = 2.1, y0 = -.5, x1 = 11.9, y1 = -.5, lwd = 3)
segments(x0 = 12.1, y0 = -.5, x1 = 20, y1 = -.5, lwd = 3)
text(x = c(1,7,16), y = c(-1,-1,-1), labels = c('L=0','L=2','L=0'),cex=1)
text(x = 1, y = 14, labels = 'C',cex = 2)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(d_output3[,1],d_output3[,i],col=var_color[i],lwd=line_tickness[i])
	}
#plot(0,0,type='l',col='white',)
legend(3,15,
	legend = c('X1 - free AhR',
		'X2 - AhR + ligand',
		'X3 - free ARNT',
		'X4 - free AhRR',
		'X5 - AhR + ligand + ARNT',
		'X6 - AhRR + ARNT',
		'X7 - mRNA_TP',
		'X8 - mRNA_AhRR',
		'X9 - TP'
		),
	lty = c(1,1,1,1,1,1,1,1,1),
	lwd = c(2,2,2,4,2,2,4,4,4),
	col = var_color[2:10],
	cex = 1,
	bty="n")
# dev.off()

##### discrete with delay and noise
######################################################################################################################################





##############################################################################
### Big Pic

# setwd('C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\2021\\20211105_Eberhard_discrete_BST\\pictures')
# tiff("Fig_caseStudy_bigpic.tiff", height = 30, width = 15, units = 'cm',compression = "lzw", res = 300)
line_tickness = c(0,2,2,2,4,2,2,4,4,4)
var_color = c('white','pink','purple','brown','red','grey','cyan','green','darkgreen','blue')
par(mfrow = c(3,1),mai = c(.5,.5, .1, .1))
plot(d_output1[,1],d_output1[,2],type='l',col='darkgreen',lwd=1,ylim=c(0,15),axes=FALSE,xlab='',ylab='',cex.lab=1.5,cex.axis=1.5)
axis(2,col='black',col.ticks = 'black', col.axis = 'black',cex.lab=1.5,cex.axis=1.5)
text(x = 1, y = 14, labels = 'A',cex = 2)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(d_output1[,1],d_output1[,i],col=var_color[i],lwd=line_tickness[i])
	}
#plot(0,0,type='l',col='white',)
legend(3,15,
	legend = c('X1 - free AhR',
		'X2 - AhR + ligand',
		'X3 - free ARNT',
		'X4 - free AhRR',
		'X5 - AhR + ligand + ARNT',
		'X6 - AhRR + ARNT',
		'X7 - mRNA_TP',
		'X8 - mRNA_AhRR',
		'X9 - TP'
		),
	lty = c(1,1,1,1,1,1,1,1,1),
	lwd = c(2,2,2,4,2,2,4,4,4),
	col = var_color[2:10],
	cex = 1.3,
	bty="n")

plot(d_output2[,1],d_output2[,2],type='l',col='darkgreen',lwd=1,ylim=c(0,15),axes=FALSE,xlab='',ylab='',cex.lab=1.5,cex.axis=1.5)
axis(2,col='black',col.ticks = 'black', col.axis = 'black',cex.lab=1.5,cex.axis=1.5)
text(x = 1, y = 14, labels = 'B',cex = 2)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(d_output2[,1],d_output2[,i],col=var_color[i],lwd=line_tickness[i])
	}

plot(d_output3[,1],d_output3[,2],type='l',col='darkgreen',lwd=1,ylim=c(-1,15),axes=FALSE,xlab='',ylab='',cex.lab=1.5,cex.axis=1.5)
axis(1,col='black',col.ticks = 'black', col.axis = 'black',cex.lab=1.5,cex.axis=1.5)
axis(2,col='black',col.ticks = 'black', col.axis = 'black',cex.lab=1.5,cex.axis=1.5)
segments(x0 = 0, y0 = -.5, x1 = 1.9, y1 = -.5, lwd = 3)
segments(x0 = 2.1, y0 = -.5, x1 = 11.9, y1 = -.5, lwd = 3)
segments(x0 = 12.1, y0 = -.5, x1 = 20, y1 = -.5, lwd = 3)
text(x = c(1,7,16), y = c(-1,-1,-1), labels = c('L=0','L=2','L=0'),cex=1)
text(x = 1, y = 14, labels = 'C',cex = 2)
for (i in c(2,3,4,6,7,8,9,5,10))
	{
	lines(d_output3[,1],d_output3[,i],col=var_color[i],lwd=line_tickness[i])
	}

# dev.off()




















