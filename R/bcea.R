#####################################################################################################
## Define Classes & Methods
## (C) GianlucaBaio 4 January, 2012
##
bcea <- function(e,c,ref=1,interventions=NULL,Kmax=50000) UseMethod("bcea")

## Functions included
# bcea.default
# CEanalysis
# summary.bcea
# ceplane.plot
# ib.plot
# eib.plot
# ceac.plot
# evi.plot
# plot.bcea
# contour.bcea
# contour2
# sim.table

# mixedAn
# CEriskav 
#####################################################################################################


#####################################################################################################
## Default function
bcea.default <- function(e,c,ref=1,interventions=NULL,Kmax=50000) {
## Compute a Bayesian cost-effectiveness analysis of two or more interventions
## INPUTS:
## 1. Two objects (e,c). These can be directly computed in a simulation object "sim" from JAGS/BUGS, 
##    or derived by postprocessing of "sim" in R. The objects (e,c) have dimension (n.sim x number of 
##    interventions) and contain n.sim simulated values for the measures of effectiveness and costs 
##    for each intervention being compared. 
## 2. The reference intervention as a numeric value. Each intervention is a column in the matrixes e 
##    and c so if ref=1 the first column is assumed to be associated with the reference intervention. 
##    Intervention 1 is assumed the default reference. All others are considered comparators.
## 3. A string vector "interventions" including the names of the interventions. If none is provided 
##    then labels each as "intervention1",...,"interventionN".
## 4. The value Kmax which represents the maximum value for the willingness to pay parameter. If none 
##    is provided, then it is assumed Kmax=50000.
## 5. The value Ktable is the wtp choosen to produce a summary table (default value = 25000)
##
## OUTPUTS:
## Graphs & computed values for CE Plane, ICER, EIB, CEAC, EVPI 

# Calls the main function to make the analysis
he <- CEanalysis(e,c,ref,interventions,Kmax)
class(he) <- "bcea"
plot(he)
he
}
#####################################################################################################


#####################################################################################################
## Main analysis function
CEanalysis <- function(e,c,ref=1,interventions=NULL,Kmax=50000){
## Compute a Bayesian cost-effectiveness analysis of two or more interventions
## INPUTS:
## 1. Two objects (e,c). These can be directly computed in a simulation object "sim" from JAGS/BUGS, 
##    or derived by postprocessing of "sim" in R. The objects (e,c) have dimension (n.sim x number of 
##    interventions) and contain n.sim simulated values for the measures of effectiveness and costs 
##    for each intervention being compared. 
## 2. The reference intervention as a numeric value. Each intervention is a column in the matrixes e 
##    and c so if ref=1 the first column is assumed to be associated with the reference intervention. 
##    Intervention 1 is assumed the default reference. All others are considered comparators.
## 3. A string vector "interventions" including the names of the interventions. If none is provided 
##    then labels each as "intervention1",...,"interventionN".
## 4. The value Kmax which represents the maximum value for the willingness to pay parameter. If none 
##    is provided, then it is assumed Kmax=50000.
##

# Set the working directory to wherever the user is working, if not externally set
if(!exists("working.dir")){working.dir <- getwd()}

# Number of simulations & interventions analysed
n.sim <- dim(e)[1]
n.comparators <- dim(e)[2]

# Define reference & comparator intervention (different labels can be given here if available!)
if(is.null(interventions)){interventions <- paste("intervention",1:n.comparators)}
ints <- 1:n.comparators

# Define intervention i (where i can be a number in [1,...,n.comparators]) as the reference 
# and the other(s) as comparator(s). Default is the first intervention (first column of e or c)
comp <- ints[-ref]
n.comparisons <- n.comparators-1


# Compute Effectiveness & Cost differentials (wrt to reference intervention)
if(n.comparisons==1) {
	delta.e <- delta.c <- numeric()
	delta.e <- e[,ref]-e[,comp]
	delta.c <- c[,ref]-c[,comp] 
}
if(n.comparisons>1) {
	delta.e <- delta.c <- matrix(NA,n.sim,n.comparisons)
	for (i in 1:length(comp)) {
		delta.e[,i] <- e[,ref]-e[,comp[i]]
		delta.c[,i] <- c[,ref]-c[,comp[i]]
	}
}


# Compute the ICER
if(n.comparisons==1) {
	ICER <- mean(delta.c)/mean(delta.e)
}
if(n.comparisons>1) {
	ICER <- numeric()
	for (i in 1:n.comparisons) {
		ICER[i] <- mean(delta.c[,i])/mean(delta.e[,i])
	}
}



# Compute and plot CEAC & EIB
if(!exists("Kmax")){Kmax<-50000}
npoints <- 500 #1000
step <- Kmax/npoints
k <- seq(0,Kmax,step)
K <- length(k)

if(n.comparisons==1) {
	ceac <- numeric()
	ib <- matrix(NA,K,n.sim)
	for (i in 1:K) {
		ceac[i] <- sum(k[i]*delta.e - delta.c >0)/n.sim
		ib[i,1:n.sim] <- k[i]*delta.e - delta.c
	}
}
if(n.comparisons>1) { 
	ceac <- matrix(NA,K,n.comparisons)
	ib <- array(NA,c(K,n.sim,n.comparisons))
	for (i in 1:K) {
		for (j in 1:n.comparisons) {
			ceac[i,j] <- sum(k[i]*delta.e[,j] - delta.c[,j] >0)/n.sim
			ib[i,,j] <- k[i]*delta.e[,j] - delta.c[,j]
		}
	}
}


# Select the best option for each value of the willingness to pay parameter
if(n.comparisons==1) {
	eib <- apply(ib,1,mean)
	t1 <- numeric()
	tmp1 <- ifelse(eib>0,0,eib)
	for (i in 1:K){
		t1[i] <- ifelse((tmp1[i] & tmp1[i]<0),comp,0)
	}
	best <- t1
	best[best==0] <- ref
	## Finds the k for which the optimal decision changes
	check <- numeric()
	check[1] <- 0
	for (i in 2:K) {
		ifelse(best[i]==best[i-1],check[i] <- 0,check[i] <- 1)
	}
	kstar <- k[check==1]
}
if(n.comparisons>1) {
	eib <- matrix(NA,K,n.comparisons)
	for (j in 1:n.comparisons) {
		eib[,j] <- apply(ib[,,j],1,mean)
	}
	t <- tmp <- matrix(NA,K,n.comparisons)
	for (j in 1:n.comparisons) {
		tmp[,j] <- ifelse(eib[,j]>0,0,eib[,j])
	}
	for (i in 1:K) {
		for (j in 1:n.comparisons) {
			t[i,j] <- ifelse((min(tmp[i,])==tmp[i,j] & tmp[i,j]<0),comp[j],0)
		}
	}
	best <- apply(t,1,sum)
	best[best==0] <- ref
	# Finds the k for which the optimal decision changes
	check <- numeric()
	check[1] <- 0
	for (i in 2:K) {
		ifelse(best[i]==best[i-1],check[i] <- 0,check[i] <- 1)
	}
	kstar <- k[check==1]
}



# Compute EVPI 
U <- array(NA,c(n.sim,K,n.comparators))
vi <- Ustar <- ol <- matrix(NA,n.sim,K) 
for (i in 1:K) {
	for (j in 1:n.comparators) {
		U[,i,j] <- k[i]*e[,j]-c[,j]
	}
	Ustar[,i] <- apply(U[,i,],1,max)
	cmd <- paste("ol[,i] <- Ustar[,i] - U[,i,",best[i],"]",sep="")
	eval(parse(text=cmd))
	vi[,i] <- Ustar[,i] - max(apply(U[,i,],2,mean))
}
evi <- apply(ol,2,mean)





## Outputs of the function
list(
n.sim=n.sim,n.comparators=n.comparators,n.comparisons=n.comparisons,delta.e=delta.e,
delta.c=delta.c,ICER=ICER,Kmax=Kmax,k=k,ceac=ceac,ib=ib,eib=eib,kstar=kstar,
best=best,U=U,vi=vi,Ustar=Ustar,ol=ol,evi=evi,interventions=interventions,
ref=ref,comp=comp,step=step
)
}
#####################################################################################################



##############################################################################################################
## Summary of the results
summary.bcea <- function(object,wtp=25000,...) {

if(max(object$k)<wtp) {
wtp <- max(object$k)
cat(paste("NB: k (wtp) is defined in the interval [",min(object$k)," - ",wtp,"]\n",sep=""))
}
if(length(which(object$k==wtp))==0) {
stop(paste("The willingness to pay parameter is defined in the interval [0-",object$Kmax,
	"], with increments of ",object$step,"\n",sep=""))
}
ind.table <- which(object$k==wtp)
cols.u <- 1:object$n.comparators
cols.ustar <- max(cols.u)+1
cols.ib <- (cols.ustar+1):(cols.ustar+object$n.comparisons)
cols.ol <- max(cols.ib)+1
cols.vi <- cols.ol+1
n.cols <- cols.vi

Table <- matrix(NA,(object$n.sim+1),n.cols)
Table[1:object$n.sim,cols.u] <- object$U[,ind.table,]
Table[1:object$n.sim,cols.ustar] <- object$Ustar[,ind.table]
if(length(dim(object$ib))==2){Table[1:object$n.sim,cols.ib] <- object$ib[ind.table,]}
if(length(dim(object$ib))>2){Table[1:object$n.sim,cols.ib] <- object$ib[ind.table,,]}
Table[1:object$n.sim,cols.ol] <- object$ol[,ind.table]
Table[1:object$n.sim,cols.vi] <- object$vi[,ind.table]
if(length(dim(object$ib))==2){
	Table[(object$n.sim+1),] <- c(apply(object$U[,ind.table,],2,mean),mean(object$Ustar[,ind.table]),
		mean(object$ib[ind.table,]),mean(object$ol[,ind.table]),mean(object$vi[,ind.table]))	
}
if(length(dim(object$ib))>2){
	Table[(object$n.sim+1),] <- c(apply(object$U[,ind.table,],2,mean),mean(object$Ustar[,ind.table]),
		apply(object$ib[ind.table,,],2,mean),mean(object$ol[,ind.table]),mean(object$vi[,ind.table]))
}

names.cols <- c(paste("U",seq(1:object$n.comparators),sep=""),"U*",
		paste("IB",object$ref,"_",object$comp,sep=""),"OL","VI")
colnames(Table) <- names.cols

tab1 <- matrix(NA,object$n.comparators,1)
tab1[,1] <- Table[object$n.sim+1,(paste("U",seq(1:object$n.comparators),sep=""))]
colnames(tab1) <- "Expected utility"
rownames(tab1) <- object$interventions

tab2 <- matrix(NA,object$n.comparisons,3)
tab2[,1] <- Table[object$n.sim+1,paste("IB",object$ref,"_",object$comp,sep="")]
if (object$n.comparisons==1) {
tab2[,2] <- sum(Table[1:object$n.sim,paste("IB",object$ref,"_",object$comp,sep="")]>0)/object$n.sim
tab2[,3] <- object$ICER
}
if (object$n.comparisons>1) {
for (i in 1:object$n.comparisons) {
	tab2[i,2] <- sum(Table[1:object$n.sim,paste("IB",object$ref,"_",object$comp[i],sep="")]>0)/object$n.sim
	tab2[i,3] <- object$ICER[i]
}
}
colnames(tab2) <- c("EIB","CEAC","ICER")
rownames(tab2) <- paste(object$interventions[object$ref]," vs ",object$interventions[object$comp],sep="")

tab3 <- matrix(NA,1,1)
tab3[,1] <- Table[object$n.sim+1,"VI"]
rownames(tab3) <- "EVPI"
colnames(tab3) <- ""

## Prints the summary table
cat("\n")
cat("Cost-effectiveness analysis summary \n")
cat("\n")
cat(paste("Reference intervention: ",object$interventions[object$ref],"\n",sep=""))
if(object$n.comparisons==1) {
text.temp <- paste("Comparator intervention: ",object$interventions[object$comp],"\n",sep="")
cat(text.temp)
}

if(object$n.comparisons>1) {
text.temp <- paste("Comparator intervention(s): ",object$interventions[object$comp[1]],"\n",sep="")
cat(text.temp)
for (i in 2:object$n.comparisons) {
	cat(paste("                          : ",object$interventions[object$comp[i]],"\n",sep=""))
}
}
cat("\n")
if(length(object$kstar)==0){
cat(paste(object$interventions[object$best[1]]," dominates for all k in [",
	min(object$k)," - ",max(object$k),"] \n",sep=""))
}
if(length(object$kstar)==1){
cat(paste("Optimal decision: choose ",object$interventions[object$best[object$k==object$kstar-object$step]],
	" for k<",object$kstar," and ",object$interventions[object$best[object$k==object$kstar]],
	" for k>=",object$kstar,"\n",sep=""))
}
if(length(object$kstar)>1){
cat(paste("Optimal decision: choose ",object$interventions[object$best[object$k==object$kstar[1]-object$step]],
	" for k < ",object$kstar[1],"\n",sep=""))
for (i in 2:length(object$kstar)) {
	cat(paste("                         ",object$interventions[object$best[object$k==object$kstar[i]-object$step]],
		" for ",object$kstar[i-1]," <= k < ",object$kstar[i],"\n",sep=""))
}
cat(paste("                         ",object$interventions[object$best[object$k==object$kstar[length(object$kstar)]]],
	" for k >= ",object$kstar[length(object$kstar)],"\n",sep=""))
}
cat("\n\n")
cat(paste("Analysis for willingness to pay parameter k = ",wtp,"\n",sep=""))
cat("\n")
print(tab1,quote=F,digits=5,justify="center")
cat("\n")
print(tab2,quote=F,digits=5,justify="center")
cat("\n")
cat(paste("Optimal intervention (max expected utility) for k=",wtp,": ",
	object$interventions[object$best][object$k==wtp],"\n",sep=""))
print(tab3,quote=F,digits=5,justify="center")
}
##############################################################################################################



##############################################################################################################
# Produce a summary table with the results of simulations for the health economic variables of interest
sim.table <- function(he,wtp=25000) {

if(wtp>he$Kmax){wtp=he$Kmax}
if(length(which(he$k==wtp))==0) {
	stop(paste("The willingness to pay parameter is defined in the interval [0-",he$Kmax,"], 
		with increments of ",step,"\n",sep=""))
}
ind.table <- which(he$k==wtp)
cols.u <- 1:he$n.comparators
cols.ustar <- max(cols.u)+1
cols.ib <- (cols.ustar+1):(cols.ustar+he$n.comparisons)
cols.ol <- max(cols.ib)+1
cols.vi <- cols.ol+1
n.cols <- cols.vi

Table <- matrix(NA,(he$n.sim+1),n.cols)
Table[1:he$n.sim,cols.u] <- he$U[,ind.table,]
Table[1:he$n.sim,cols.ustar] <- he$Ustar[,ind.table]
if(length(dim(he$ib))==2){Table[1:he$n.sim,cols.ib] <- he$ib[ind.table,]}
if(length(dim(he$ib))>2){Table[1:he$n.sim,cols.ib] <- he$ib[ind.table,,]}
Table[1:he$n.sim,cols.ol] <- he$ol[,ind.table]
Table[1:he$n.sim,cols.vi] <- he$vi[,ind.table]
if(length(dim(he$ib))==2){
	Table[(he$n.sim+1),] <- c(apply(he$U[,ind.table,],2,mean),mean(he$Ustar[,ind.table]),
		mean(he$ib[ind.table,]),mean(he$ol[,ind.table]),mean(he$vi[,ind.table]))	
}
if(length(dim(he$ib))>2){
	Table[(he$n.sim+1),] <- c(apply(he$U[,ind.table,],2,mean),mean(he$Ustar[,ind.table]),
		apply(he$ib[ind.table,,],2,mean),mean(he$ol[,ind.table]),mean(he$vi[,ind.table]))
}

names.cols <- c(paste("U",seq(1:he$n.comparators),sep=""),"U*",paste("IB",he$ref,"_",he$comp,sep=""),"OL","VI")
colnames(Table) <- names.cols
rownames(Table) <- c(1:he$n.sim,"Average")

## Outputs of the function
list(Table=Table,names.cols=names.cols,wtp=wtp,ind.table=ind.table)
}
##############################################################################################################




##############################################################################################################
## Plots the CE Plane
ceplane.plot <- function(he,comparison=NULL,wtp=25000) {
# Forces R to avoid scientific format for graphs labels
options(scipen=10)
# Encodes characters so that the graph can be saved as ps or pdf
ps.options(encoding="CP1250")
pdf.options(encoding="CP1250")

if(he$n.comparisons==1) {
	m.e <- range(he$delta.e)[1]
	M.e <- range(he$delta.e)[2]
	m.c <- range(he$delta.c)[1]
	M.c <- range(he$delta.c)[2]
	step <- (M.e-m.e)/10
	m.e <- ifelse(m.e<0,m.e,-m.e)
	m.c <- ifelse(m.c<0,m.c,-m.c)
	x.pt <- .95*m.e
	y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
	xx <- seq(100*m.c/wtp,100*M.c/wtp,step)
	yy <- xx*wtp
	xx[1] <- ifelse(min(xx)<m.e,xx[1],2*m.e)
	yy[1] <- ifelse(min(yy)<m.c,yy[1],2*m.c)
	xx[length(xx)] <- ifelse(xx[length(xx)]<M.e,1.5*M.e,xx[length(xx)])
	plot(xx,yy,col="white",xlim=c(m.e,M.e),ylim=c(m.c,M.c),
		xlab="Effectiveness differential",ylab="Cost differential",
		main=paste("Cost effectiveness plane \n",
		he$interventions[he$ref]," vs ",he$interventions[he$comp],sep=""),axes=F)
	polygon(c(min(xx),seq(min(xx),max(xx),step),max(xx)),
		c(min(yy),wtp*seq(min(xx),max(xx),step),min(yy)),
		col="grey95",border="black")
#	polygon(c(xx,xx),c(yy,rev(yy)),col="grey95",border="black")
	axis(1); axis(2); box()
	points(he$delta.e,he$delta.c,pch=20,cex=.35,col="grey55")
	abline(h=0,col="dark grey")
	abline(v=0,col="dark grey")
	text(M.e,M.c,paste("\U2022"," ICER=",format(he$ICER,digits=6,nsmall=2),sep=""),cex=.95,pos=2,col="red")
	points(mean(he$delta.e),mean(he$delta.c),pch=20,col="red",cex=1)
	t1 <- paste("k==",format(wtp,digits=3,nsmall=2,scientific=F),sep="")
	text(x.pt,y.pt,parse(text=t1),cex=.8,pos=4)
}
if(he$n.comparisons>1 & is.null(comparison)==TRUE) { 
	# select the possible colors --- maximum number of comparisons equal to 8 
	# (can be extended, though --- and other might be chosen)
	cl <- colors()
	# choose colors: "black","blue","red","magenta","orange","purple","salmon","pink"
	### color <- cl[c(24,26,552,450,498,547,568,536)]
	color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	# gray scale
	plot(he$delta.e[,1],he$delta.c[,1],pch=20,cex=.35,xlim=range(he$delta.e),ylim=range(he$delta.c),
		xlab="Effectiveness differential",ylab="Cost differential",
		main="Cost-effectiveness plane")
	for (i in 2:he$n.comparisons) {
		points(he$delta.e[,i],he$delta.c[,i],pch=20,cex=.35,col=color[i])
	}
	abline(h=0,col="dark grey")
	abline(v=0,col="dark grey")
	text <- paste(he$interventions[he$ref]," vs ",he$interventions[he$comp])
	legend("topright",text,col=color,cex=.7,bty="n",lty=1)
}
if(he$n.comparisons>1 & is.null(comparison)==FALSE) { 
	m.e <- range(he$delta.e[,comparison])[1]
	M.e <- range(he$delta.e[,comparison])[2]
	m.c <- range(he$delta.c[,comparison])[1]
	M.c <- range(he$delta.c[,comparison])[2]
	step <- (M.e-m.e)/10
	m.e <- ifelse(m.e<0,m.e,-m.e)
	m.c <- ifelse(m.c<0,m.c,-m.c)
	x.pt <- .95*m.e
	y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
	xx <- seq(100*m.c/wtp,100*M.c/wtp,step)
	yy <- xx*wtp
	xx[which(xx==min(xx))] <- ifelse(min(xx)<m.e,xx[which(xx==min(xx))],2*m.e)
	yy[which(yy==min(yy))] <- ifelse(min(yy)<m.c,yy[which(yy==min(yy))],2*m.c)
	plot(xx,yy,col="white",xlim=c(m.e,M.e),ylim=c(m.c,M.c),
		xlab="Effectiveness differential",ylab="Cost differential",
		main=paste("Cost effectiveness plane \n",
		he$interventions[he$ref]," vs ",he$interventions[he$comp[comparison]],sep=""),axes=F)
	polygon(c(min(xx),seq(min(xx),max(xx),step),max(xx)),
		c(min(yy),wtp*seq(min(xx),max(xx),step),min(yy)),
		col="grey95",border="black")
#	polygon(c(xx,xx),c(yy,rev(yy)),col="grey95",border="black")
	axis(1); axis(2); box()
	points(he$delta.e[,comparison],he$delta.c[,comparison],pch=20,cex=.35,col="grey55")
	abline(h=0,col="dark grey")
	abline(v=0,col="dark grey")
	text(M.e,M.c,paste("\U2022"," ICER=",format(he$ICER[comparison],digits=6,nsmall=2),sep=""),cex=.95,pos=2,col="red")
	points(mean(he$delta.e[,comparison]),mean(he$delta.c[,comparison]),pch=20,col="red",cex=1)
	t1 <- paste("k==",format(wtp,digits=3,nsmall=2,scientific=F),sep="")
	text(x.pt,y.pt,parse(text=t1),cex=.8,pos=4)
}
}
##############################################################################################################



##############################################################################################################
## Plots the IB
ib.plot <- function(he,comparison=NULL,wtp=25000,bw=nbw,n=512,xlim=NULL){
# comparison controls which comparator is used when more than 2 interventions are present
# bw and n control the level of smoothness of the kernel density estimation
options(scipen=10)
if(max(he$k)<wtp) {
wtp <- max(he$k)
cat(paste("NB: k (wtp) is defined in the interval [",min(he$k)," - ",wtp,"]\n",sep=""))
}
w <- which(he$k==wtp)
if(he$n.comparisons==1) {
	nbw <- sd(he$ib[w,])/1.5
	d <- density(he$ib[w,],bw=bw,n=n)
	txt <- paste("Incremental Benefit distribution\n",he$interventions[he$ref],
		" vs ",he$interventions[he$comp],sep="")
}
if(he$n.comparisons>1) {
	if(is.null(comparison)){
		comparison <- 1
	}
	nbw <- sd(he$ib[w,,comparison])/1.5
	d <- density(he$ib[w,,comparison],bw=bw,n=n)
	txt <- paste("Incremental Benefit distribution\n",he$interventions[he$ref],
		" vs ",he$interventions[he$comp[comparison]],sep="")
}
if(is.null(xlim)==TRUE){
	xlim<-range(d$x)
}
plot(d$x,d$y,t="l",ylab="Density",xlab=expression(paste("IB(",bold(theta),")",sep="")),main=txt,axes=F,col="white",xlim=xlim)
box()
axis(1)
ypt <- .95*max(d$y)
xpt <- d$x[max(which(d$y>=ypt))]
text(xpt,ypt,parse(text=paste("p(IB(",expression(bold(theta)),"),k==",format(wtp,digits=8,nsmall=2),")",sep="")),
	cex=.85,pos=4)
xplus <- d$x[d$x>=0]
yplus <- d$y[d$x>=0]
polygon(c(0,xplus),c(0,yplus),density=20,border="white")
points(d$x,d$y,t="l")
abline(v=0,col="black")
}



##############################################################################################################
## Plots the EIB
eib.plot <- function(he) {
options(scipen=10)
if(he$n.comparisons==1) {
	plot(he$k,he$eib,t="l",xlab="Willingness to pay", ylab="EIB", main="Expected Incremental Benefit")
	abline(h=0,col="grey")
	if(length(he$kstar)>0) {
		abline(v=he$kstar,col="dark grey",lty="dotted")
		text(he$kstar,min(range(he$eib)),paste("k* = ",he$kstar,sep=""))
	}
}
if(he$n.comparisons>1) {
	# select the possible colors --- maximum number of comparisons equal to 8 
	# (can be extended, though --- and other might be chosen)
	cl <- colors()
	# choose colors: "black","blue","red","magenta","orange","purple","salmon","pink"
	### color <- cl[c(24,26,552,450,498,547,568,536)]
	color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	# gray scale
	
	plot(he$k,he$eib[,1],t="l",xlab="Willingness to pay", ylab="EIB",ylim=range(he$eib),
		main="Expected Incremental Benefit")
	for (j in 2:he$n.comparisons) {
		points(he$k,he$eib[,j],t="l",col=color[j])
	}
	abline(h=0,col="grey")
	if(length(he$kstar)>0) {
		abline(v=he$kstar,col="dark grey",lty="dotted")
		text(he$kstar,min(range(he$eib)),paste("k* = ",he$kstar,sep=""))
	}
	text <- paste(he$interventions[he$ref]," vs ",he$interventions[he$comp])
	legend("topleft",text,col=color,cex=.7,bty="n",lty=1)
}
}
##############################################################################################################



##############################################################################################################
## Plots the CEAC
ceac.plot <- function(he) {
options(scipen=10)
if(he$n.comparisons==1) {
	plot(he$k,he$ceac,t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
		ylim=c(0,1),main="Cost Effectiveness Acceptability Curve")
}
if(he$n.comparisons>1) {
	# select the possible colors --- maximum number of comparisons equal to 8 
	# (can be extended, though --- and other might be chosen)
	cl <- colors()
	# choose colors: "black","blue","red","magenta","orange","purple","salmon","pink"
	### color <- cl[c(24,26,552,450,498,547,568,536)]
	color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	# gray scale
	plot(he$k,he$ceac[,1],t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
		ylim=c(0,1),main="Cost Effectiveness Acceptability Curve")
	for (j in 2:he$n.comparisons) {
		points(he$k,he$ceac[,j],t="l",col=color[j])
	}
	text <- paste(he$interventions[he$ref]," vs ",he$interventions[he$comp])
	legend("bottomright",text,col=color,cex=.7,bty="n",lty=1)
}
}
##############################################################################################################



##############################################################################################################
## Plots the EVI
evi.plot <- function(he) {
options(scipen=10)
plot(he$k,he$evi,t="l",xlab="Willingness to pay",ylab="EVPI",
	main="Expected Value of Information")
if(length(he$kstar)==1) {
points(rep(he$kstar,3),c(-10000,he$evi[he$k==he$kstar]/2,he$evi[he$k==he$kstar]),t="l",lty=2,col="dark grey")
points(c(-10000,he$kstar/2,he$kstar),rep(he$evi[he$k==he$kstar],3),t="l",lty=2,col="dark grey")
}
if(length(he$kstar)>1) {
for (i in 1:length(he$kstar)) {
points(rep(he$kstar[i],3),c(-10000,he$evi[he$k==he$kstar[i]]/2,he$evi[he$k==he$kstar[i]]),
	t="l",lty=2,col="dark grey")
points(c(-10000,he$kstar[i]/2,he$kstar[i]),rep(he$evi[he$k==he$kstar[i]],3),t="l",lty=2,col="dark grey")
}
}
}
##############################################################################################################



##############################################################################################################
## Plots the main health economics outcomes in just one graph
plot.bcea <- function(x,...) {
options(scipen=10)
par(mfrow=c(2,2))
#CE Plane
if(x$n.comparisons==1) {
	plot(x$delta.e,x$delta.c,pch=20,cex=.35,xlim=range(x$delta.e),ylim=range(x$delta.c),
		xlab="Effectiveness differential",ylab="Cost differential",main="Cost-effectiveness plane")
	abline(h=0,col="dark grey")
	abline(v=0,col="dark grey")
}
if(x$n.comparisons>1) { 
	# select the possible colors --- maximum number of comparisons equal to 8 
	# (can be extended, though --- and other might be chosen)
	cl <- colors()
	# choose colors: "black","blue","red","magenta","orange","purple","salmon","pink"
	### color <- cl[c(24,26,552,450,498,547,568,536)]
	color <- cl[floor(seq(262,340,length.out=x$n.comparators))]	# gray scale
	plot(x$delta.e[,1],x$delta.c[,1],pch=20,cex=.35,xlim=range(x$delta.e),ylim=range(x$delta.c),
		xlab="Effectiveness differential",ylab="Cost differential",main="Cost-effectiveness plane")
	for (i in 2:x$n.comparisons) {
		points(x$delta.e[,i],x$delta.c[,i],pch=20,cex=.35,col=color[i])
	}
	abline(h=0,col="dark grey")
	abline(v=0,col="dark grey")
	text <- paste(x$interventions[x$ref]," vs ",x$interventions[x$comp])
	legend("topright",text,col=color,cex=.7,bty="n",lty=1)
}
#EIB
if(x$n.comparisons==1) {
	plot(x$k,x$eib,t="l",xlab="Willingness to pay", ylab="EIB", main="Expected Incremental Benefit")
	abline(h=0,col="grey")
	if(length(x$kstar)>0) {
		abline(v=x$kstar,col="dark grey",lty="dotted")
		text(x$kstar,min(range(x$eib)),paste("k = ",x$kstar,sep=""))
	}
}
if(x$n.comparisons>1) {
	plot(x$k,x$eib[,1],t="l",xlab="Willingness to pay", ylab="EIB",ylim=range(x$eib),
		main="Expected Incremental Benefit")
	for (j in 2:x$n.comparisons) {
		points(x$k,x$eib[,j],t="l",col=color[j])
	}
	abline(h=0,col="grey")
	if(length(x$kstar)>0) {
		abline(v=x$kstar,col="dark grey",lty="dotted")
		text(x$kstar,min(range(x$eib)),paste("k = ",x$kstar,sep=""))
	}
	legend("topleft",text,col=color,cex=.7,bty="n",lty=1)
}
#CEAC
if(x$n.comparisons==1) {
	plot(x$k,x$ceac,t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
		ylim=c(0,1),main="Cost Effectiveness \n Acceptability Curve")
}
if(x$n.comparisons>1) {
	plot(x$k,x$ceac[,1],t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
		ylim=c(0,1),main="Cost Effectiveness \n Acceptability Curve")
	for (j in 2:x$n.comparisons) {
		points(x$k,x$ceac[,j],t="l",col=color[j])
	}
	legend("bottomright",text,col=color,cex=.7,bty="n",lty=1)
}
#EVPI
plot(x$k,x$evi,t="l",xlab="Willingness to pay",ylab="EVPI",
	main="Expected Value of \n Information")
if(length(x$kstar)==1) {
points(rep(x$kstar,3),c(-10000,x$evi[x$k==x$kstar]/2,x$evi[x$k==x$kstar]),t="l",lty=2,col="dark grey")
points(c(-10000,x$kstar/2,x$kstar),rep(x$evi[x$k==x$kstar],3),t="l",lty=2,col="dark grey")
}
if(length(x$kstar)>1) {
for (i in 1:length(x$kstar)) {
points(rep(x$kstar[i],3),c(-10000,x$evi[x$k==x$kstar[i]]/2,x$evi[x$k==x$kstar[i]]),
	t="l",lty=2,col="dark grey")
points(c(-10000,x$kstar[i]/2,x$kstar[i]),rep(x$evi[x$k==x$kstar[i]],3),t="l",lty=2,col="dark grey")
}
}
}
##############################################################################################################


##############################################################################################################
## Contour plots for the cost-effectiveness plane
contour.bcea <- function(x,comparison=1,scale=0.5,levels=NULL,nlevels=4,...) {
options(scipen=10)
# comparison selects which plot should be made
# by default it is the first possible

require(MASS)
if (x$n.comparisons==1) {
density <- kde2d(x$delta.e,x$delta.c,n=300,h=c(sd(x$delta.e)/scale,sd(x$delta.c)/scale))
offset <- 1.0

p.ne <- sum(x$delta.e>0 & x$delta.c>0)/x$n.sim
p.nw <- sum(x$delta.e<=0 & x$delta.c>0)/x$n.sim
p.sw <- sum(x$delta.e<=0 & x$delta.c<=0)/x$n.sim
p.se <- sum(x$delta.e>0 & x$delta.c<=0)/x$n.sim

m.c <- range(x$delta.c)[1]; M.c <- range(x$delta.c)[2]
m.e <- range(x$delta.e)[1]; M.e <- range(x$delta.e)[2]

# Changes the range so that the plot always shows the x and y axes
ch1 <- ifelse(m.e>0,m.e<--m.e,m.e<-m.e)
ch2 <- ifelse(M.e<0,M.e<--M.e,M.e<-M.e)
ch3 <- ifelse(m.c>0,m.c<--m.c,m.c<-m.c)
ch4 <- ifelse(M.c<0,M.c<--M.c,M.c<-M.c)

plot(x$delta.e,x$delta.c,pch=20,cex=.3,col="dark grey",xlab="Effectiveness differential",
	ylab="Cost differential",main=paste("Cost effectiveness plane contour plot \n",
	x$interventions[x$ref]," vs ",x$interventions[x$comp],sep=""),xlim=c(m.e,M.e),ylim=c(m.c,M.c))
abline(h=0,col="dark grey")
abline(v=0,col="dark grey")
if (is.null(levels)==FALSE){
	# Normalise the density and use levels in the contour
	density$z <- (density$z-min(density$z))/(max(density$z)-min(density$z))
	contour(density$x,density$y,density$z,add=TRUE,levels=levels,drawlabels=TRUE)
}
if (is.null(levels)==TRUE) {
	contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE)
}
t1 <- paste("Pr(Delta[e]>0, Delta[c]>0)==",format(p.ne,digits=4,nsmall=3),sep="")
text(offset*M.e,offset*M.c,parse(text=t1),cex=.8,pos=2)
t2 <- paste("Pr(Delta[e]<=0, Delta[c]>0)==",format(p.nw,digits=4,nsmall=3),sep="")
text(offset*m.e,offset*M.c,parse(text=t2),cex=.8,pos=4)
t3 <- paste("Pr(Delta[e]<=0, Delta[c]<=0)==",format(p.sw,digits=4,nsmall=3),sep="")
text(offset*m.e,offset*m.c,parse(text=t3),cex=.8,pos=4)
t4 <- paste("Pr(Delta[e]>0, Delta[c]<=0)==",format(p.se,digits=4,nsmall=3),sep="")
text(offset*M.e,offset*m.c,parse(text=t4),cex=.8,pos=2)
}

if(x$n.comparisons>1) {
density <- kde2d(x$delta.e[,comparison],x$delta.c[,comparison],n=300,h=c(sd(x$delta.e[,comparison])/scale,sd(x$delta.c[,comparison])/scale))
offset <- 1.0

p.ne <- sum(x$delta.e[,comparison]>0 & x$delta.c[,comparison]>0)/x$n.sim
p.nw <- sum(x$delta.e[,comparison]<=0 & x$delta.c[,comparison]>0)/x$n.sim
p.sw <- sum(x$delta.e[,comparison]<=0 & x$delta.c[,comparison]<=0)/x$n.sim
p.se <- sum(x$delta.e[,comparison]>0 & x$delta.c[,comparison]<=0)/x$n.sim

m.c <- range(x$delta.c[,comparison])[1]; M.c <- range(x$delta.c[,comparison])[2]
m.e <- range(x$delta.e[,comparison])[1]; M.e <- range(x$delta.e[,comparison])[2]

# Changes the range so that the plot always shows the x and y axes
ch1 <- ifelse(m.e>0,m.e<--m.e,m.e<-m.e)
ch2 <- ifelse(M.e<0,M.e<--M.e,M.e<-M.e)
ch3 <- ifelse(m.c>0,m.c<--m.c,m.c<-m.c)
ch4 <- ifelse(M.c<0,M.c<--M.c,M.c<-M.c)

plot(x$delta.e[,comparison],x$delta.c[,comparison],pch=20,cex=.3,col="dark grey",
	xlab="Effectiveness differential",ylab="Cost differential",
	main=paste("Cost effectiveness plane contour plot \n",x$interventions[x$ref]," vs ",
		x$interventions[x$comp[comparison]],sep=""),xlim=c(m.e,M.e),ylim=c(m.c,M.c))	
abline(h=0,col="dark grey")
abline(v=0,col="dark grey")
contour(density$x,density$y,density$z,add=TRUE,drawlabels=TRUE)
if (is.null(levels)==FALSE){
	# Normalise the density and use levels in the contour
	density$z <- (density$z-min(density$z))/(max(density$z)-min(density$z))
	contour(density$x,density$y,density$z,add=TRUE,levels=levels,drawlabels=TRUE)
}
if (is.null(levels)==TRUE) {
	contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE)
}
t1 <- paste("Pr(Delta[e]>0, Delta[c]>0)==",format(p.ne,digits=4,nsmall=3),sep="")
text(offset*M.e,offset*M.c,parse(text=t1),cex=.8,pos=2)
t2 <- paste("Pr(Delta[e]<=0, Delta[c]>0)==",format(p.nw,digits=4,nsmall=3),sep="")
text(offset*m.e,offset*M.c,parse(text=t2),cex=.8,pos=4)
t3 <- paste("Pr(Delta[e]<=0, Delta[c]<=0)==",format(p.sw,digits=4,nsmall=3),sep="")
text(offset*m.e,offset*m.c,parse(text=t3),cex=.8,pos=4)
t4 <- paste("Pr(Delta[e]>0, Delta[c]<=0)==",format(p.se,digits=4,nsmall=3),sep="")
text(offset*M.e,offset*m.c,parse(text=t4),cex=.8,pos=2)
}
}



##############################################################################################################
contour2 <- function(he,wtp=25000,xl=NULL,yl=NULL) {
# Forces R to avoid scientific format for graphs labels
options(scipen=10)
# Encodes characters so that the graph can be saved as ps or pdf
ps.options(encoding="CP1250")
pdf.options(encoding="CP1250")

m.e <- range(he$delta.e)[1]
M.e <- range(he$delta.e)[2]
m.c <- range(he$delta.c)[1]
M.c <- range(he$delta.c)[2]
step <- (M.e-m.e)/10

m.e <- ifelse(m.e<0,m.e,-m.e)
m.c <- ifelse(m.c<0,m.c,-m.c)
x.pt <- .95*m.e
y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
xx <- seq(100*m.c/wtp,100*M.c/wtp,step)
yy <- xx*wtp
xx[1] <- ifelse(min(xx)<m.e,xx[1],2*m.e)
yy[1] <- ifelse(min(yy)<m.c,yy[1],2*m.c)
xx[length(xx)] <- ifelse(xx[length(xx)]<M.e,1.5*M.e,xx[length(xx)])
ifelse(is.null(xl),xl2<-c(m.e,M.e),xl2<-xl)
ifelse(is.null(yl),yl2<-c(m.c,M.c),yl2<-yl)
plot(xx,yy,col="white",xlim=xl2,ylim=yl2, ####xlim=c(m.e,M.e),ylim=c(m.c,M.c),
	xlab="Effectiveness differential",ylab="Cost differential",
	main=paste("Cost effectiveness plane \n",
	he$interventions[he$ref]," vs ",he$interventions[he$comp],sep=""),axes=F)
polygon(c(min(xx),seq(min(xx),max(xx),step),max(xx)),
	c(min(yy),wtp*seq(min(xx),max(xx),step),min(yy)),
	col="grey95",border="black")
#	polygon(c(xx,xx),c(yy,rev(yy)),col="grey95",border="black")
axis(1); axis(2); box()
points(he$delta.e,he$delta.c,pch=20,cex=.35,col="grey55")
abline(h=0,col="dark grey")
abline(v=0,col="dark grey")
text(M.e,M.c,paste("\U2022"," ICER=",format(he$ICER,digits=6,nsmall=2),sep=""),cex=.95,pos=2,col="red")
points(mean(he$delta.e),mean(he$delta.c),pch=20,col="red",cex=1)
t1 <- paste("k==",format(wtp,digits=3,nsmall=2,scientific=F),sep="")
text(x.pt,y.pt,parse(text=t1),cex=.8,pos=4)

# And then plots the contour
require(MASS)
offset <- 1.0
nlevels <- 4
scale <- 0.5

density <- kde2d(he$delta.e,he$delta.c,n=300,h=c(sd(he$delta.e)/scale,sd(he$delta.c)/scale))

m.c <- range(he$delta.c)[1]; M.c <- range(he$delta.c)[2]
m.e <- range(he$delta.e)[1]; M.e <- range(he$delta.e)[2]

# Changes the range so that the plot always shows the x and y axes
ch1 <- ifelse(m.e>0,m.e<--m.e,m.e<-m.e)
ch2 <- ifelse(M.e<0,M.e<--M.e,M.e<-M.e)
ch3 <- ifelse(m.c>0,m.c<--m.c,m.c<-m.c)
ch4 <- ifelse(M.c<0,M.c<--M.c,M.c<-M.c)

par(new=TRUE)
contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE,lwd=1.5)
}



##############################################################################################################
CEriskav <- function(he,r=NULL,comparison=1) UseMethod("CEriskav")


##############################################################################################################
CEriskav.default <- function(he,r=NULL,comparison=1) {
cr <- CEriskav.fn(he,r,comparison=1)
class(cr) <- "CEriskav"
cr
}


##############################################################################################################
CEriskav.fn <- function(he,r=NULL,comparison=1) {
### COMPARISON IS USED TO SELECT THE COMPARISON FOR WHICH THE ANALYSIS IS CARRIED OUT!!!
# Reference: Baio G, Dawid AP (2011).
# Default vector of risk aversion parameters
if(is.null(r)==TRUE){
	r <- c(0.000000000001,0.0000025,.000005)
}

# Computes expected utilities & EVPI for the risk aversion cases
K <- length(he$k)
R <- length(r)
Ur <- array(NA,c(dim(he$U),R))
Urstar <- array(NA,c(dim(he$Ustar),R))
for (i in 1:K) {
	for (l in 1:R) {
		for (j in 1:he$n.comparators) {
			Ur[,i,j,l] <- (1/r[l])*(1-exp(-r[l]*he$U[,i,j]))
		}
		Urstar[,i,l] <- apply(Ur[,i,,l],1,max)
	}
}

if (he$n.comparisons==1){
IBr <- Ur[,,he$ref,] - Ur[,,he$comp,]
}
if (he$n.comparisons>1){
IBr <- Ur[,,he$ref,] - Ur[,,he$comp[comparison],]
}

eibr <- apply(IBr,c(2,3),mean)
vir <- array(NA,c(he$n.sim,K,R))
for (i in 1:K) {
	for (l in 1:R) {
		vir[,i,l] <- Urstar[,i,l] - max(apply(Ur[,i,,l],2,mean))
	}
}
evir <- apply(vir,c(2,3),mean)

## Outputs of the function
list(
Ur=Ur,Urstar=Urstar,IBr=IBr,eibr=eibr,vir=vir,evir=evir,R=R,r=r,k=he$k
)
}
##############################################################################################################



##############################################################################################################
# Plots the EIB for the risk aversion case
plot.CEriskav <- function(x,...) {
options(scipen=10)
plot(x$k,x$eibr[,1],t="l",xlab="Willingness to pay", ylab=" ", main="EIB as a function of the risk aversion parameter",ylim=range(x$eibr))
linetype=seq(1,x$R)
for (l in 2:x$R) {
	points(x$k,x$eibr[,l],t="l",lty=linetype[l])
}
text <- paste("r = ",x$r,sep=""); 
# If the first value for r is small enough, consider it close to 0 and print the label accordingly
if (x$r[1]<1e-8) {
text[1] <- expression(r%->%0)
}
legend("topleft",text,lty=seq(1:x$R),cex=.9,box.lty=0)
abline(h=0,col="grey")

# Plots the EVPI for the risk aversion case
x11()
plot(x$k,x$evir[,1],t="l",ylim=range(x$evir),xlab="Willingness to pay",ylab=" ",main="EVPI as a function of the risk aversion parameter")
for (l in 2:x$R) {
	points(x$k,x$evir[,l],t="l",lty=linetype[l])
}
legend("topleft",text,lty=seq(1:x$R),cex=.9,box.lty=0)
abline(h=0,col="grey")
}


##############################################################################################################
mixedAn <- function(he,mkt.shares=NULL) UseMethod("mixedAn")


##############################################################################################################
mixedAn.default <- function(he,mkt.shares=NULL) {
ma <- mixedAn.fn(he,mkt.shares)
class(ma) <- "mixedAn"
ma
}


##############################################################################################################
mixedAn.fn <- function(he,mkt.shares=NULL) {
# mkt.shares is a vector of market shares for each comparators
# if no value is provided, then assumes uniform distribution
# dev is the device to which the graph should be printed
# default is x11() --- on screen. Possibilities are jpeg and postscript
# Reference: Baio G, Russo P (2009).


Ubar <- OL.star <- evi.star <- NULL
if(is.null(mkt.shares)==TRUE){
	mkt.shares <- rep(1,he$n.comparators)/he$n.comparators
}
temp <- array(NA,c(he$n.sim,length(he$k),he$n.comparators))
for (j in 1:he$n.comparators) {
	temp[,,j] <- mkt.shares[j]*he$U[,,j]
}
Ubar <- apply(temp,c(1,2),sum)
OL.star <- he$Ustar - Ubar
evi.star <- apply(OL.star,2,mean)


## Plot the EVPI and the mixed strategy
options(scipen=10)
plot(he$k,he$evi,t="l",xlab="Willingness to pay",ylab="EVPI",
	main="Expected Value of Information",ylim=range(he$evi,evi.star))
polygon(c(he$k,rev(he$k)),c(evi.star,rev(he$evi)),density=20,col="grey")
points(he$k,evi.star,t="l",col="red")
points(he$k,he$evi,t="l",col="black")
txt <- c("Optimal strategy","Mixed strategy:",
	paste("   ",he$interventions,"=",format(100*mkt.shares,digits=3,nsmall=2),"%",sep=""))
cols <- c("black","red",rep("white",length(he$interventions)))
legend("topleft",txt,col=cols,cex=.6,bty="n",lty=1)


## Outputs of the function
list(
Ubar=Ubar,OL.star=OL.star,evi.star=evi.star,k=he$k,Kmax=he$Kmax,step=he$step,
ref=he$ref,comp=he$comp,mkt.shares=mkt.shares,n.comparisons=he$n.comparisons,
interventions=he$interventions,evi=he$evi
)
}



##############################################################################################################
summary.mixedAn <- function(object,wtp=25000,...) {
if(max(object$k)<wtp) {
wtp <- max(object$k)
}
if(length(which(object$k==wtp))==0) {
stop(paste("The willingness to pay parameter is defined in the interval [0-",object$Kmax,"], 
	with increments of ",object$step,"\n",sep=""))
}

n.digits <- 2
n.small <- 2
cat("\n")
cat(paste("Analysis of mixed strategy for willingness to pay parameter k = ",
	wtp,"\n",sep=""))
cat("\n")
cat(paste("Reference intervention: ",object$interventions[object$ref]," (",
	format(100*object$mkt.shares[object$ref],digits=n.digits,nsmall=n.small),"% market share)\n",sep=""))
if(object$n.comparisons==1) {
text.temp <- paste("Comparator intervention: ",object$interventions[object$comp]," (",
	format(100*object$mkt.shares[object$comp],digits=n.digits,nsmall=n.small),"% market share)\n",sep="")
cat(text.temp)
}

if(object$n.comparisons>1) {
text.temp <- paste("Comparator intervention(s): ",object$interventions[object$comp[1]]," (",
	format(100*object$mkt.shares[object$comp[1]],digits=n.digits,nsmall=n.small),"% market share)\n",sep="")
cat(text.temp)
for (i in 2:object$n.comparisons) {
	cat(paste("                          : ",object$interventions[object$comp[i]]," (", 
	format(100*object$mkt.shares[object$comp[i]],digits=n.digits,nsmall=n.small),"% market share)\n",sep=""))
}
}
cat("\n")
cat(paste("Loss in the expected value of information = ",
	format(object$evi.star[object$k==wtp]-object$evi[object$k==wtp],digits=n.digits,nsmall=n.small),"\n",sep=""))
cat("\n")
}



##############################################################################################################
plot.mixedAn <- function(x,y.limits=NULL,...) {
## Plot the EVPI and the mixed strategy
options(scipen=10)
if(is.null(y.limits)){
y.limits=range(x$evi,x$evi.star)
}
plot(x$k,x$evi,t="l",xlab="Willingness to pay",ylab="EVPI",
	main="Expected Value of Information",ylim=y.limits)
polygon(c(x$k,rev(x$k)),c(x$evi.star,rev(x$evi)),density=20,col="grey")
points(x$k,x$evi.star,t="l",col="red")
points(x$k,x$evi,t="l",col="black")
txt <- c("Optimal strategy","Mixed strategy:",
	paste("   ",x$interventions,"=",format(100*x$mkt.shares,digits=3,nsmall=2),"%",sep=""))
cols <- c("black","red",rep("white",length(x$interventions)))
legend("topleft",txt,col=cols,cex=.6,bty="n",lty=1)
}
##############################################################################################################




##############################################################################################################
multi.ce <- function(he){
# Cost-effectiveness analysis for multiple comparison 
# Identifies the probability that each comparator is the most cost-effective as well as the
# cost-effectiveness acceptability frontier
	cl <- colors()
	# choose colors on gray scale
	color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	

	rank <- most.ce <- array(NA,c(he$n.sim,length(he$k),he$n.comparators))
	for (t in 1:he$n.comparators) {
		for (j in 1:length(he$k)) {
			rank[,j,t] <- apply(he$U[,j,]<=he$U[,j,t],1,sum)
			most.ce[,j,t] <- rank[,j,t]==he$n.comparators
		}
	}

	m.ce <- apply(most.ce,c(2,3),mean)		# Probability most cost-effective
	ceaf <- apply(m.ce,1,max)			# Cost-effectiveness acceptability frontier

# Output of the function
	list(
	m.ce=m.ce,ceaf=ceaf,n.comparators=he$n.comparators,k=he$k,interventions=he$interventions
	)
}
##############################################################################################################


mce.plot <- function(mce){
	cl <- colors()
	# choose colors on gray scale
	color <- cl[floor(seq(262,340,length.out=mce$n.comparators))]	
	plot(mce$k,mce$m.ce[,1],t="l",col=color[1],lwd=2,,lty=1,xlab="Willingness to pay",
		ylab="Probability of most cost effectiveness",ylim=c(0,1),
		main="Cost-effectiveness acceptability curve \nfor multiple comparisons")
	for (i in 2:mce$n.comparators) {
		points(mce$k,mce$m.ce[,i],t="l",col=color[i],lwd=2,lty=i)
	}
	legend("topright",mce$interventions,col=color,cex=.7,bty="n",lty=1:mce$n.comparators)
}


ceaf.plot <- function(mce){
	cl <- colors()
	# choose colors on gray scale
	color <- cl[floor(seq(262,340,length.out=mce$n.comparators))]	

	plot(mce$k,mce$ceaf,,t="l",lty=1,lwd=2,
		ylim=c(0,1),xlab="Willingness to pay",col=color[1],
		ylab="Probability of most cost effectiveness",
		main="Cost-effectiveness acceptability frontier")
}
