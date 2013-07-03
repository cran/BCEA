#####################################################################################################
## Define Classes & Methods
## v1.0. 4 January, 2012
## v1.1. 14 September, 2012
## v1.2. 17 September 2012
## v1.3-0 June, 2013
## (C) GianlucaBaio
##
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
# multi.ce
# evppi, evppi.func, evppi.plot


#####################################################################################################
bcea <- function(e,c,ref=1,interventions=NULL,Kmax=50000,plot=FALSE) UseMethod("bcea")
#####################################################################################################


#####################################################################################################
## Default function
bcea.default <- function(e,c,ref=1,interventions=NULL,Kmax=50000,plot=FALSE) {
## Compute a Bayesian cost-effectiveness analysis of two or more interventions
## INPUTS:
## 1. Two objects (e,c). These can be directly computed in a simulation object "sim" from JAGS/BUGS, 
##    or derived by postprocessing of "sim" in R. The objects (e,c) have dimension (n.sim x number of 
##    interventions) and contain n.sim simulated values for the measures of effectiveness and costs 
##    for each intervention being compared. 
## 2. The reference intervention as a numeric value. Each intervention is a column in the matrices e 
##    and c so if ref=1 the first column is assumed to be associated with the reference intervention. 
##    Intervention 1 is assumed the default reference. All others are considered comparators.
## 3. A string vector "interventions" including the names of the interventions. If none is provided 
##    then labels each as "intervention1",...,"interventionN".
## 4. The value Kmax which represents the maximum value for the willingness to pay parameter. If none 
##    is provided, then it is assumed Kmax=50000.
##
## OUTPUTS:
## Graphs & computed values for CE Plane, ICER, EIB, CEAC, EVPI 

# Calls the main function to make the analysis
he <- CEanalysis(e,c,ref,interventions,Kmax)
class(he) <- "bcea"
if(plot)
  plot(he)
return(he)
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
	ICER <- apply(delta.c,2,mean)/apply(delta.e,2,mean)
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



##########################################ceplane.plot########################################################
## Plots the CE Plane
ceplane.plot <- function(he,comparison=NULL,wtp=25000,pos=c(1,1),size=NULL,graph=c("base","ggplot2"),...) {
  ### hidden options for ggplot2 ###
  # ICER.size =                    # changes ICER point size
  # label.pos = FALSE              # uses alternate position for wtp label (old specification)
  base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
  alt.legend <- pos
  
  # Forces R to avoid scientific format for graphs labels
  options(scipen=10)
  
if(base.graphics) {
  if(!is.null(size))
    message("option size will be ignored using base graphics")
  
  if(is.numeric(alt.legend)&length(alt.legend)==2){
    temp <- ""
    if(alt.legend[2]==0)
      temp <- paste0(temp,"bottom")
    else
      temp <- paste0(temp,"top")
    if(alt.legend[1]==0)
      temp <- paste0(temp,"left")
    else
      temp <- paste0(temp,"right")
    alt.legend <- temp
    if(length(grep("^(bottom|top)(left|right)$",temp))==0)
      alt.legend <- FALSE
  }
  if(is.logical(alt.legend)){
    if(!alt.legend)
      alt.legend="topright"
    else
      alt.legend="topleft"
  }
  
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
	cl <- colors()
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
	legend(alt.legend,text,col=color,cex=.7,bty="n",lty=1)
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
	xx[1] <- ifelse(min(xx)<m.e,xx[1],2*m.e)
	yy[1] <- ifelse(min(yy)<m.c,yy[1],2*m.c)
	xx[length(xx)] <- ifelse(xx[length(xx)]<M.e,1.5*M.e,xx[length(xx)])
	plot(xx,yy,col="white",xlim=c(m.e,M.e),ylim=c(m.c,M.c),
		xlab="Effectiveness differential",ylab="Cost differential",
		main=paste("Cost effectiveness plane \n",
		he$interventions[he$ref]," vs ",he$interventions[he$comp[comparison]],sep=""),axes=F)
	polygon(c(min(xx),seq(min(xx),max(xx),step),max(xx)),
		c(min(yy),wtp*seq(min(xx),max(xx),step),min(yy)),
		col="grey95",border="black")
	axis(1); axis(2); box()
	points(he$delta.e[,comparison],he$delta.c[,comparison],pch=20,cex=.35,col="grey55")
	abline(h=0,col="dark grey")
	abline(v=0,col="dark grey")
	text(M.e,M.c,paste("\U2022"," ICER=",format(he$ICER[comparison],digits=6,nsmall=2),sep=""),cex=.95,pos=2,col="red")
	points(mean(he$delta.e[,comparison]),mean(he$delta.c[,comparison]),pch=20,col="red",cex=1)
	t1 <- paste("k==",format(wtp,digits=3,nsmall=2,scientific=F),sep="")
	text(x.pt,y.pt,parse(text=t1),cex=.8,pos=4)
}
} #if(base.graphics)
else{
  if(!isTRUE(require(ggplot2)&require(grid))){
    message("falling back to base graphics\n")
    ceplane.plot(he,comparison=comparison,wtp=wtp,pos=alt.legend,graph="base"); return(invisible(NULL))
  }
  # no visible binding note
  delta.e <- delta.c <- NULL
  
  k <- he
  rm(he)
  
  if(is.null(size))
    size <- rel(3.5)
  
  label.pos <- TRUE
  opt.theme <- theme()
  ICER.size <- 2
  exArgs <- list(...)
  if(length(exArgs)>=1){
    if(exists("ICER.size",where=exArgs))
      ICER.size <- exArgs$ICER.size
    if(exists("label.pos",where=exArgs))
      if(is.logical(exArgs$label.pos))
        label.pos <- exArgs$label.pos
    for(obj in exArgs)
      if(is.theme(obj))
        opt.theme <- opt.theme + obj
  }
  
  if(k$n.comparisons==1) {
    kd <- data.frame(k$delta.e,k$delta.c)
    names(kd) <- c("delta.e","delta.c")
    # for scale_x_continuous(oob=)
    do.nothing=function(x,limits) return(x)
    # plot limits
    range.e <- range(kd$delta.e)
    range.c <- range(kd$delta.c)
    range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
    range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
    # ce plane data
    x1 <- range.e[1]-2*abs(diff(range.e))
    x2 <- range.e[2]+2*abs(diff(range.e))
    x3 <- x2
    x <- c(x1,x2,x3)
    y <- x*wtp; y[3] <- x1*wtp
    plane <- data.frame(x=x,y=y)
    
    # build a trapezoidal plane instead of a triangle if the y value is less than the minimum difference on costs
    if(y[1]>1.2*range.c[1]) {
      plane <- rbind(plane,
                     c(x2,2*range.c[1]), #new bottom-right vertex
                     c(x1,2*range.c[1])) #new bottom-left vertex
    }
    
    # actual plot
    ceplane <- ggplot(kd, aes(delta.e,delta.c)) +
      theme_bw() +
      scale_x_continuous(limits=range.e,oob=do.nothing) + scale_y_continuous(limits=range.c,oob=do.nothing) +
      scale_color_manual("",labels=paste0("ICER = ",format(k$ICER,digits=6,nsmall=2),"  "),values="red") +     
      geom_line(data=plane[1:2,],aes(x=x,y=y),color="black",linetype=1) +
      geom_polygon(data=plane,aes(x=x,y=y),fill="light gray",alpha=.3) +
      geom_hline(aes(yintercept=0),colour="grey") + geom_vline(aes(xintercept=0),colour="grey") +
      geom_point(size=1,colour="grey33") +
      geom_point(aes(mean(delta.e),mean(delta.c),color=as.factor(1)),size=ICER.size)
    
    if(!label.pos) {
      # sposta la label di wtp se la retta di wtp intercetta o no le ordinate
      ceplane <- ceplane + annotate(geom="text",x=ifelse(range.c[1]/wtp>range.e[1],range.c[1]/wtp,range.e[1]),
                                    #y=(1+.085)*range.c[1],
                                    y=range.c[1],
                                    label=paste0("k = ",format(wtp,digits=6)),hjust=-.15,size=size)
    }
    else{
      m.e <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
      m.c <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
      x.pt <- .95*m.e
      y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
      ceplane <- ceplane + annotate(geom="text",x=x.pt,y=y.pt,
                                    label=paste0("k = ",format(wtp,digits=6)),hjust=.15,size=size)
    }
  }
  
  if(k$n.comparisons>1&is.null(comparison)==TRUE) {
    
    # create dataframe for plotting
    kd <- data.frame(c(k$delta.e),c(k$delta.c))
    names(kd) <- c("delta.e","delta.c")
    kd$comparison <- as.factor(sort(rep(1:k$n.comparisons,dim(k$delta.e)[1])))
    
    # dataset for ICERs
    means.mat <- matrix(NA,nrow=k$n.comparisons,ncol=2)
    for(i in 1:k$n.comparisons)
      means.mat[i,] <- colMeans(kd[kd$comparison==i,-3])
    means <- data.frame(means.mat)
    means$comparison <- factor(1:k$n.comparisons)
    names(means) <- c("lambda.e","lambda.c","comparison")
    
    # labels for legend
    comparisons.label <- paste0(k$interventions[k$ref]," vs ",k$interventions[k$comp])
    # vector of values for color, take out white, get integer values
    colors.label <- paste0("gray",round(seq(0,100,length.out=(k$n.comparisons+1))[-(k$n.comparisons+1)]))
    
    # polygon stuff
    do.nothing=function(x,limits) return(x)
    # plot limits
    range.e <- range(kd$delta.e)
    range.c <- range(kd$delta.c)
    range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
    range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
    # ce plane data
    x1 <- range.e[1]-2*abs(diff(range.e))
    x2 <- range.e[2]+2*abs(diff(range.e))
    x3 <- x2
    x <- c(x1,x2,x3)
    y <- x*wtp; y[3] <- x1*wtp
    plane <- data.frame(x=x,y=y,comparison=factor(rep(k$n.comparisons+1,3)))
    
    # build a trapezoidal plane instead of a triangle if the y value is less than the minimum difference on costs
    if(y[1]>min(kd$delta.c)) {
      plane <- rbind(plane,
                     c(x2,2*min(kd$delta.c),k$n.comparisons+1), #new bottom-right vertex
                     c(x1,2*min(kd$delta.c),k$n.comparisons+1)) #new bottom-left vertex
    }
    
    ceplane <-
      ggplot(kd,aes(x=delta.e,y=delta.c,col=comparison)) +
      theme_bw() +
      scale_color_manual(labels=comparisons.label,values=colors.label,na.value="black") +
      scale_x_continuous(limits=range.e,oob=do.nothing) +
      scale_y_continuous(limits=range.c,oob=do.nothing) +
      annotate("line",x=plane[1:2,1],y=plane[1:2,2],colour="black") +
      annotate("polygon",plane$x,plane$y,fill="light grey",alpha=.3) +
      geom_hline(aes(yintercept=0),colour="grey") + geom_vline(aes(xintercept=0),colour="grey") +
      geom_point(size=1)
    
    # wtp label
    if(!label.pos){
      ceplane <- ceplane +
        annotate(geom="text",
                 x=ifelse(range.c[1]/wtp>range.e[1],range.c[1]/wtp,range.e[1]),
                 y=range.c[1],
                 label=paste0("k = ",format(wtp,digits=6),"  "),hjust=.15,size=size
                 )
    }
    else{
      m.e <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
      m.c <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
      x.pt <- .95*m.e
      y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
      ceplane <- ceplane + annotate(geom="text",x=x.pt,y=y.pt,
                                    label=paste0("k = ",format(wtp,digits=6)),hjust=.15,size=size)
    }
  }
  
  if(k$n.comparisons>1&is.null(comparison)==FALSE) {
    # adjusts bcea object for the correct number of dimensions and comparators
    k$comp <- k$comp[comparison]
    k$delta.e <- k$delta.e[,comparison]
    k$delta.c <- k$delta.c[,comparison]
    k$n.comparators=length(comparison)+1
    k$n.comparisons=length(comparison)
    k$interventions=k$interventions[sort(c(k$ref,k$comp))]
    k$ICER=k$ICER[comparison]
    k$ib=k$ib[,,comparison]
    k$eib=k$eib[,comparison]
    k$U=k$U[,,sort(c(k$ref,comparison+1))]
    k$ceac=k$ceac[,comparison]
    k$ref=rank(c(k$ref,k$comp))[1]
    k$comp=rank(c(k$ref,k$comp))[-1]
    k$mod <- TRUE #
    
    return(ceplane.plot(k,wtp=wtp,pos=alt.legend,graph="ggplot2",size=size,...))
  }
  
  labs.title <- "Cost-Effectiveness Plane"
  labs.title <- paste0(labs.title,
                       ifelse(k$n.comparisons==1,
                              paste0("\n",k$interventions[k$ref]," vs ",k$interventions[-k$ref]),
                              paste0(
                                ifelse(isTRUE(k$mod),
                                       paste0("\n",k$interventions[k$ref]," vs ",
                                              paste0(k$interventions[k$comp],collapse=", ")),
                                       ""))))
  
  ceplane <- ceplane + labs(title=labs.title,x="Effectiveness differential",y="Cost differential")
  
  jus <- NULL
  if(isTRUE(alt.legend)) {
    alt.legend="bottom"
    ceplane <- ceplane + theme(legend.direction="vertical")
  }
  else{
    if(is.character(alt.legend)) {
      choices <- c("left", "right", "bottom", "top")
      alt.legend <- choices[pmatch(alt.legend,choices)]
      jus="center"
      if(is.na(alt.legend))
        alt.legend=FALSE
    }
    if(length(alt.legend)>1)
      jus <- alt.legend
    if(length(alt.legend)==1 & !is.character(alt.legend)) {
      alt.legend <- c(1,1)
      jus <- alt.legend
    }
  }
  
  ceplane <- ceplane + 
    theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(face="bold")) +
  theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),legend.text.align=0) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3)) +
    opt.theme
  
  if(k$n.comparisons==1)
    ceplane  <- ceplane + theme(legend.key.size=unit(.1,"lines")) + opt.theme
  
  return(ceplane)
} # !base.graphics

}
##############################################################################################################



##############################################################################################################
## Plots the IB
ib.plot <- function(he,comparison=NULL,wtp=25000,bw=nbw,n=512,xlim=NULL,graph=c("base","ggplot2")){
  base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE)
# comparison controls which comparator is used when more than 2 interventions are present
# bw and n control the level of smoothness of the kernel density estimation
options(scipen=10)

  if(!is.null(comparison))
    stopifnot(comparison<=he$n.comparison)
  
if(base.graphics) {
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
abline(h=0,col="black")
}                   # ! base graphics
else{
  if(!isTRUE(require(ggplot2)&require(grid))) {
    message("falling back to base graphics\n")
    ib.plot(he,comparison=comparison,wtp=wtp,bw=bw,n=n,xlim=xlim,graph="base")
    return(invisible(NULL))
  }
  
  if(is.null(comparison)) {
    comparison <- 1
  }
  
  # no visible bindings note
  x <- y <- NULL
  
  if(max(he$k)<wtp) {
    wtp <- max(he$k)
    message(paste0("NB: k (wtp) is defined in the interval [",min(he$k)," - ",wtp,"]\n"))
  }
  w <- which(he$k==wtp)
  if(he$n.comparisons==1) {
    nbw <- sd(he$ib[w,])/1.5
    density <- density(he$ib[w,],bw=bw,n=n)
    df <- data.frame(x=density$x,y=density$y)
  }
  if(he$n.comparisons>1) {
    nbw <- sd(he$ib[w,,comparison])/1.5
    density <- density(he$ib[w,,comparison],bw=bw,n=n)
    df <- data.frame(x=density$x,y=density$y)
  }
  if(is.null(xlim)){
    xlim<-range(df$x)
  }
  ib <- ggplot(df,aes(x,y)) + theme_bw() +
    geom_vline(xintercept=0,colour="grey50",size=0.5) + geom_hline(yintercept=0,colour="grey50",size=0.5) +
    geom_ribbon(data=subset(df,x>0),aes(ymax=y),ymin=0,fill="grey50",alpha=.2) + geom_line() +
    annotate(geom="text",label=paste0("p(IB(theta),k=",wtp,")"),parse=T,x=df$x[which.max(df$y)],y=max(df$y),hjust=-.5,vjust=1,size=3.5) + coord_cartesian(xlim=xlim)
  
  labs.title <- paste0("Incremental Benefit Distribution\n",he$interventions[he$ref]," vs ",
                       he$interventions[he$comp[comparison]],"")
  
  ib <- ib + 
    theme(plot.title=element_text(face="bold"),text=element_text(size=11),panel.grid=element_blank()) +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    labs(title=labs.title,x=parse(text="IB(theta)"),y="Density") +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
  return(ib)
  
} #! base.graphics
}



##############################################################################################################
## Plots the EIB
eib.plot <- function(he,comparison=NULL,pos=c(1,0),size=NULL,graph=c("base","ggplot2"),...) {
  
options(scipen=10)
alt.legend <- pos
base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 


if(base.graphics) {
  if(!is.null(comparison))
    message("option comparison will be ignored using base graphics")
  if(!is.null(size))
    message("option size will be ignored using bsae graphics")
  
  if(is.numeric(alt.legend)&length(alt.legend)==2){
    temp <- ""
    if(alt.legend[2]==0)
      temp <- paste0(temp,"bottom")
    else
      temp <- paste0(temp,"top")
    if(alt.legend[1]==1)
      temp <- paste0(temp,"right")
    else
      temp <- paste0(temp,"left")
    alt.legend <- temp
    if(length(grep("^(bottom|top)(left|right)$",temp))==0)
      alt.legend <- FALSE
  }
  if(is.logical(alt.legend)){
    if(!alt.legend)
      alt.legend="topleft"
    else
      alt.legend="topright"
  }
  
if(he$n.comparisons==1) {
	plot(he$k,he$eib,t="l",xlab="Willingness to pay", ylab="EIB", main="Expected Incremental Benefit")
	abline(h=0,col="grey")
	if(length(he$kstar)>0) {
		abline(v=he$kstar,col="dark grey",lty="dotted")
		text(he$kstar,min(range(he$eib)),paste("k* = ",he$kstar,sep=""))
	}
}
if(he$n.comparisons>1) {
	color <- rep(1,he$n.comparisons); lwd <- 1
	if (he$n.comparisons>6) {
		cl <- colors()
		color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	# gray scale
		lwd <- 1.5
	}
	
	plot(he$k,he$eib[,1],t="l",xlab="Willingness to pay", ylab="EIB",ylim=range(he$eib),
		main="Expected Incremental Benefit",lty=1,lwd=lwd)
	for (j in 2:he$n.comparisons) {
		points(he$k,he$eib[,j],t="l",col=color[j],lty=j,lwd=lwd)
	}
	abline(h=0,col="grey")
	if(length(he$kstar)>0) {
		abline(v=he$kstar,col="dark grey",lty="dotted")
		text(he$kstar,min(range(he$eib)),paste("k* = ",he$kstar,sep=""))
	}
	text <- paste(he$interventions[he$ref]," vs ",he$interventions[he$comp])
	legend(alt.legend,text,col=color,cex=.7,bty="n",lty=1:he$n.comparisons)
}
} # base.graphics
else{
  if(!isTRUE(require(ggplot2)&require(grid))){
    message("falling back to base graphics\n")
    eib.plot(he,pos=alt.legend,graph="base"); return(invisible(NULL))}
  
  # no visible binding note
  k.kstar <- NULL
  
  k <- he
  rm(he)
  
  if(is.null(size))
    size=rel(3.5)
  
  opt.theme <- theme()
  exArgs <- list(...)
  for(obj in exArgs)
    if(is.theme(obj))
      opt.theme <- opt.theme + obj
  
  if(k$n.comparisons==1) {
    # data frame
    data.psa <- data.frame(k$k,k$eib)
    names(data.psa) <- c("k","eib")
    
    if(!length(k$kstar)==0) {
      # label
      label <- paste0("k* = ",format(k$kstar,digits=6))
      eib <- ggplot(data.psa, aes(k,eib)) + geom_line() + theme_bw() +
        geom_hline(aes(yintercept=0),colour="grey") + 
        geom_vline(aes(xintercept=k.kstar),data=data.frame(k$kstar),colour="grey50",linetype=2,size=.5) +
        annotate("text",label=label,x=k$kstar,y=min(k$eib),hjust=ifelse((max(k$k)-k$kstar)/max(k$k)>1/6,-.1,1.1),size=size)
    }
    else {
      eib <- ggplot(data.psa, aes(k,eib)) + geom_line() + theme_bw() +geom_hline(aes(yintercept=0),colour="grey")
    }
  }
  
  if(k$n.comparisons>1&is.null(comparison)==TRUE) {
    data.psa <- data.frame(c(k$k),c(k$eib))
    names(data.psa) <- c("k","eib")
    data.psa$comparison <- as.factor(sort(rep(1:k$n.comparisons,length(k$k))))
    
    # labels for legend
    comparisons.label <- paste0(k$interventions[k$ref]," vs ",k$interventions[k$comp])
    
    # linetype is the indicator of the comparison.
    # 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
    linetypes <- rep(c(1,2,3,4,5,6),ceiling(k$n.comparisons/6))[1:k$n.comparisons]
    
    eib <- 
      ggplot(data.psa,aes(k,eib,linetype=comparison)) + geom_line() + theme_bw() +
      geom_hline(yintercept=0,linetype=1,color="grey") + 
      scale_linetype_manual("",labels=comparisons.label,values=linetypes)
    
    if(!length(k$kstar)==0) {
      # label
      label <- paste0("k* = ",format(k$kstar,digits=6))
      
      eib <-eib +
      geom_vline(aes(xintercept=k.kstar),data=data.frame(k$kstar),colour="grey50",linetype=2,size=.5) + 
        annotate("text",label=label,x=k$kstar,y=min(k$eib),hjust=ifelse((max(k$k)-k$kstar)/max(k$k)>1/6,-.1,1.1),size=size,vjust=1)
    }
  }
  
  if(k$n.comparisons>1&is.null(comparison)==FALSE) {
    # adjusts bcea object for the correct number of dimensions and comparators
    k$comp <- k$comp[comparison]
    k$delta.e <- k$delta.e[,comparison]
    k$delta.c <- k$delta.c[,comparison]
    k$n.comparators=length(comparison)+1
    k$n.comparisons=length(comparison)
    k$interventions=k$interventions[sort(c(k$ref,k$comp))]
    k$ICER=k$ICER[comparison]
    k$ib=k$ib[,,comparison]
    k$eib=k$eib[,comparison]
    k$U=k$U[,,sort(c(k$ref,comparison+1))]
    k$ceac=k$ceac[,comparison]
    k$ref=rank(c(k$ref,k$comp))[1]
    k$comp=rank(c(k$ref,k$comp))[-1]
    k$mod <- TRUE #
    
    return(eib.plot(k,pos=alt.legend,graph="ggplot2",size=size,...))
  }
  
  eib <- eib + labs(title="Expected Incremental Benefit",x="Willingness to pay",y="EIB")
  
  jus <- NULL
  if(isTRUE(alt.legend)) {
    alt.legend="bottom"
    eib <- eib + theme(legend.direction="vertical")
  }
  else{
    if(is.character(alt.legend)) {
      choices <- c("left", "right", "bottom", "top")
      alt.legend <- choices[pmatch(alt.legend,choices)]
      jus="center"
      if(is.na(alt.legend))
        alt.legend=FALSE
    }
    if(length(alt.legend)>1)
      jus <- alt.legend
    if(length(alt.legend)==1 & !is.character(alt.legend)) {
      alt.legend <- c(1,0)
      jus <- alt.legend
    }
  }
  eib <- eib + theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank()) + 
    theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),legend.text.align=0) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3)) +
    opt.theme
  
  return(eib)
} # !base.graphics
}
##############################################################################################################



##############################################################################################################
## Plots the CEAC
ceac.plot <- function(he,comparison=NULL,pos=c(1,0),graph=c("base","ggplot2")) {
options(scipen=10)

alt.legend <- pos
base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 

if(base.graphics) {
  
  if(!is.null(comparison))
    message("option comparison will be ignored in base graphics\n")
  
    if(is.numeric(alt.legend)&length(alt.legend)==2){
      temp <- ""
      if(alt.legend[2]==1)
        temp <- paste0(temp,"top")
      else
        temp <- paste0(temp,"bottom")
      if(alt.legend[1]==0)
        temp <- paste0(temp,"left")
      else
        temp <- paste0(temp,"right")
      alt.legend <- temp
      if(length(grep("^(bottom|top)(left|right)$",temp))==0)
        alt.legend <- FALSE
    }
    if(is.logical(alt.legend)){
      if(!alt.legend)
        alt.legend="bottomright"
      else
        alt.legend="bottomleft"
    }
  
if(he$n.comparisons==1) {
	plot(he$k,he$ceac,t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
		ylim=c(0,1),main="Cost Effectiveness Acceptability Curve")
}
if(he$n.comparisons>1) {
	color <- rep(1,he$n.comparisons); lwd <- 1
	if (he$n.comparisons>6) {
		cl <- colors()
		color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	# gray scale
		lwd <- 1.5
	}

	plot(he$k,he$ceac[,1],t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
		ylim=c(0,1),main="Cost Effectiveness Acceptability Curve",lty=1,lwd=lwd)
	for (j in 2:he$n.comparisons) {
		points(he$k,he$ceac[,j],t="l",col=color[j],lty=j,lwd=lwd)
	}
	text <- paste(he$interventions[he$ref]," vs ",he$interventions[he$comp])
	legend(alt.legend,text,col=color,cex=.7,bty="n",lty=1:he$n.comparisons)
}
} # base.graphics
else{
  if(!isTRUE(require(ggplot2)&require(grid))){
    message("falling back to base graphics\n")
    ceac.plot(he,pos=alt.legend,graph="base"); return(invisible(NULL))
  }
  
  ### no visible binding problems
  k.k <- k.ceac <- NULL
  
  k <- he
  rm(he)
  if(k$n.comparisons==1) {
    data.psa <- data.frame(k$k,k$ceac)
    
    ceac <- ggplot(data.psa, aes(k.k,k.ceac)) + geom_line() 
  }
  if(k$n.comparisons>1 & is.null(comparison)==TRUE) {
    data.psa <- data.frame(c(k$k),c(k$ceac))
    names(data.psa) <- c("k","ceac")
    data.psa$comparison <- as.factor(sort(rep(1:k$n.comparisons,length(k$k))))
    
    # labels for legend
    comparisons.label <- paste0(k$interventions[k$ref]," vs ",k$interventions[k$comp])
    
    # linetype is the indicator
    linetypes <- rep(c(1,2,3,4,5,6),ceiling(k$n.comparisons/6))[1:k$n.comparisons]
    
    ceac <- ggplot(data.psa,aes(k,ceac,linetype=comparison)) +
      geom_line() +
      scale_linetype_manual("",labels=comparisons.label,values=linetypes)
  }
  
  if(k$n.comparisons>1&is.null(comparison)==FALSE) {
    # adjusts bcea object for the correct number of dimensions and comparators
    k$comp <- k$comp[comparison]
    k$delta.e <- k$delta.e[,comparison]
    k$delta.c <- k$delta.c[,comparison]
    k$n.comparators=length(comparison)+1
    k$n.comparisons=length(comparison)
    k$interventions=k$interventions[sort(c(k$ref,k$comp))]
    k$ICER=k$ICER[comparison]
    k$ib=k$ib[,,comparison]
    k$eib=k$eib[,comparison]
    k$U=k$U[,,sort(c(k$ref,comparison+1))]
    k$ceac=k$ceac[,comparison]
    k$ref=rank(c(k$ref,k$comp))[1]
    k$comp=rank(c(k$ref,k$comp))[-1]
    k$mod <- TRUE #
    
    return(ceac.plot(k,pos=alt.legend,graph="ggplot2"))
  }
  
  ceac <- ceac + theme_bw() + 
    scale_y_continuous(limits=c(0,1)) +
    labs(title="Cost-Effectiveness Acceptability Curve",x="Willingness to pay",y="Probability of cost-effectiveness") 
  
  jus <- NULL
  if(isTRUE(alt.legend)) {
    alt.legend="bottom"
    ceac <- ceac + theme(legend.direction="vertical")
  }
  else{
    if(is.character(alt.legend)) {
      choices <- c("left", "right", "bottom", "top")
      alt.legend <- choices[pmatch(alt.legend,choices)]
      jus="center"
      if(is.na(alt.legend))
        alt.legend=FALSE
    }
    if(length(alt.legend)>1)
      jus <- alt.legend
    if(length(alt.legend)==1 & !is.character(alt.legend)) {
      alt.legend <- c(1,0)
      jus <- alt.legend
    }
  }
  
  ceac <- ceac + 
    theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank()) +
    theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),legend.text.align=0) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
  
  return(ceac)
  
} # !base.graphics

}
##############################################################################################################



##############################################################################################################
## Plots the EVI
evi.plot <- function(he,graph=c("base","ggplot2")) {
options(scipen=10)
base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 

if(base.graphics){
  
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
} # base.graphics
else{
  if(!isTRUE(require(ggplot2)&require(grid))){
    message("falling back to base graphics\n")
    evi.plot(he,graph="base"); return(invisible(NULL))
  }
  
  k <- he
  rm(he)
  
  data.psa <- data.frame(c(k$k),c(k$evi))
  names(data.psa) <- c("k","evi")
  
  evi <- ggplot(data.psa, aes(k,evi)) + geom_line() + theme_bw() +
    labs(title="Expected Value of Information",x="Willingness to pay",y="EVPI") +
    theme(plot.title=element_text(face="bold"))
  
  if(length(k$kstar)!=0) {
    kstars=length(k$kstar)
    evi.at.kstar <- numeric(kstars)
    for(i in 1:kstars) {
      evi.at.kstar[i] <- k$evi[which.min(abs(k$k-k$kstar[i]))]
    }
    
    for(i in 1:kstars) {
      evi <- evi +
        annotate("segment",x=k$kstar[i],xend=k$kstar[i],y=evi.at.kstar[i],yend=-Inf,linetype=2,colour="grey50") +
        annotate("segment",x=k$kstar[i],xend=-Inf,y=evi.at.kstar[i],yend=evi.at.kstar[i],linetype=2,colour="grey50")
    }
  }
  
  evi <- evi +
    theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank()) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
  
  return(evi)
}
}
##############################################################################################################

##############################################################################################################
## Plots the main health economics outcomes in just one graph
plot.bcea <- function(x,comparison=NULL,wtp=25000,pos=FALSE,graph=c("base","ggplot2"),...) {
	options(scipen=10)
	base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
  
  if(base.graphics) {
	op <- par(mfrow=c(2,2))
	ceplane.plot(x,comparison=comparison,wtp=wtp,pos=pos,graph="base")
	eib.plot(x,comparison=comparison,pos=pos,graph="base")
	ceac.plot(x,comparison=comparison,pos=pos,graph="base")
	evi.plot(x,graph="base")
  par(op)
  }
  else{
    
    if(!require(ggplot2) & !require(grid)){
      message("falling back to base graphics\n")
      plot.bcea(x,comparison=comparison,wtp=wtp,pos=pos,graph="base",...)
      return(invisible(NULL))
    }
    
    ####### multiplot ###### 
    # source: R graphics cookbook
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      require(grid)
      plots <- c(list(...),plotlist)
      numPlots = length(plots)
      if(is.null(layout)) {
        layout <- matrix(seq(1,cols*ceiling(numPlots/cols)),
                         ncol=cols, nrow=ceiling(numPlots/cols))
      }
      if(numPlots==1) {
        print(plots[[1]])
      } else {
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(nrow(layout),ncol(layout))))
        
        for(i in 1:numPlots) {
          matchidx <- as.data.frame(which(layout==i,arr.ind=TRUE))
          print(plots[[i]],vp=viewport(layout.pos.row=matchidx$row,
                                       layout.pos.col=matchidx$col))
        }
      }
    } #### multiplot end ####
    
    k <- x
    rm(x)
    
    theme.multiplot <- theme(text=element_text(size=9),legend.key.size=unit(.5,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),plot.title=element_text(lineheight=1,face="bold",size=11.5))
    
    exArgs <- list(...)
    for(obj in exArgs)
      if(is.theme(obj))
        theme.multiplot <- theme.multiplot + obj
    
    ceplane.pos <- pos
    if(isTRUE(pos==FALSE)){
      ceplane.pos <- c(1,1.025)
    }
    ceplane <- ceplane.plot(k,wtp=wtp,pos=ceplane.pos,comparison=comparison,graph="ggplot2",...) +
      theme.multiplot
    eib <- eib.plot(k,pos=pos,comparison=comparison,graph="ggplot2",...) +
      theme.multiplot
    ceac <- ceac.plot(k,pos=pos,comparison=comparison,graph="ggplot2") +
      theme.multiplot
    evi <- evi.plot(k,graph="ggplot2") +
      theme.multiplot
    # then call multiplot
    multiplot(ceplane,ceac,eib,evi,cols=2)
  } # !base.graphics
}
##############################################################################################################



##############################################################################################################
## Contour plots for the cost-effectiveness plane
contour.bcea <- function(x,comparison=1,scale=0.5,nlevels=4,levels=NULL,pos=c(1,0),graph=c("base","ggplot2"),...) {
require(MASS)
options(scipen=10)
# comparison selects which plot should be made
# by default it is the first possible

alt.legend <- pos
base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 

if(base.graphics){
  
  if(is.null(comparison))
    comparison <- 1
  
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
if(any(is.na(density$z))==FALSE) {
	if (is.null(levels)==FALSE){
		# Normalise the density and use levels in the contour
		density$z <- (density$z-min(density$z))/(max(density$z)-min(density$z))
		contour(density$x,density$y,density$z,add=TRUE,levels=levels,drawlabels=TRUE)
	}
	if (is.null(levels)==TRUE) {
		contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE)
	}
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
if(any(is.na(density$z))==FALSE) {
	contour(density$x,density$y,density$z,add=TRUE,drawlabels=TRUE)
	if (is.null(levels)==FALSE){
		# Normalise the density and use levels in the contour
		density$z <- (density$z-min(density$z))/(max(density$z)-min(density$z))
		contour(density$x,density$y,density$z,add=TRUE,levels=levels,drawlabels=TRUE)
	}
	if (is.null(levels)==TRUE) {
		contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE)
	}
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
} # if base.graphics
else{
  if(!isTRUE(require(ggplot2)&require(grid))){
    message("falling back to base graphics\n")
    contour.bcea(x,comparison=comparison,scale=scale,nlevels=nlevels,pos=alt.legend,levels=levels,graph="base",...)
    return(invisible(NULL))
  }
  
  if(!is.null(levels))
    message("option level will be ignored using ggplot2 graphics")

  # no visible binding note
  delta.e <- delta.c <- e <- z <- y <- hjust <- label <- NULL
  
  k <- x
  rm(x)
  
  if(!is.null(nlevels)){
    nlevels <- round(nlevels)
    if(nlevels<0)
      nlevels <- 10
    if(nlevels==0)
      nlevels <- 1
  }
  
  
  if(k$n.comparisons==1) {
    kd <- data.frame(e=k$delta.e,c=k$delta.c)
    names(kd) <- c("e","c")
    # for scale_x_continuous(oob=)
    do.nothing=function(x,limits) return(x)
    # plot limits
    range.e <- range(kd$e)
    range.c <- range(kd$c)
    range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
    range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
    
    # labels
    p.ne <- sum(k$delta.e > 0 & k$delta.c > 0)/k$n.sim
    p.ne <- paste0("Pr(Delta[e]>0, Delta[c]>0)==",format(p.ne,digits=4,nsmall=3))
    p.nw <- sum(k$delta.e <= 0 & k$delta.c > 0)/k$n.sim
    p.nw <- paste0("Pr(Delta[e]<=0, Delta[c]>0)==",format(p.nw,digits=4,nsmall=3))
    p.sw <- sum(k$delta.e <= 0 & k$delta.c <= 0)/k$n.sim
    p.sw <- paste0("Pr(Delta[e]<=0, Delta[c]<=0)==",format(p.sw,digits=4,nsmall=3))
    p.se <- sum(k$delta.e > 0 & k$delta.c <= 0)/k$n.sim
    p.se <- paste0("Pr(Delta[e]>0, Delta[c]<=0)==",format(p.se,digits=4,nsmall=3))
    
    # labels dataframe
    labels.df <- data.frame(x=c(range.e[2],range.e[1],range.e[1],range.e[2]),y=c(rep(range.c[2],2),rep(range.c[1],2)),label=c(p.ne,p.nw,p.sw,p.se),hjust=as.factor(c(1,0,0,1)))
    
    # actual plot
    points.colour="grey"
    if(nlevels==1)
      points.colour="black"
    ceplane <- ggplot(kd, aes(e,c)) + geom_hline(aes(yintercept=0),colour="grey") + geom_vline(aes(xintercept=0),colour="grey") +
      theme_bw() + geom_point(size=1,color=points.colour) +
      scale_x_continuous(limits=range.e,oob=do.nothing) + scale_y_continuous(limits=range.c,oob=do.nothing)
    if(!is.null(scale)&require(MASS)){
      density <- kde2d(k$delta.e,k$delta.c,n=300,h=c(sd(k$delta.e)/scale,sd(k$delta.c)/scale))
      densitydf <- data.frame(expand.grid(e=density$x,c=density$y),z=as.vector(density$z))
      ceplane <- ceplane + geom_contour(aes(z=z),data=densitydf,colour="black",bins=nlevels)
    }
    else{
      ceplane <- ceplane + stat_density2d(color="black")
    }
    
    
    ceplane <- ceplane + 
      geom_text(data=labels.df,aes(x=x,y=y,hjust=hjust,label=label),parse=TRUE,size=rel(3.5))
  }
  if(k$n.comparisons>1&is.null(comparison)==TRUE) {
    # creates dataframe for plotting
    kd <- data.frame(c(k$delta.e),c(k$delta.c))
    names(kd) <- c("delta.e","delta.c")
    kd$comparison <- as.factor(sort(rep(1:k$n.comparisons,dim(k$delta.e)[1])))
    
    # vector of values for color, take out white, get integer values
    colors.label <- paste0("gray",round(seq(0,100,length.out=(k$n.comparisons+1))[-(k$n.comparisons+1)]))
    comparisons.label <- paste0(k$interventions[k$ref]," vs ",k$interventions[k$comp])
    do.nothing=function(x,limits) return(x)
    # plot limits
    range.e <- range(kd$delta.e)
    range.c <- range(kd$delta.c)
    range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
    range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
    
    ceplane <-
      ggplot(kd,aes(x=delta.e,y=delta.c,col=comparison)) +
      geom_hline(yintercept=0,colour="grey") + geom_vline(xintercept=0,colour="grey") + theme_bw() +
      geom_point(size=1) +
      scale_color_manual(label=comparisons.label,values=colors.label,na.value="black") +
      scale_x_continuous(limits=range.e,oob=do.nothing) +
      scale_y_continuous(limits=range.c,oob=do.nothing)
    
    if(!is.null(scale)&require(MASS)) {
      require(MASS)
      densitydf <- data.frame()
      for(i in 1:k$n.comparison) {
        temp <- kde2d(k$delta.e[,i],k$delta.c[,i],n=300,h=c(sd(k$delta.e[,i])/scale,sd(k$delta.c[,i])/scale))
        temp <- data.frame(expand.grid(e=temp$x,c=temp$y),z=as.vector(temp$z))
        densitydf <- rbind(densitydf,cbind(temp,rep(i,dim(temp)[[1]])))
      }
      names(densitydf) <- c("delta.e","delta.c","z","comparison")
      densitydf$comparison <- as.factor(densitydf$comparison)
      ceplane <- ceplane + geom_contour(aes(z=z,colour=comparison),data=densitydf,bins=nlevels) +
        guides(colour=guide_legend(override.aes=list(linetype=0)))
    }
    else{
      ceplane <- ceplane + stat_density2d() + 
        guides(colour=guide_legend(override.aes=list(linetype=0)))
    }
  }
  
  if(k$n.comparisons>1&is.null(comparison)==FALSE) {
    # adjusts bcea object for the correct number of dimensions and comparators
    k$comp <- k$comp[comparison]
    k$delta.e <- k$delta.e[,comparison]
    k$delta.c <- k$delta.c[,comparison]
    k$n.comparators=length(comparison)+1
    k$n.comparisons=length(comparison)
    k$interventions=k$interventions[sort(c(k$ref,k$comp))]
    k$ICER=k$ICER[comparison]
    k$ib=k$ib[,,comparison]
    k$eib=k$eib[,comparison]
    k$U=k$U[,,sort(c(k$ref,comparison+1))]
    k$ceac=k$ceac[,comparison]
    k$ref=rank(c(k$ref,k$comp))[1]
    k$comp=rank(c(k$ref,k$comp))[-1]
    k$mod <- TRUE #
    
    return(contour.bcea(k,scale=scale,pos=alt.legend,nlevels=nlevels,graph="ggplot2",comparison=NULL))
  }
  
  labs.title <- "Cost-Effectiveness Plane"
  labs.title <- paste0(labs.title,
                       ifelse(k$n.comparisons==1,
                              paste0("\n",k$interventions[k$ref]," vs ",k$interventions[-k$ref]),
                              paste0(
                                ifelse(isTRUE(k$mod),
                                       paste0("\n",k$interventions[k$ref]," vs ",
                                              paste0(k$interventions[k$comp],collapse=", ")),
                                       ""))))
  
  ceplane <- ceplane + labs(title=labs.title,x="Effectiveness differential",y="Cost differential")
  
  jus <- NULL
  if(isTRUE(alt.legend)) {
    alt.legend="bottom"
    ceplane <- ceplane + theme(legend.direction="vertical")
  }
  else{
    if(is.character(alt.legend)) {
      choices <- c("left", "right", "bottom", "top")
      alt.legend <- choices[pmatch(alt.legend,choices)]
      jus="center"
      if(is.na(alt.legend))
        alt.legend=FALSE
    }
    if(length(alt.legend)>1)
      jus <- alt.legend
    if(length(alt.legend)==1 & !is.character(alt.legend)) {
      alt.legend <- c(1,0)
      jus <- alt.legend
    }
  }
  
  ceplane <- ceplane + 
    theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank()) +
    theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),legend.text.align=0) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
  
  return(ceplane)
} # !base.graphics
}

##############################################################################################################
contour2 <- function(he,wtp=25000,xl=NULL,yl=NULL,comparison=NULL,graph=c("base","ggplot2"),...) {
# Forces R to avoid scientific format for graphs labels
options(scipen=10)

base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 

if(base.graphics) {
# Encodes characters so that the graph can be saved as ps or pdf
ps.options(encoding="CP1250")
pdf.options(encoding="CP1250")

# Selects the first comparison by default if not selected
if(is.null(comparison))
  comparison <- 1

## Selects the correct values for the case of multiple comparisons
if (he$n.comparisons>1) {
	he$delta.e <- he$delta.e[,comparison]
	he$delta.c <- he$delta.c[,comparison]
	he$comp <- he$comp[comparison]
	he$ICER <- he$ICER[comparison]
}

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
return(invisible(NULL))
} # end if base.graphics
else{
  
  if(!isTRUE(require(ggplot2)&require(grid))){
    message("falling back to base graphics\n")
    contour2(he,comparison=comparison,xl=xl,yl=yl,wtp=wtp,graph="base"); return(invisible(NULL))
  }
  
  # no visible binding note
  z <- e <- NULL
  
  k <- he
  rm(he)
  scale=0.5
  nlevels=5
  
  require(MASS)
  
  if(k$n.comparisons==1) {
  density <- kde2d(k$delta.e,k$delta.c,n=300,h=c(sd(k$delta.e)/scale,sd(k$delta.c)/scale))
  densitydf <- data.frame(expand.grid(e=density$x,c=density$y),z=as.vector(density$z))
  contour <- ceplane.plot(k,wtp=wtp,graph="ggplot2",...) +
    geom_contour(aes(z=z,x=e,y=c),data=densitydf,colour="black",bins=nlevels)
  }
  if(k$n.comparisons>1&is.null(comparison)) {
    densitydf <- data.frame()
    for(i in 1:k$n.comparisons) {
      density <- kde2d(k$delta.e[,i],k$delta.c[,i],n=300,h=c(sd(k$delta.e)/scale,sd(k$delta.c)/scale))
      densitydf <- rbind(densitydf,cbind(expand.grid(density$x,density$y),as.vector(density$z)))
    }
    names(densitydf) <- c("e","c","z")
    densitydf <- cbind(densitydf,"comparison"=as.factor(sort(rep(1:k$n.comparisons,dim(densitydf)[1]/k$n.comparisons))))
    contour <- ceplane.plot(k,wtp=wtp,graph="ggplot2",...) +
      geom_contour(data=densitydf,aes(x=e,y=c,z=z,colour=comparison),bins=nlevels,linetype=1)
  }
  if(k$n.comparisons>1&!is.null(comparison)) {
    # adjusts bcea object for the correct number of dimensions and comparators
    k$comp <- k$comp[comparison]
    k$delta.e <- k$delta.e[,comparison]
    k$delta.c <- k$delta.c[,comparison]
    k$n.comparators=length(comparison)+1
    k$n.comparisons=length(comparison)
    k$interventions=k$interventions[sort(c(k$ref,k$comp))]
    k$ICER=k$ICER[comparison]
    k$ib=k$ib[,,comparison]
    k$eib=k$eib[,comparison]
    k$U=k$U[,,sort(c(k$ref,comparison+1))]
    k$ceac=k$ceac[,comparison]
    k$ref=rank(c(k$ref,k$comp))[1]
    k$comp=rank(c(k$ref,k$comp))[-1]
    k$mod <- TRUE #
    
    return(contour2(k,wtp=wtp,xl=xl,yl=yl,comparison=NULL,graph="ggplot2",...))
  }
    
  contour <- contour + coord_cartesian(xlim=xl,ylim=yl)
  return(contour)
} # end if !base.graphics

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
plot.CEriskav <- function(x,pos=c(0,1),graph=c("base","ggplot2"),...) {
options(scipen=10)

alt.legend <- pos
base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 

howplot <- NULL
exArgs <- list(...)
if(length(exArgs)>=1)
  if(exists("plot",where=exArgs))
    howplot <- exArgs$plot

if(base.graphics) {
  
  if(is.numeric(alt.legend)&length(alt.legend)==2){
    temp <- ""
    if(alt.legend[2]==0)
      temp <- paste0(temp,"bottom")
    else
      temp <- paste0(temp,"top")
    if(alt.legend[1]==1)
      temp <- paste0(temp,"right")
    else
      temp <- paste0(temp,"left")
    alt.legend <- temp
    if(length(grep("^(bottom|top)(left|right)$",temp))==0)
      alt.legend <- FALSE
  }
  if(is.logical(alt.legend)){
    if(!alt.legend)
      alt.legend="topright"
    else
      alt.legend="topleft"
  }
  
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
legend(alt.legend,text,lty=seq(1:x$R),cex=.9,box.lty=0)
abline(h=0,col="grey")

# Plots the EVPI for the risk aversion case
  if(is.null(howplot)) {
    if(!isTRUE(Sys.getenv("RSTUDIO")==1))
      dev.new()
  }
  else{
    opt <- c("x11","ask","dev.new")
    howplot <- ifelse(is.na(pmatch(howplot,opt)),"dev.new",opt[pmatch(howplot,opt)])
    if(howplot=="x11")
      dev.new()
    if(howplot=="dev.new")
      dev.new()
    if(howplot=="ask")
      devAskNewPage(ask=TRUE)
  }
  
plot(x$k,x$evir[,1],t="l",ylim=range(x$evir),xlab="Willingness to pay",ylab=" ",main="EVPI as a function of the risk aversion parameter")
for (l in 2:x$R) {
	points(x$k,x$evir[,l],t="l",lty=linetype[l])
}
legend(alt.legend,text,lty=seq(1:x$R),cex=.9,box.lty=0)
abline(h=0,col="grey")
  
  if(!is.null(howplot))
    if(howplot=="ask")
      devAskNewPage(ask=FALSE)
  
} # base.graphics
else{
  if(!isTRUE(require(ggplot2)&require(grid))) {
    message("falling back to base graphics\n")
    plot.CEriskav(x,graph="base",pos=pos,...)
    return(invisible(NULL))
  }
  # no visible bindings note
  k <- r <- NULL
  
  linetypes <- rep(c(1,2,3,4,5,6),ceiling(x$R/6))[1:x$R]
  df <- data.frame(cbind(rep(x$k,x$R),c(x$eibr),c(x$evir)),as.factor(sort(rep(1:x$R,length(x$k)))))
  names(df) <- c("k","eibr","evir","r")
  
  # labels
  text <- paste0("r = ",x$r)
  # if the first value for r is small enough, consider it close to 0 and print the label accordingly
  if(x$r[1]<1e-8) {
    text[1] <- expression(r%->%0)
  }
  
  eibr <- ggplot(df,aes(x=k,y=eibr,linetype=r))+geom_hline(yintercept=0,linetype=1,colour="grey50")+geom_line()+scale_linetype_manual("",labels=text,values=linetypes)+theme_bw() +
    labs(title="EIB as a function of the risk aversion parameter",x="Willingness to pay",y="EIB") +
    theme(text=element_text(size=11),legend.key.size=unit(.66,"line"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank())
  
  ### evir ###
  evir <- ggplot(df,aes(x=k,y=evir,linetype=r))+geom_hline(yintercept=0,linetype=1,colour="grey50")+geom_line()+scale_linetype_manual("",labels=text,values=linetypes)+theme_bw() +
    labs(title="EVPI as a function of the risk aversion parameter",x="Willingness to pay",y="EVPI") +
    theme(text=element_text(size=11),legend.key.size=unit(.66,"line"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank())
    
  jus <- NULL
  if(isTRUE(alt.legend)) {
    alt.legend="bottom"
    eibr <- eibr + theme(legend.direction="vertical")
    evir <- evir + theme(legend.direction="vertical")
  }
  else{
    if(is.character(alt.legend)) {
      choices <- c("left", "right", "bottom", "top")
      alt.legend <- choices[pmatch(alt.legend,choices)]
      jus="center"
      if(is.na(alt.legend))
        alt.legend=FALSE
    }
    if(length(alt.legend)>1)
      jus <- alt.legend
    if(length(alt.legend)==1 & !is.character(alt.legend)) {
      alt.legend <- c(0,1)
      jus <- alt.legend
    }
  }
  
  eibr <- eibr + 
    theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(face="bold"),legend.text.align=0) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
  
  evir <- evir + 
    theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(face="bold"),legend.text.align=0) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
  
  plot(eibr)
  
  if(is.null(howplot)) {
    if(!isTRUE(Sys.getenv("RSTUDIO")==1))
      dev.new()
  }
  else{
    opt <- c("x11","ask","dev.new")
    howplot <- ifelse(is.na(pmatch(howplot,opt)),"dev.new",opt[pmatch(howplot,opt)])
    if(howplot=="x11")
      dev.new()
    if(howplot=="dev.new")
      dev.new()
    if(howplot=="ask")
      devAskNewPage(ask=TRUE)
  }
  
  plot(evir)
  
  if(!is.null(howplot))
    if(howplot=="ask")
      devAskNewPage(ask=FALSE)
  
  return(invisible(list("eib"=eibr,"evi"=evir)))
}

}


##############################################################################################################
mixedAn <- function(he,mkt.shares=NULL,plot=FALSE) UseMethod("mixedAn")


##############################################################################################################
mixedAn.default <- function(he,mkt.shares=NULL,plot=FALSE) {
ma <- mixedAn.fn(he,mkt.shares)
class(ma) <- "mixedAn"
if(plot) {
  plot.mixedAn(ma)
}
return(ma)
}


##############################################################################################################
mixedAn.fn <- function(he,mkt.shares=NULL,plot=FALSE) {
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

## Outputs of the function
return(list(
Ubar=Ubar,OL.star=OL.star,evi.star=evi.star,k=he$k,Kmax=he$Kmax,step=he$step,
ref=he$ref,comp=he$comp,mkt.shares=mkt.shares,n.comparisons=he$n.comparisons,
interventions=he$interventions,evi=he$evi
))
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
plot.mixedAn <- function(x,y.limits=NULL,pos=c(0,1),graph=c("base","ggplot2"),...) {
## Plot the EVPI and the mixed strategy
options(scipen=10)

alt.legend <- pos
base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 

if(is.null(y.limits)){
y.limits=range(x$evi,x$evi.star)
}

if(base.graphics) {
  
  if(is.numeric(alt.legend)&length(alt.legend)==2){
    temp <- ""
    if(alt.legend[2]==0)
      temp <- paste0(temp,"bottom")
    else
      temp <- paste0(temp,"top")
    if(alt.legend[1]==1)
      temp <- paste0(temp,"right")
    else
      temp <- paste0(temp,"left")
    alt.legend <- temp
    if(length(grep("^(bottom|top)(left|right)$",temp))==0)
      alt.legend <- FALSE
  }
  if(is.logical(alt.legend)){
    if(!alt.legend)
      alt.legend="topright"
    else
      alt.legend="topleft"
  }
  
plot(x$k,x$evi,t="l",xlab="Willingness to pay",ylab="EVPI",
	main="Expected Value of Information",ylim=y.limits)
polygon(c(x$k,rev(x$k)),c(x$evi.star,rev(x$evi)),density=20,col="grey")
points(x$k,x$evi.star,t="l",col="red")
points(x$k,x$evi,t="l",col="black")
txt <- c("Optimal strategy","Mixed strategy:",
	paste("   ",x$interventions,"=",format(100*x$mkt.shares,digits=3,nsmall=2),"%",sep=""))
cols <- c("black","red",rep("white",length(x$interventions)))
legend(alt.legend,txt,col=cols,cex=.6,bty="n",lty=1)
} # base.graphics
else {
  if(!isTRUE(require(ggplot2)&require(grid))) {
    message("falling back to base graphics\n")
    plot.mixedAn(x,y.limits=y.limits,pos=pos,graph="base")
    return(invisible(NULL))
  } 
  # no visible bindings note
  k <- evi.star <- NULL
  
  # legend
  txt <- c("Optimal strategy",paste0("Mixed strategy:", paste0("\n   ",x$interventions,"=",format(100*x$mkt.shares,digits=3,nsmall=2),"%",collapse="")))
  colours=c("black","red")
  
  df <- data.frame(k=x$k,evi=x$evi,evi.star=x$evi.star)
  evi <- ggplot(df) + theme_bw() +
    geom_ribbon(aes(x=k,ymin=evi,ymax=evi.star),colour="lightgrey",alpha=.2) +
    geom_line(aes(x=k,y=evi,colour=as.factor(1))) +
    geom_line(aes(x=k,y=evi.star,colour=as.factor(2))) +
    coord_cartesian(ylim=y.limits) +
    scale_colour_manual("",labels=txt,values=colours) +
    labs(title="Expected Value of Information",x="Willingness to pay",y="EVPI") +
    theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),plot.title=element_text(face="bold"))
  
  jus <- NULL
  if(isTRUE(alt.legend)) {
    alt.legend="bottom"
    evi <- evi + theme(legend.direction="vertical")
  }
  else{
    if(is.character(alt.legend)) {
      choices <- c("left", "right", "bottom", "top")
      alt.legend <- choices[pmatch(alt.legend,choices)]
      jus="center"
      if(is.na(alt.legend))
        alt.legend=FALSE
    }
    if(length(alt.legend)>1)
      jus <- alt.legend
    if(length(alt.legend)==1 & !is.character(alt.legend)) {
      alt.legend <- c(0,1)
      jus <- alt.legend
    }
  }
  
  evi <- evi + 
    theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(face="bold"),legend.text.align=0) +
    theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
  
  return(evi)
}
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

############mce.plot###################################
mce.plot <- function(mce,pos=c(1,0.5),graph=c("base","ggplot2")){
  
  alt.legend <- pos
  base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
  
if(base.graphics) {
  
  if(is.numeric(alt.legend)&length(alt.legend)==2){
    temp <- ""
    if(alt.legend[2]==0)
      temp <- paste0(temp,"bottom")
    else if(alt.legend[2]!=0.5)
      temp <- paste0(temp,"top")
    if(alt.legend[1]==1)
      temp <- paste0(temp,"right")
    else
      temp <- paste0(temp,"left")
    alt.legend <- temp
    if(length(grep("^((bottom|top)(left|right)|right)$",temp))==0)
      alt.legend <- FALSE
  }
  if(is.logical(alt.legend)){
    if(!alt.legend)
      alt.legend="topright"
    else
      alt.legend="right"
  }
  
	color <- rep(1,(mce$n.comparators+1)); lwd <- 1
	if (mce$n.comparators>7) {
		cl <- colors()
		color <- cl[floor(seq(262,340,length.out=mce$n.comparators))]	# gray scale
		lwd <- 1.5
	}

	plot(mce$k,mce$m.ce[,1],t="l",col=color[1],lwd=lwd,lty=1,xlab="Willingness to pay",
		ylab="Probability of most cost effectiveness",ylim=c(0,1),
		main="Cost-effectiveness acceptability curve \nfor multiple comparisons")
	for (i in 2:mce$n.comparators) {
		points(mce$k,mce$m.ce[,i],t="l",col=color[i],lwd=lwd,lty=i)
	}
	legend(alt.legend,mce$interventions,col=color,cex=.7,bty="n",lty=1:mce$n.comparators)
} # base graphics
  else{
    if(!isTRUE(require(ggplot2)&require(grid))) {
      message("falling back to base graphics\n")
      mce.plot(mce,pos=pos,graph="base")
      return(invisible(NULL))
    }
    # no visible bindings note
    k <- ce <- comp <- ceplane <- NULL
    
    alt.legend <- pos
    lty <- rep(1:6,ceiling(mce$n.comparators/6))[1:mce$n.comparators]
    col <- paste0("gray",round(seq(0,100,length.out=mce$n.comparators+1))[-(mce$n.comparators+1)])
    label <- paste0(mce$interventions)
    
    df <- cbind(rep(mce$k,mce$n.comparators),c(mce$m.ce))
    df <- data.frame(df,as.factor(sort(rep(1:mce$n.comparators,length(mce$k)))))
    names(df) <- c("k","ce","comp")
    
    mceplot <- ggplot(df) + theme_bw() +
      geom_line(aes(x=k,y=ce,colour=comp,linetype=comp)) + 
      scale_colour_manual("",labels=label,values=col) +
      scale_linetype_manual("",labels=label,values=lty) +
      labs(title="Cost-effectiveness acceptability curve\nfor multiple comparisons",x="Willingness to pay",y="Probabiity of most cost effectiveness") +
      theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),plot.title=element_text(face="bold"))
    
    jus <- NULL
    if(isTRUE(alt.legend)) {
      alt.legend="bottom"
      mceplot <- ceplane + theme(legend.direction="vertical")
    }
    else{
      if(is.character(alt.legend)) {
        choices <- c("left", "right", "bottom", "top")
        alt.legend <- choices[pmatch(alt.legend,choices)]
        jus="center"
        if(is.na(alt.legend))
          alt.legend=FALSE
      }
      if(length(alt.legend)>1)
        jus <- alt.legend
      if(length(alt.legend)==1 & !is.character(alt.legend)) {
        alt.legend <- c(1,0.5)
        jus <- alt.legend
      }
    }
    
    mceplot <- mceplot + coord_cartesian(ylim=c(0,1)) +
      theme(legend.position=alt.legend,legend.justification=jus,legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(face="bold"),legend.text.align=0) +
      theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
      
    return(mceplot)
  }
}

#######################ceaf.plot##################################
ceaf.plot <- function(mce,graph=c("base","ggplot2")){
  base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
  if(base.graphics) {
	plot(mce$k,mce$ceaf,t="l",lty=1,
		ylim=c(0,1),xlab="Willingness to pay",
		ylab="Probability of most cost effectiveness",
		main="Cost-effectiveness acceptability frontier")
  }
  else{
    if(!isTRUE(require(ggplot2)&require(grid))){
      message("falling back to base graphics\n")
      ceaf.plot(mce,graph="base")
      return(invisible(NULL))
    }
    # no visible binding note
    k <- NULL
    
    df <- data.frame(k=mce$k,ceaf=mce$ceaf)
    ceaf <- ggplot(df) + theme_bw() +
      geom_line(aes(x=k,y=ceaf)) +
      coord_cartesian(ylim=c(0,1)) +
      theme(text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank(),plot.title=element_text(face="bold")) +
      labs(title="Cost-effectiveness acceptability frontier",x="Willingness to pay",y="Probability of most cost effectiveness") +
      theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
    
    return(ceaf)
  }
}

######evppi###################################################################################################
evppi <- function(param,he) {
  ## Computes the Expected Value of Partial Information with respect to a given parameter
  ## using the algorithm for fast computation of Sadatsafavi et al (2013)
  ## source: http://webservices.core.ubc.ca/research/other/1evppi/
  f <- numeric()  	# will use to store the values of EVPI for each k 

  nb <- he$U		# extracts the net benefits from the full economic model
#  tic()			# useful to check the time used to run this computation 
  evppi.graphics=FALSE	# avoids graphs being printed while running the function
  for (k in 1:length(he$k)) {
    nbs <- nb[,k,]	# this is a n.sims x n.interventions matrix with the NB for each value of k
    f[k] <- evppi.func(param,nbs)$evppi	# EVPPI for each value of k
  }
  
  ### time measurement disabled
#  time2run2 <- toc()				
#  time2run2 <- time2run2[3]	# This gives the time needed to run (in seconds)
  
  # Outputs of the function
  list(evppi=f)#,mins=time2run/60)
}

######evppi.func##############################################################################################
evppi.func <- function(param,nbs,...) {
  evppi.graphics <- FALSE
  exArgs <- list(...)
  if(length(exArgs)>0)
    if(exists("evppi.graphics",where=exArgs))
      evppi.graphics=exArgs$evppi.graphics
  
  ## Computes the Expected Value of Partial Information with respect to a given parameter
  ## using the algorithm for fast computation of Sadatsafavi et al (2013)
  ## source: http://webservices.core.ubc.ca/research/other/1evppi/
  # param = the parameter of interest (must be in the format of a n.sim x 1 simulations from its posterior
  # nbs = a matrix of net benefits for the interventions being compared 
  if(length(unique(param))==1)
    stop("Error: input parameter is not stochastic")
  n<-length(param)
  o<-order(param)
  param<-param[o]
  d<-dim(nbs)[2]
  if(is.vector(nbs))
  {
    nbs<-cbind(nbs[o],0)
    d<-2
  }
  if(d==1)
  {
    nbs<-cbind(nbs[o,],0)
    d<-2
  }
  else
    nbs<-nbs[o,]
  nSegs<-matrix(1,d,d)
  exArgs<-list(...)
  for(obj in exArgs)
  {
    if(is.null(names(obj)))
      if (length(obj)==1)
      {nSegs[1,2]<-obj; nSegs[2,1]<-obj}
    else
    {nSegs[obj[1],obj[2]]<-obj[3]; nSegs[obj[2],obj[1]]<-obj[3];}
  }
  segPoints<-c()
  for(i in 1:(d-1))
    for(j in (i+1):d)
    {
      ###			message(paste('Fitting ',nSegs[i,j],' segmentation points for decisions ',i,j))
      cm<-cumsum(nbs[,i]-nbs[,j])/n
      if(nSegs[i,j]==1)
      {
        l<-which.min(cm)
        u<-which.max(cm)
        if(cm[u]-max(cm[1],cm[n])>min(cm[1],cm[n])-cm[l])
          segPoint<-u
        else
          segPoint<-l
        if (segPoint>1 && segPoint<n)
          segPoints<-c(segPoints, segPoint)
      }
      if(nSegs[i,j]==2)
      {
        distMaxMin<-0
        distMinMax<-0
        minL<-Inf
        maxL<--Inf
        for(k in 1:n)
        {
          #max-min pattern
          if(cm[k]>maxL)
          {
            maxLP<-k
            maxL<-cm[k]
          }
          else
            if(maxL-cm[k]>distMaxMin)
            {
              distMaxMin<-maxL-cm[k]
              segMaxMinL<-maxLP	
              segMaxMinR<-k
            }	
          #min-max pattern
          if(cm[k]<minL)
          {
            minLP<-k
            minL<-cm[k]
          }
          else
            if(cm[k]-minL>distMinMax)
            {
              distMinMax<-cm[k]-minL
              segMinMaxL<-minLP	
              segMinMaxR<-k
            }	
        }
        siMaxMin<-cm[segMaxMinL]+distMaxMin+(cm[n]-cm[segMaxMinR])
        siMinMax<--cm[segMaxMinL]+distMinMax-(cm[n]-cm[segMinMaxR])
        if(siMaxMin>siMinMax)
        {
          segPoint<-c(segMaxMinL,segMaxMinR)
        }
        else
        {
          segPoint<-c(segMinMaxL,segMinMaxR)
        }
        if (segPoint[1]>1 && segPoint[1]<n)
          segPoints<-c(segPoints, segPoint[1])
        if (segPoint[2]>1 && segPoint[2]<n)
          segPoints<-c(segPoints, segPoint[2])
      }
      if(exists('evppi.graphics') && evppi.graphics==TRUE){
        dev.new()
        plot(param,cm)
        for(k in 1:length(segPoint))
        {
          if (segPoint[k]>1 && segPoint[k]<n)
            lines(c(param[segPoint[k]],param[segPoint[k]]),c(min(cm),max(cm)))
        }
        title(paste("Decision",i,"vs.",j))
      }
      else
      {
        evppi.graphics = FALSE
        ###  			message("Graphic output is disabled. Set evppi.graphics=TRUE to enable graphics")
      }
    }
  if(length(segPoints)>0)
  {
    segPoints2<-unique(c(0, segPoints[order(segPoints)], n))
    evppi<-0
    for(j in 1:(length(segPoints2)-1))
      evppi<-evppi+max(colSums(matrix(nbs[(1+segPoints2[j]):segPoints2[j+1],],ncol=d)))/n    
    evppi<-evppi-max(colMeans(nbs))
    return(list(evppi=evppi,segPoints=param[segPoints2[-c(1,length(segPoints2))]]))
  }
  else
    evppi<-0
  return(list(evppi=evppi,segPoints=c()))     
}
######evppi.plot##############################################################################################

evppi.plot <- function(he,evppi=NULL,param=NULL,graph=c("base","ggplot2")) {
  options(scipen=10)
  base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
  
  stopifnot(isTRUE(class(he)=="bcea"))
  stopifnot(!isTRUE(is.null(evppi)&is.null(param)))
    if(!is.null(evppi)){
      if(!is.null(param))
        warning("Both evppi and param supplied. Only evppi will be used.")
      values <- evppi
    }
    else stopifnot(!is.null(param))
      values <- unlist(evppi(param=param,he=he))
  
  if(base.graphics) {
    plot(x=he$k,y=values,t="l",xlab="Willingness to pay",ylab="EVPPI",main="Expected Value of Perfect Partial Information")
    return(invisible(NULL))
  }
  else{
    if(!isTRUE(require(ggplot2)&require(grid))) {
      message("falling back to base graphics\n")
      evppi.plot(he,param=param,graph="base")
      return(invisible(NULL))
    }
    # no visible binding note
    x <- y <- NULL
    
    evppi <- ggplot(data.frame(x=he$k,y=values)) + geom_line(aes(x=x,y=y)) +
      labs(x="Willingness to pay",y="EVPPI",title="Expected Value of Perfect Partial Information") +
      theme_bw() + theme(plot.title=element_text(face="bold"),text=element_text(size=11),legend.key.size=unit(.66,"lines"),legend.margin=unit(-1.25,"line"),panel.grid=element_blank(),legend.key=element_blank()) +
      theme(plot.title = element_text(lineheight=1.05, face="bold",size=14.3))
    return(evppi)
  } 
}
