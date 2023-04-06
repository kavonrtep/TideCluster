#!/usr/bin/env Rscript

## FUNCTIONS:
letterA <- function(x.pos,y.pos,ht,wt,id=NULL){
    
    x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
    y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
    x <- 0.1*x
    y <- 0.1*y
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    if (is.null(id)){
        id <- c(rep(1,9),rep(2,4))
    }else{
        id <- c(rep(id,9),rep(id+1,4))
    }
    
    fill <- c("green","white")
    
    list(x=x,y=y,id=id,fill=fill)
}

## T
letterT <- function(x.pos,y.pos,ht,wt,id=NULL){
    
    x <- c(0,10,10,6,6,4,4,0)
    y <- c(10,10,9,9,0,0,9,9)
    x <- 0.1*x
    y <- 0.1*y
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    if (is.null(id)){
        id <- rep(1,8)
    }else{
        id <- rep(id,8)
    }
    
    fill <- "red"
    
    list(x=x,y=y,id=id,fill=fill)
}

## C
letterC <- function(x.pos,y.pos,ht,wt,id=NULL){
    angle1 <- seq(0.3+pi/2,pi,length=100)
    angle2 <- seq(pi,1.5*pi,length=100)
    x.l1 <- 0.5 + 0.5*sin(angle1)
    y.l1 <- 0.5 + 0.5*cos(angle1)
    x.l2 <- 0.5 + 0.5*sin(angle2)
    y.l2 <- 0.5 + 0.5*cos(angle2)
    
    x.l <- c(x.l1,x.l2)
    y.l <- c(y.l1,y.l2)
    
    x <- c(x.l,rev(x.l))
    y <- c(y.l,1-rev(y.l))
    
    x.i1 <- 0.5 +0.35*sin(angle1)
    y.i1 <- 0.5 +0.35*cos(angle1)
    x.i1 <- x.i1[y.i1<=max(y.l1)]
    y.i1 <- y.i1[y.i1<=max(y.l1)]
    y.i1[1] <- max(y.l1)
    
    x.i2 <- 0.5 +0.35*sin(angle2)
    y.i2 <- 0.5 +0.35*cos(angle2)
    
    x.i <- c(x.i1,x.i2)
    y.i <- c(y.i1,y.i2)
    
    x1 <- c(x.i,rev(x.i))
    y1 <- c(y.i,1-rev(y.i))
    
    x <- c(x,rev(x1))
    y <- c(y,rev(y1))
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    if (is.null(id)){
        id <- rep(1,length(x))
    }else{
        id <- rep(id,length(x))
    }
    
    fill <- "blue"
    
    list(x=x,y=y,id=id,fill=fill)
}


## G
letterG <- function(x.pos,y.pos,ht,wt,id=NULL){
    angle1 <- seq(0.3+pi/2,pi,length=100)
    angle2 <- seq(pi,1.5*pi,length=100)
    x.l1 <- 0.5 + 0.5*sin(angle1)
    y.l1 <- 0.5 + 0.5*cos(angle1)
    x.l2 <- 0.5 + 0.5*sin(angle2)
    y.l2 <- 0.5 + 0.5*cos(angle2)
    
    x.l <- c(x.l1,x.l2)
    y.l <- c(y.l1,y.l2)
    
    x <- c(x.l,rev(x.l))
    y <- c(y.l,1-rev(y.l))
    
    x.i1 <- 0.5 +0.35*sin(angle1)
    y.i1 <- 0.5 +0.35*cos(angle1)
    x.i1 <- x.i1[y.i1<=max(y.l1)]
    y.i1 <- y.i1[y.i1<=max(y.l1)]
    y.i1[1] <- max(y.l1)
    
    x.i2 <- 0.5 +0.35*sin(angle2)
    y.i2 <- 0.5 +0.35*cos(angle2)
    
    x.i <- c(x.i1,x.i2)
    y.i <- c(y.i1,y.i2)
    
    x1 <- c(x.i,rev(x.i))
    y1 <- c(y.i,1-rev(y.i))
    
    x <- c(x,rev(x1))
    y <- c(y,rev(y1))
    
    r1 <- max(x.l1)
    
    h1 <- 0.4
    x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
    y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)
    
    
    
    if (is.null(id)){
        id <- c(rep(1,length(x)),rep(2,length(x.add)))
    }else{
        id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
    }
    
    x <- c(rev(x),x.add)
    y <- c(rev(y),y.add)
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    
    fill <- c("orange","orange")
    
    list(x=x,y=y,id=id,fill=fill)
    
}

Letter <- function(which,x.pos,y.pos,ht,wt){
    
    if (which == "A"){
        letter <- letterA(x.pos,y.pos,ht,wt)
    }else if (which == "C"){
        letter <- letterC(x.pos,y.pos,ht,wt)    
    }else if (which == "G"){
        letter <- letterG(x.pos,y.pos,ht,wt)    
    }else if (which == "T"){
        letter <- letterT(x.pos,y.pos,ht,wt)    
    }else{
        stop("which must be one of A,C,G,T")
    }
    
    letter
}


plot_multiline_logo <- function(cons.logo, read=NULL, W=50, setpar=TRUE, gaps = NULL){
    ## logo - base order  - A C G T
    if (ncol(cons.logo)==5){
        gaps_prob <- cons.logo[, 5]
    }else{
        gaps_prob <- NULL
    }
    tm <- 4
    pwm <- as.matrix(cons.logo[, 1:4])
    N <- nrow(pwm)
    Nori <- N
    if (N<W){
        W <- N
    }
    s1 <- seq(1, N, by=W)
    s2 <- seq(W, N, by=W)
    if (length(s2)<length(s1)){
	pwm <- rbind(pwm, matrix(0, nrow=W*length(s1)-N, ncol=4, dimnames=list(NULL, c('A', 'C', 'G', 'T'))))
        if (!is.null(read)){
            pwm_read <- rbind(read, matrix(0, nrow=W*length(s1)-N, ncol=4, dimnames=list(NULL, c('A', 'C', 'G', 'T'))))
        }
	N <- nrow(pwm)
	s2 <- seq(W, N, by=W)
    }
    if (setpar){
        par(mfrow = c(ceiling(N/W),1), mar=c(1,4,1,0))
    }
    for (i in seq_along(s1)){
        if (!is.null(read)){
            plot.logo(pwm_read[s1[i]:s2[i],],maxh=2)
        }
        plot.logo(pwm[s1[i]:s2[i],],maxh=max(rowSums(cons.logo)))
        if(!is.null(gaps)){
            ## transparent rectangles
            rect((gaps[ ,'start']-s1[i]+1),0, (gaps[,'end']-s1[i]+2), max(pwm), col="#00000005")
            
        }
        if(!is.null(gaps_prob)){
            rect(seq_along(s1[i]:s2[i]),
                 max(rowSums(cons.logo)),
                 seq_along(s1[i]:s2[i])+1,
                 max(rowSums(cons.logo)) - gaps_prob[s1[i]:s2[i]],
                 col="#00000030")

            
        }
        ticks <- intersect(intersect(pretty(pretty(s1[i]:s2[i])+1), s1[i]:s2[i]), 1:Nori)
        axis(1,at=ticks+1.5-s1[i],label=ticks,tick=FALSE)
        y <- pretty(c(0, max(pwm)), n=tm)
        axis(2,at=y,label=y,las=2,cex.axis=.7)
    }
}

plot.logo <- function(pwm, maxh=NULL){
    acgt <- c("A", "C", "G", "T")
    pwm <- pwm[, acgt]
    nbp <- dim(pwm)[1]
    if (is.null(maxh)) {maxh <- max(rowSums(pwm))}
    
    plot(0,0,xlim=c(0,nbp + 1),ylim=c(0,maxh),type="n",axes=F,xlab="",ylab="")
    for ( i in 1:nbp){
        S <- order(pwm[i,])
        hgts <- pwm[i, S]
        nts <- acgt[S]
        ypos <- c(0, cumsum(hgts)[1:3])
        for (j in 1:4){
            if (hgts[j]==0) next
            L <- Letter(which=nts[j], x.pos=i, y.pos=ypos[j], ht=hgts[j], wt=1)
            Id <- L$id==1
            polygon(L$x[Id],L$y[Id],lty=0,col=L$fill[1])
            if (sum(L$id==2)>0) {
                polygon(L$x[!Id],L$y[!Id],lty=0,col=L$fill[2])
            }
        }
    }
}



