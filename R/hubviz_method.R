
# generic methods for "hubviz" class

setMethod(
  f="show",
  signature="hubviz",
  definition=function( object ) {
    
    colname <- colnames(object@data)
    position <- object@result$w.estimate
    
    cat( "Summary: hub-centric visualization (class: hubviz)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Model settings:\n")
    cat( "Number of samples: ", object@init$nsample, "\n", sep="" )
    cat( "Number of variables: ", object@init$nitem, "\n", sep="" )
    cat( "Dimension of latent space: ", object@init$ndim, "\n", sep="" )
    cat( "Number of iterations: ", object@init$niter, "\n", sep="" )
    cat( "Number of burn-in: ", object@init$nburn, "\n", sep="" )
    cat( "Number of thining: ", object@init$nthin, "\n", sep="" )
    cat( "--------------------------------------------------\n" )
    cat( "Model results: Latent position \n")
    for ( i in 1:nrow(position) ) { 
      cat( "   ", paste(colname[i],": ","(", paste( format( round(position[i,], digits=2), nsmall=2 ), collapse="," ),")" ,sep=""), "\n", sep="" )
    }
    cat( "--------------------------------------------------\n" )
  }
)


 
setGeneric(name="estimate",
           def=function(object)
           {
             standardGeneric("estimate")
           }
)

setMethod(
  f="estimate",
  signature="hubviz",
  definition=function( object ) {
    
     #extract objects
    
    result <- object@result
    
    theta <- result$theta.estimate
    w <- result$w.estimate
    sigma <- result$sigma.w.estimate
    
    re1 <- list( theta = theta, w = w, sigma = sigma )
    
    if (length(result)==3) {
    return(re1)
    } else {
    return(result)  
    }  
  }
)


setMethod(
  f="plot",
  signature=c("hubviz", "missing"),
  definition=function( x, y, xlim=NA, ylim=NA ) {

    # extract objects
    colname <- colnames(x@data)
    ndim <- x@init$ndim
    
    position <- as.data.frame(x@result$w.estimate)
    rownames(position) <- colname
    colnames(position) <- paste("position",1:ndim,sep="")
    
    if (any(is.na(xlim))) {
      x1 <- -max(abs(position[,1]))-2
      x2 <- max(abs(position[,1]))+2
    } else {
       x1 <- xlim[1]
       x2 <- xlim[2]
    }
    if (any(is.na(ylim))) {
      y1 <- -max(abs(position[,2]))-2
      y2 <- max(abs(position[,2]))+2  
    } else {
      y1 <- ylim[1]
      y2 <- ylim[2]
    }
    
    # plot
      ggplot(position,aes(x=position1,y=position2))+
      geom_point(color="red")+
      xlim(x1,x2) + ylim(y1,y2) +
	  xlab("Position 1") + ylab("Position 2") +
      geom_hline(yintercept = 0, color = "gray70", linetype=2) +
      geom_vline(xintercept = 0, color = "gray70", linetype=2) + 
      geom_text_repel(aes(y = position2 + 0.25), label=rownames(position), segment.color = "grey50")    

  }
)


setGeneric(name="hpdplot",
           def=function(object, xlim=NA, ylim=NA)
           {
             standardGeneric("hpdplot")
           }
)


setMethod(
  f="hpdplot",
  signature=c("hubviz"),
  definition=function( object, xlim=NA, ylim=NA ) {
    
    # extract objects
    colname <- colnames(object@data)
    ndim <- object@init$ndim
    
    position <- as.data.frame(object@result$w.estimate)
    rownames(position) <- colname
    colnames(position) <- paste("position",1:ndim,sep="")
    
    if (any(is.na(xlim))) {
      x1 <- -max(abs(position[,1]))-2
      x2 <- max(abs(position[,1]))+2
    } else {
      x1 <- xlim[1]
      x2 <- xlim[2]
    }
    if (any(is.na(ylim))) {
      y1 <- -max(abs(position[,2]))-2
      y2 <- max(abs(position[,2]))+2  
    } else {
      y1 <- ylim[1]
      y2 <- ylim[2]
    }
    
    
    range <- object@result$w.hpd
    
    x.min <- apply(abs(t(range[,,1])-position[,1]),1,min)
    y.min <- apply(abs(t(range[,,2])-position[,2]),1,min)
    min.radius <- apply(cbind(x.min,y.min),1,min)
    
    position <- cbind(position,radius=min.radius)
    
    ggplot(position, aes(x = position1, y = position2)) + 
      xlim(x1,x2) + ylim(y1,y2) +
      geom_point(aes(size = radius*2), alpha = 0.5,color='darkblue') +
      theme(legend.title=element_blank())+
      scale_size(range = c(0.5, 12))+
    geom_text_repel(aes(y = position2 + 0.25), label=rownames(position), segment.color = "grey50")    
    
    
    
  }
)

