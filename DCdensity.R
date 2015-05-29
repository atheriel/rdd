#' McCrary Sorting Test
#'
#' \code{DCdensity} implements the McCrary (2008) sorting test.
#'
#' @param runvar numerical vector of the running variable
#' @param cutpoint the cutpoint (defaults to 0)
#' @param bin the binwidth (defaults to \code{2*sd(runvar)*length(runvar)^(-.5)})
#' @param bw the bandwidth to use (by default uses bandwidth selection calculation from McCrary (2008))
#' @param verbose logical flag specifying whether to print diagnostic information to the terminal. (defaults to \code{FALSE})
#' @param plot logical flag indicating whether to plot the histogram and density estimations (defaults to \code{TRUE}). The user may wrap this function in additional graphical options to modify the plot.
#' @param ext.out logical flag indicating whether to return extended output. When \code{FALSE} (the default) \code{DCdensity} will return only the p-value of the test. When \code{TRUE}, \code{DCdensity} will return the additional information documented below.
#' @param htest logical flag indicating whether to return an \code{"htest"} object compatible with base R's hypothesis test output.
#' @return If \code{ext.out} is \code{FALSE}, only the p value will be returned. Additional output is enabled when \code{ext.out} is \code{TRUE}. In this case, a list will be returned with the following elements:
#' \item{theta}{the estimated log difference in heights at the cutpoint}
#' \item{se}{the standard error of \code{theta}}
#' \item{z}{the z statistic of the test}
#' \item{p}{the p-value of the test. A p-value below the significance threshhold indicates that the user can reject the null hypothesis of no sorting.}
#' \item{binsize}{the calculated size of bins for the test}
#' \item{bw}{the calculated bandwidth for the test}
#' \item{cutpoint}{the cutpoint used}
#' \item{data}{a dataframe for the binning of the histogram. Columns are \code{cellmp} (the midpoints of each cell) and \code{cellval} (the normalized height of each cell)}
#' @references McCrary, Justin. (2008) "Manipulation of the running variable in the regression discontinuity design: A density test," \emph{Journal of Econometrics}. 142(2): 698-714. \url{http://dx.doi.org/10.1016/j.jeconom.2007.05.005}
#' @include kernelwts.R
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>
#' @examples
#' #No discontinuity
#' x<-runif(1000,-1,1)
#' DCdensity(x,0)
#'
#' #Discontinuity
#' x<-runif(1000,-1,1)
#' x<-x+2*(runif(1000,-1,1)>0&x<0)
#' DCdensity(x,0)

DCdensity <- function(runvar,cutpoint,bin=NULL,bw=NULL,verbose=FALSE,plot=FALSE) {
  runvar <- runvar[complete.cases(runvar)]
  #Grab some summary vars
  rn <- length(runvar)
  rsd <- sd(runvar)
  rmin <- min(runvar)
  rmax <- max(runvar)
  if(missing(cutpoint)) {
    if(verbose) cat("Assuming cutpoint of zero.\n")
    cutpoint<-0
  }

  if(cutpoint<=rmin | cutpoint>=rmax){
   stop("Cutpoint must lie within range of runvar")
  }

  if(is.null(bin)) {
    bin <- 2*rsd*rn^(-1/2)
    if(verbose) cat("Using calculated bin size: ",sprintf("%.3f",bin),"\n")
  }

  l <- floor((rmin - cutpoint)/bin)*bin + bin/2 + cutpoint #Midpoint of lowest bin
  r <- floor((rmax - cutpoint)/bin)*bin + bin/2 + cutpoint #Midpoint of highest bin
  lc <- cutpoint-(bin/2) #Midpoint of bin just left of breakpoint
  rc <- cutpoint+(bin/2) #Midpoint of bin just right of breakpoint
  j <- floor((rmax - rmin)/bin) + 2

  binnum <- round((((floor((runvar - cutpoint)/bin)*bin + bin/2 + cutpoint) - l)/bin) + 1)

  cellval <- rep(0,j)
  for(i in seq(1,rn)){
   cnum <- binnum[i]
   cellval[cnum] <- cellval[cnum]+1
  }
  cellval <- ( cellval / rn ) / bin

  cellmp <- seq(from=1,to=j,by=1)
  cellmp <- floor(((l + (cellmp - 1)*bin ) - cutpoint)/bin)*bin + bin/2 + cutpoint

  #If no bandwidth is given, calc it
  if(is.null(bw)){
    #bin number just left of breakpoint
    leftofc <-  round((((floor((lc - cutpoint)/bin)*bin + bin/2 + cutpoint) - l)/bin) + 1)
    #bin number just right of breakpoint
    rightofc <- round((((floor((rc - cutpoint)/bin)*bin + bin/2 + cutpoint) - l)/bin) + 1)
    if ( rightofc - leftofc != 1) {
      stop("Error occurred in bandwidth calculation")
    }
    cellmpleft <- cellmp[1:leftofc]
    cellmpright <- cellmp[rightofc:j]

    #Estimate 4th order polynomial to the left
    P.lm <- lm(
      cellval ~ poly(cellmp,degree=4,raw=T),
      subset=cellmp<cutpoint
    )
    mse4 <- summary(P.lm)$sigma^2
    lcoef <- coef(P.lm)
    fppleft <- 2*lcoef[3] +
      6*lcoef[4]*cellmpleft +
      12*lcoef[5]*cellmpleft*cellmpleft
    hleft <- 3.348*(mse4*( cutpoint - l ) / sum(fppleft*fppleft))^(1/5)

    #And to the right
    P.lm <- lm(
      cellval ~ poly(cellmp,degree=4,raw=T),
      subset=cellmp>=cutpoint
    )
    mse4 <- summary(P.lm)$sigma^2
    rcoef <- coef(P.lm)
    fppright <- 2*rcoef[3] +
      6*rcoef[4]*cellmpright +
      12*rcoef[5]*cellmpright*cellmpright
    hright <- 3.348*(mse4*( r - cutpoint ) / sum(fppright*fppright))^(1/5)


    bw = .5*( hleft + hright )
    if(verbose) cat("Using calculated bandwidth: ",sprintf("%.3f",bw),"\n")
  }
  if( sum(runvar>cutpoint-bw & runvar<cutpoint) ==0 |
    sum(runvar<cutpoint+bw & runvar>=cutpoint) ==0)
    stop("Insufficient data within the bandwidth.")
  if(plot){
    #estimate density to either side of the cutpoint using a triangular kernel
    d.l<-data.frame(cellmp=cellmp[cellmp<cutpoint],cellval=cellval[cellmp<cutpoint],dist=NA,est=NA,lwr=NA,upr=NA)
    pmin<-cutpoint-2*rsd
    pmax<-cutpoint+2*rsd
    for(i in 1:nrow(d.l)) {
      d.l$dist<-d.l$cellmp-d.l[i,"cellmp"]
      w<-kernelwts(d.l$dist,0,bw,kernel="triangular")
      newd<-data.frame(dist=0)
      pred<-predict(lm(cellval~dist,weights=w,data=d.l),interval="confidence",newdata=newd)
      d.l$est[i]<-pred[1]
      d.l$lwr[i]<-pred[2]
      d.l$upr[i]<-pred[3]
    }
    d.r<-data.frame(cellmp=cellmp[cellmp>=cutpoint],cellval=cellval[cellmp>=cutpoint],dist=NA,est=NA,lwr=NA,upr=NA)
    for(i in 1:nrow(d.r)) {
      d.r$dist<-d.r$cellmp-d.r[i,"cellmp"]
      w<-kernelwts(d.r$dist,0,bw,kernel="triangular")
      newd<-data.frame(dist=0)
      pred<-predict(lm(cellval~dist,weights=w,data=d.r),interval="confidence",newdata=newd)
      d.r$est[i]<-pred[1]
      d.r$lwr[i]<-pred[2]
      d.r$upr[i]<-pred[3]
    }
    #plot to the left
    #return(list(d.l,d.r))
    plot(d.l$cellmp,d.l$est,
       lty=1,lwd=2,col="black",type="l",
       xlim=c(pmin,pmax),
       ylim=c(min(cellval[cellmp<=pmax&cellmp>=pmin]),
              max(cellval[cellmp<=pmax&cellmp>=pmin])),
       xlab=NA,
       ylab=NA,
       main=NA
    )

    lines(d.l$cellmp,d.l$lwr,
         lty=2,lwd=1,col="black",type="l"
    )
    lines(d.l$cellmp,d.l$upr,
          lty=2,lwd=1,col="black",type="l"
    )

    #plot to the right
    lines(d.r$cellmp,d.r$est,
        lty=1,lwd=2,col="black",type="l"
    )
    lines(d.r$cellmp,d.r$lwr,
          lty=2,lwd=1,col="black",type="l"
    )
    lines(d.r$cellmp,d.r$upr,
          lty=2,lwd=1,col="black",type="l"
    )

    #plot the histogram as points
    points(cellmp,cellval,type="p",pch=20)
  }
  cmp<-cellmp
  cval<-cellval
  padzeros <- ceiling(bw/bin)
  jp <- j + 2*padzeros
  if(padzeros>=1) {
    cval <- c(rep(0,padzeros),
               cellval,
               rep(0,padzeros)
             )
    cmp <- c(seq(l-padzeros*bin,l-bin,bin),
              cellmp,
              seq(r+bin,r+padzeros*bin,bin)
            )
  }

  #Estimate to the left
  dist <- cmp - cutpoint
  w <- 1-abs(dist/bw)
  w <- ifelse(w>0, w*(cmp<cutpoint), 0)
  w <- (w/sum(w))*jp
  fhatl <- predict(lm(cval~dist,weights=w),newdata=data.frame(dist=0))[[1]]

  #Estimate to the right
  w <- 1-abs(dist/bw)
  w <- ifelse(w>0, w*(cmp>=cutpoint), 0)
  w <- (w/sum(w))*jp
  fhatr<-predict(lm(cval~dist,weights=w),newdata=data.frame(dist=0))[[1]]

  #Calculate and display dicontinuity estimate
  thetahat <- log(fhatr) - log(fhatl)
  sethetahat <- sqrt( (1/(rn*bw)) * (24/5) * ((1/fhatr) + (1/fhatl)) )
  z<-thetahat/sethetahat
  p<-2*pnorm(abs(z),lower.tail=FALSE)

  # Return an htest object, for compatibility with base R test output.
  structure(list(
    statistic   = c(`z` = z),
    p.value     = p,
    method      = "McCrary (2008) sorting test",
    parameter   = c(`binwidth`  = bin,
                    `bandwidth` = bw,
                    `cutpoint`  = cutpoint),
    alternative = "no apparent sorting",
    data        = data.frame(cellmp, cellval)),
    class = c("mccrary_htest", "htest"))
}

plot.mccrary_htest <- function(x, kernel = "triangular", ...) {
    left.data  <- data.frame(subset(x$data, cellmp < x$parameter["cutpoint"]),
                             dist = NA, est = NA, lwr = NA, upr = NA)
    right.data <- data.frame(subset(x$data, cellmp >= x$parameter["cutpoint"]),
                             dist = NA, est = NA, lwr = NA, upr = NA)

    for(i in 1:nrow(left.data)) {
      left.data$dist <- left.data$cellmp - left.data[i, "cellmp"]
      w <- kernelwts(left.data$dist, 0, x$parameter["bandwidth"], kernel = kernel)
      newd <- data.frame(dist = 0)
      pred <- predict(lm(cellval ~ dist, weights = w, data = left.data),
                      interval = "confidence", newdata = newd)
      left.data$est[i] <- pred[1]
      left.data$lwr[i] <- pred[2]
      left.data$upr[i] <- pred[3]
    }

    for(i in 1:nrow(right.data)) {
      right.data$dist <- right.data$cellmp - right.data[i, "cellmp"]
      w <- kernelwts(right.data$dist, 0, x$parameter["bandwidth"], kernel = kernel)
      newd <- data.frame(dist = 0)
      pred <- predict(lm(cellval ~ dist, weights = w, data = right.data),
                      interval = "confidence", newdata = newd)
      right.data$est[i] <- pred[1]
      right.data$lwr[i] <- pred[2]
      right.data$upr[i] <- pred[3]
    }

    # Create a new, empty plot of a reasonable size.
    plot(x$data$cellmp, x$data$cellval, type = "n",
         xlim = range(x$data$cellmp),
         ylim = c(0, 2.55), #range(x$data$cellval),
         xlab = NA, ylab = NA, main = NA)

    # Plot the fit and confidence interval on the left side.
    lines(left.data$cellmp, left.data$est, type = "l",
          lty = 1, lwd = 2, col = "black")
    lines(left.data$cellmp, left.data$lwr, type = "l",
          lty = 2, lwd = 1, col = "black")
    lines(left.data$cellmp, left.data$upr, type = "l",
          lty = 2, lwd = 1, col = "black")

    # And the right.
    lines(right.data$cellmp, right.data$est, type = "l",
          lty = 1, lwd = 2, col = "black")
    lines(right.data$cellmp, right.data$lwr, type = "l",
          lty = 2, lwd = 1, col = "black")
    lines(right.data$cellmp, right.data$upr, type = "l",
          lty = 2, lwd = 1, col = "black")

    # Plot the histogram as points.
    points(x$data$cellmp, x$data$cellval, type = "p", pch = 20)
}

autoplot.mccrary_htest <- function(x, kernel = "triangular", cutpoint = TRUE,
                                   line.color = "red",
                                   point.shape = 20, point.alpha = 0.8,
                                   ...) {

    left.data  <- data.frame(subset(x$data, cellmp < x$parameter["cutpoint"]),
                             dist = NA, est = NA, lwr = NA, upr = NA)
    right.data <- data.frame(subset(x$data, cellmp >= x$parameter["cutpoint"]),
                             dist = NA, est = NA, lwr = NA, upr = NA)

    for(i in 1:nrow(left.data)) {
      left.data$dist <- left.data$cellmp - left.data[i, "cellmp"]
      w <- kernelwts(left.data$dist, 0, x$parameter["bandwidth"], kernel = kernel)
      newd <- data.frame(dist = 0)
      pred <- predict(lm(cellval ~ dist, weights = w, data = left.data),
                      interval = "confidence", newdata = newd)
      left.data$est[i] <- pred[1]
      left.data$lwr[i] <- pred[2]
      left.data$upr[i] <- pred[3]
    }

    for(i in 1:nrow(right.data)) {
      right.data$dist <- right.data$cellmp - right.data[i, "cellmp"]
      w <- kernelwts(right.data$dist, 0, x$parameter["bandwidth"], kernel = kernel)
      newd <- data.frame(dist = 0)
      pred <- predict(lm(cellval ~ dist, weights = w, data = right.data),
                      interval = "confidence", newdata = newd)
      right.data$est[i] <- pred[1]
      right.data$lwr[i] <- pred[2]
      right.data$upr[i] <- pred[3]
    }

    plt <- ggplot2::ggplot(x$data, ggplot2::aes(x = cellmp)) +
        # Adds a vertical line at the cut point.
        ggplot2::geom_vline(x = x$parameter["cutpoint"], linetype = 5,
                            colour = "gray25") +
        # The left-hand fit and interval.
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.8,
                             data = left.data, colour = NA, fill = "gray65") +
        ggplot2::geom_line(ggplot2::aes(y = est), data = left.data,
                           colour = line.color, size = 0.75,
                           lineend = "round") +
        # The right-hand fit and interval.
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.8,
                             data = right.data, colour = NA, fill = "gray65") +
        ggplot2::geom_line(ggplot2::aes(y = est), data = right.data,
                           colour = line.color, size = 0.75,
                           lineend = "round") +
        # A scatter plot of the bins themselves.
        ggplot2::geom_point(ggplot2::aes(y = cellval),
                            alpha = point.alpha, shape = point.shape) +
        # Suppress labelling.
        ggplot2::labs(x = NULL, y = NULL)

    plt
}
