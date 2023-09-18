#####
## Helper functions for implementing a 'dumb' k-means algorithm that shows the process
####

#implement own k-means function
doKMeans = function(data, startCentroids, k = nrow(startCentroids), maxIterations = 10, distance = "euclidean", stopThreshold = 0.001) {
  centroidList = vector("list", (maxIterations+1))
  clusterMembershipList  = vector("list", (maxIterations+1))
  currentCentroids = startCentroids
  
  dimensions = ncol(data)
  dataWithClusterAssignmentPerIteration = data
  
  for (i in seq(1:maxIterations)) {
    
    
    
    distancesList = vector("list", length = k)
    for (clusterKDistance in seq(1:k)) { 
      
      
      distVec = numeric(length = nrow(data))
      if (distance == 'euclidean') {
        
        for (row in seq(1:nrow(data))) {
          distThisRow  = rowSums((data[row, ] - currentCentroids[clusterKDistance, ])^2)
          distVec[row] = distThisRow
        }
        
      }
      
      if (distance != "euclidean") {
        
        stop("Sorry, I didn't implement any other distance measure than euclidean for this function. Feel free to experiment with R's built-in k-means function later. Be sure to set distance to euclidean")
      }
      
      distancesList[[clusterKDistance]] = distVec
      
    }
    #assign all points to the centroid they are closest to and save this assignment in a list
    distancesToAllCentroids = dplyr::bind_cols(distancesList)
    clusterAssignment       = apply(distancesToAllCentroids, 1, which.min)
    
    #add in initial situation
    if (i == 1) {
      clusterMembershipList[[1]] = clusterAssignment
      centroidList[[1]]          = startCentroids
      dataWithClusterAssignmentPerIteration[[paste0("clusterMembershipIteration", 0)]] = factor(clusterAssignment, levels = as.character(seq(1:k)))
    }
    
    
    dataWithClusterAssignmentPerIteration[[paste0("clusterMembershipIteration", i)]] = factor(clusterAssignment, levels = as.character(seq(1:k)))
    clusterMembershipList[[i+1]] = clusterAssignment
    
    
    #recalculate the cluster centroids and save this assignment in a list
    newCentroids = dataWithClusterAssignmentPerIteration %>%
      group_by(across(all_of(paste0("clusterMembershipIteration", i)))) %>%
      summarize(across(.cols = 1:dimensions, mean))
    currentCentroids = newCentroids[, 2:(dimensions+1)]
    centroidList[[i+1]] = currentCentroids  
    
    
    
    
    
    #check for convergence. If total difference between centroid coordinates between iterations is
    #below certain threshold, stop.
    if (i > 1) {
      
      differenceInCentroidsBetweenIterations = abs(sum(centroidList[[i]] - centroidList[[i-1]]))
      
      if(differenceInCentroidsBetweenIterations < stopThreshold) {
        break
      }
      
    }
    
    
  }
  list(dataPointAssignedClusters = purrr::compact(clusterMembershipList),
       clusterCentroidPositions  = purrr::compact(centroidList),
       dataAllClusterAssignment  = dataWithClusterAssignmentPerIteration)
  
  
}


#plot initial situation + true cluster membership
preClusterPlot   = function(data, centers) {
  bgPlot = ggplot(data = dataforClustering, aes_string(x = "`Gene 1 Expression`", y = "`Gene 2 Expression`", colour = "`True Cluster`")) +
    geom_point(size = 3) +
    geom_point(size = 5, stroke = 2.5, shape = 21, data = centers, colour = "black",
               aes_string(fill = "`True Cluster`"),
               show.legend = FALSE) +
    scale_fill_discrete(drop = FALSE) +
    theme_bw() +
    ggtitle("True underlying data and starting centroids")
  bgPlot
}

#plot the cluster membership and centroids for each iteration of k-means
iterationPlot = function(dataPoints, dataCenters) {
  
  #make plots according to number of iterations in data
  plotsToMake = dataPoints %>% select(contains("Iteration")) %>% ncol()
  colNames = dataPoints %>% select(contains("Iteration")) %>% colnames()
  itPlotList = vector("list", length = plotsToMake)
  for (plotNr in seq(1:plotsToMake)) {
    
    dataWithCentroidIdentity = dataCenters[[plotNr]] %>% mutate(Centroid = as.factor(1:nrow(dataCenters[[plotNr]])))
    itPlot = ggplot(data = dataPoints,
                    aes_string(x = "`Gene 1 Expression`", y = "`Gene 2 Expression`",
                               colour = colNames[plotNr])
    ) +
      geom_point(size = 3) +
      geom_point(size = 5, stroke = 2.5, shape = 21, data = dataWithCentroidIdentity, colour = "black",
                 aes_string(fill = "Centroid"),
                 show.legend = NA) +
      theme_bw() +
      ggtitle(paste0("iteration: ", plotNr))
    itPlotList[[plotNr]] = itPlot
    
  }
  
  itPlotList
}

#make it animated
iterationPlotAnim = function(dataPoints, dataCenters) {
  
  #reshape data for animation
  nIterations                 = dataPoints %>% select(contains("Iteration")) %>% ncol()
  dimData                     = dataPoints %>% select(-contains("Iteration")) %>% ncol()
  dataPointsLong              = dataPoints %>% pivot_longer((dimData+1):(nIterations+dimData)) %>% arrange(name)
  dataPointsLong$name         = as.factor(readr::parse_number(dataPointsLong$name))
  colnames(dataPointsLong)[3] = "iteration"
  colnames(dataPointsLong)[4] = "cluster membership"
  #dataPointsLong %<>% select(-iteration)
  
  #reshape centroids for animation
  centroidFrame           = dplyr::bind_rows(dataCenters)
  centroidFrame$iteration = as.factor(rep((seq(1:(nIterations))-1), each = nrow(dataCenters[[1]])))
  centroidFrame$centroid  = as.factor(rep(seq(1:nrow(dataCenters[[1]])), (nIterations) ))
  
  
  itPlot = ggplot(data = dataPointsLong,
                  aes_string(x = "`Gene 1 Expression`", y = "`Gene 2 Expression`",
                             colour = "`cluster membership`")
  ) +
    geom_point(size = 3) +
    geom_point(size = 5, stroke = 2.5, shape = 21, data = centroidFrame, colour = "black",
               aes_string(fill = "centroid"),
               show.legend = NA) +
    theme_bw() +
    transition_manual(`iteration`) +
    ease_aes() +
    labs(title = "Iteration: {current_frame}")
  itPlot
}

#start with randomStarts random initialisations, picking k random centroids from the data each time
getKRandomStartPoints = function(data, k = 4, randomStarts = 4) {
  
  centersForRuns = vector("list", randomStarts)
  
  if (k <= 0) {
    stop("You can't use 0 or fewer clusters, unfortunately. Set k to >=1")
  }
  
  if (k > nrow(data)) {
    stop("You can't have more clusters than there are data points in the data. That doesn't make sense.")
  }
  
  for( i in seq_along(centersForRuns)) {
    centersForRuns[[i]] = data[sample(seq(1:nrow(data)), k), ]
  }
  centersForRuns
  
}

#calculate kMeans clustering steps for each random initialisation, output in list
#note: in this form, function won't work for multi-dim data
calculateKMeansTrajectoryEachRandomInit = function(data = dataforClustering, initialCentroids = centersForRuns) {
  
  kMeansOutcomes = vector("list", length = length(initialCentroids))
  for(randomCentroidInitialisationNr in seq_along(initialCentroids)) {
    kMeansOutcomes[[randomCentroidInitialisationNr]] = doKMeans(data[, 1:2],
                                                                initialCentroids[[randomCentroidInitialisationNr]][, 1:2])
  }
  kMeansOutcomes
}

#generate plots and animations using the plotting and animating function defined above
makePlotsAndAnimation = function(kMeansOutcomeData = kMeansOutcomes, centers = centersForRuns) {
  
  animationList = vector("list", length(kMeansOutcomeData))
  totalPlotList = vector("list", length(kMeansOutcomeData))
  for(i in seq_along(kMeansOutcomeData)) {
    
    beginPlot      = preClusterPlot(data = dataforClustering, centers = centers[[i]])
    iterationPlots = iterationPlot(dataPoints  = kMeansOutcomeData[[i]]$dataAllClusterAssignment,
                                   dataCenters = kMeansOutcomeData[[i]]$clusterCentroidPositions)
    animationPlot  = iterationPlotAnim(dataPoints  = kMeansOutcomeData[[i]]$dataAllClusterAssignment,
                                       dataCenters = kMeansOutcomeData[[i]]$clusterCentroidPositions)
    totalPlotList[[i]] = list(beginPlot)
    totalPlotList[[i]] = append(totalPlotList[[i]], iterationPlots)
    animationList[[i]] = animationPlot
    
  }
  
  return(list(plots = totalPlotList, animations = animationList))
  
}




##ggbiplot code, taken from: https://github.com/vqv/ggbiplot:


ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE, 
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = muted('red'))
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  # Label the variable axes
  if(var.axes) {
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar, y = yvar, 
                    angle = angle, hjust = hjust), 
                color = 'darkred', size = varname.size)
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}

ggscreeplot <- function(pcobj, type = c('pev', 'cev')) 
{
  type <- match.arg(type)
  d <- pcobj$sdev^2
  yvar <- switch(type, 
                 pev = d / sum(d), 
                 cev = cumsum(d) / sum(d))
  
  yvar.lab <- switch(type,
                     pev = 'proportion of explained variance',
                     cev = 'cumulative proportion of explained variance')
  
  df <- data.frame(PC = 1:length(d), yvar = yvar)
  
  ggplot(data = df, aes(x = PC, y = yvar)) + 
    xlab('principal component number') + ylab(yvar.lab) +
    geom_point() + geom_path()
}



