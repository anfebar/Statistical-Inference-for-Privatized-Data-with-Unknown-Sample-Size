library(coda)

pdf("ATUS_figures.pdf", width = 5, height = 5)
# par(mfrow=c(2,3))
iteration <- 1

for(epsilon_ss in c(1,10))
{
  if(TRUE){
    iter <- 1
    load(paste("output/epsilonSS",epsilon_ss,
               "_iteration",iteration,
               "_iter",iter,".RData",sep = ""))
    
    SD_Bounded <- summary(mcmc.list(lapply(alpha_BoundedDPMCMCp1,mcmc)))[[1]][,2]
    
    SD_Unbounded <- sapply(OUT, function(out)
      summary(mcmc.list(lapply(out, function(l) mcmc(l[[1]]))))[[1]][,2])
    
    num_shades <- 100  # Adjust the number of shades as needed
    gray_palette <- colorRampPalette(c("black", "white"))(num_shades)[c(10, 50, 75)]
    
    
    matplot(t(SD_Unbounded), x = c(.01, .1, 1, 10), log = "xy", type = "l", 
            lty = 3, col = gray_palette, lwd = 1, ylab = "", xlab = "",
            ylim = c(0.01,5))
    abline(h = SD_Bounded, col = gray_palette, lty = 1, lwd = 2)
    if(epsilon_ss == 1)
      title(main = expression(epsilon[S] == 1))
    if(epsilon_ss == 10)
      title(main = expression(epsilon[S] == 10))
    # mtext(expression(paste("SD(",alpha[j],"|Data)")), side = 2, line = 2)
    # mtext(expression(sqrt("Var(" * alpha[j] * " | " * x[1] * ",...," * x[n] * ")")), side = 2, line = 2)
    mtext(expression(sqrt("Var(" * alpha[j] * " | DP Summaries)")), side = 2, line = 2)
    mtext(expression(epsilon[ n ]), side = 1, line = 3)
    if(epsilon_ss == 10)
      legend("topright", legend = expression(alpha[1], alpha[2], alpha[3]), 
             col = gray_palette, lty = 1, lwd = 2)
    
    
    for(iter in 2:10)
    {
      load(paste("output/epsilonSS",epsilon_ss,
                 "_iteration",iteration,
                 "_iter",iter,".RData",sep = ""))
      
      SD_Bounded <- summary(mcmc.list(lapply(alpha_BoundedDPMCMCp1,mcmc)))[[1]][,2]
      
      SD_Unbounded <- sapply(OUT, function(out)
        summary(mcmc.list(lapply(out, function(l) mcmc(l[[1]]))))[[1]][,2])
      
      matplot(t(SD_Unbounded), x = c(.01, .1, 1, 10), log = "x", type = "l",
              lty = 2, col = gray_palette, lwd = 1, add = T)
    }
  }
  
  if(TRUE){
    #######################
    #######################
    # Read data
    xF <- read.csv("female.csv")
    xM <- read.csv("male.csv")
    x_true <- rbind(xF,xM)[sample(1:(nrow(xF)+nrow(xM)),6656),]
    n <- nrow(x_true)
    
    iter <- 1
    # epsilon_ss <- 1
    
    
    load(paste("output/epsilonSS",epsilon_ss,
               "_iteration",iteration,
               "_iter",iter,".RData",sep = ""))
    
    n_mcmc <- do.call(rbind, lapply(1:4, function(j)
      cbind(epsilon_n = c(.01, .1, 1, 10)[j], 
            n = as.vector(simplify2array(lapply(OUT[[j]], function(out) out[[2]]))))))
    
    boxplot(n ~ epsilon_n, data = n_mcmc, xlab = "", ylim = c(6200, 7000))
    abline(h = n, col = 1, lwd = 2)
    mtext(expression(epsilon[ n ]), side = 1, line = 3)
    if(epsilon_ss == 1)
      title(main = expression(epsilon[S] == 1))
    if(epsilon_ss == 10)
      title(main = expression(epsilon[S] == 10))
    
    # if(epsilon_ss == 10)
    #   legend("topright", legend = "True n", 
    #          col = 3, lty = 1, lwd = 2)
    
  }
  ######################
  ######################
  if(TRUE){
    library(LaplacesDemon)
    # Generate a grid of values for the simplex
    x <- seq(0.0001, 1, by = 0.01)
    y <- seq(0.0001, 1, by = 0.01)
    
    # Create a grid of points in the simplex
    grid_points <- expand.grid(x = x, y = y)
    grid_points <- cbind(grid_points, 1-rowSums(grid_points))
    
    iter <- 1
    # epsilon_ss <- 1
    
    load(paste("output/epsilonSS",epsilon_ss,
               "_iteration",iteration,
               "_iter",iter,".RData",sep = ""))
    
    alpha <- do.call(rbind, alpha_BoundedDPMCMCp1)
    
    
    # Calculate the density values for the Dirichlet distribution at each grid point
    densities_pred <- function(alpha)
    {
      densities <- rep(0, nrow(grid_points))
      Ind <- rowSums(grid_points[,1:2]) <= 1
      densities[Ind] <- ddirichlet(grid_points[Ind,], alpha)
      densities
    }
    
    # Reshape the densities into a matrix
    density_matrix <- t(matrix(rowMeans(apply(alpha, 1, densities_pred)), nrow = length(x)))
    
    # Create a contour plot
    contour(x, y, density_matrix, xlim = c(0,0.2), ylim = c(0,0.8), ylab = expression(x[1]), xlab =  expression(x[2]))
    if(epsilon_ss == 1)
      title(main = expression(paste(epsilon[S] == 1, " and ", epsilon[n] == 0.01)))
    if(epsilon_ss == 10)
      title(main = expression(paste(epsilon[S] == 10, " and ", epsilon[n] == 0.01)))
    
    J <- 1
    alpha <- do.call(rbind, lapply(OUT[[J]],function(out) out[[1]]))
    
    # Reshape the densities into a matrix
    density_matrix <- t(matrix(rowMeans(apply(alpha, 1, densities_pred)), nrow = length(x)))
    
    # Create a contour plot
    contour(x, y, density_matrix, add = T, col = "gray", lty = 2, lwd = 2)
    
    if(epsilon_ss == 10)
      legend("topright", legend = c("Bounded", "Unbounded"), 
             col = c(1,"gray"), lty = c(1,2), lwd = 2)
  }
}
dev.off()


if(TRUE)
{
  pdf("ATUS_figures_alpha_post.pdf", width = 5*(0.32/0.48), height = 5*(0.32/0.48))
  # par(mfrow=c(2,3))
  ################################
  # alpha posterior
  ################################
  epsilon_ss <- 1
  Rxx <- list(c(0,25), c(0.5,2.5),c(5,30))
  
  load(paste("output/epsilonSS",epsilon_ss,
             "_iteration",iteration,
             "_iter",iter,".RData",sep = ""))
  
  alphaB <- do.call(rbind, alpha_BoundedDPMCMCp1)
  
  J <- 1
  alphaU <- do.call(rbind, lapply(OUT[[J]],function(out) out[[1]]))
  for(j in 1:1)
  {
    Rx <- Rxx[[j]]
    plot(density(alphaB[,j], from = Rx[1], to =  Rx[2],  n = 2500), xlim = Rx,  
         ylab = "", xlab = "", main = "", lwd = 2, ylim = c(0, 0.5))
    title(main = expression(paste(epsilon[S] == 1, " and ", epsilon[n] == 0.01)))
    if(j==1) mtext(expression(paste(alpha[1])), side = 1, line = 2)
    if(j==2) mtext(expression(paste(alpha[2])), side = 1, line = 2)
    if(j==3) mtext(expression(paste(alpha[3])), side = 1, line = 2)
    # if(j==1) mtext(expression(paste("p(",alpha[1],"|Data)")), side = 2, line = 2)
    # if(j==2) mtext(expression(paste("p(",alpha[2],"|Data)")), side = 2, line = 2)
    # if(j==3) mtext(expression(paste("p(",alpha[3],"|Data)")), side = 2, line = 2)
    mtext("Density", side = 2, line = 2.5)
    lines(density(alphaU[,j]), lwd = 2, lty = 2)
  }
  
  J <- 2
  alphaU <- do.call(rbind, lapply(OUT[[J]],function(out) out[[1]]))
  for(j in 1:1)
  {
    Rx <- Rxx[[j]]
    plot(density(alphaB[,j], from = Rx[1], to =  Rx[2],  n = 2500), xlim = Rx,  
         ylab = "", xlab = "", main = "", lwd = 2, ylim = c(0, 0.5))
    title(main = expression(paste(epsilon[S] == 1, " and ", epsilon[n] == 0.1)))
    if(j==1) mtext(expression(paste(alpha[1])), side = 1, line = 2)
    if(j==2) mtext(expression(paste(alpha[2])), side = 1, line = 2)
    if(j==3) mtext(expression(paste(alpha[3])), side = 1, line = 2)
    # if(j==1) mtext(expression(paste("p(",alpha[1],"|Data)")), side = 2, line = 2)
    # if(j==2) mtext(expression(paste("p(",alpha[2],"|Data)")), side = 2, line = 2)
    # if(j==3) mtext(expression(paste("p(",alpha[3],"|Data)")), side = 2, line = 2)
    mtext("Density", side = 2, line = 2.5)
    lines(density(alphaU[,j]), lwd = 2, lty = 2)
  }
  
  J <- 3
  alphaU <- do.call(rbind, lapply(OUT[[J]],function(out) out[[1]]))
  for(j in 1:1)
  {
    Rx <- Rxx[[j]]
    plot(density(alphaB[,j], from = Rx[1], to =  Rx[2],  n = 2500), xlim = Rx,  
         ylab = "", xlab = "", main = "", lwd = 2, ylim = c(0, 0.5))
    title(main = expression(paste(epsilon[S] == 1, " and ", epsilon[n] == 1)))
    if(j==1) mtext(expression(paste(alpha[1])), side = 1, line = 2)
    if(j==2) mtext(expression(paste(alpha[2])), side = 1, line = 2)
    if(j==3) mtext(expression(paste(alpha[3])), side = 1, line = 2)
    # if(j==1) mtext(expression(paste("p(",alpha[1],"|Data)")), side = 2, line = 2)
    # if(j==2) mtext(expression(paste("p(",alpha[2],"|Data)")), side = 2, line = 2)
    # if(j==3) mtext(expression(paste("p(",alpha[3],"|Data)")), side = 2, line = 2)
    mtext("Density", side = 2, line = 2.5)
    lines(density(alphaU[,j]), lwd = 2, lty = 2)
    
    if(j==3)
      legend("topright", legend = c("Bounded", "Unbounded"), 
             lty = c(1, 2),  lwd = 2)
  }
  
  
  ################################
  epsilon_ss <- 10
  
  load(paste("output/epsilonSS",epsilon_ss,
             "_iteration",iteration,
             "_iter",iter,".RData",sep = ""))
  
  alphaB <- do.call(rbind, alpha_BoundedDPMCMCp1)
  
  J <- 1
  alphaU <- do.call(rbind, lapply(OUT[[J]],function(out) out[[1]]))
  
  for(j in 1:1)
  {
    Rx <- Rxx[[j]]
    plot(density(alphaB[,j], from = Rx[1], to =  Rx[2],  n = 2500), xlim = Rx,  
         ylab = "", xlab = "", main = "", lwd = 2, ylim = c(0, 2))
    title(main = expression(paste(epsilon[S] == 10, " and ", epsilon[n] == 0.01)))
    if(j==1) mtext(expression(paste(alpha[1])), side = 1, line = 2)
    if(j==2) mtext(expression(paste(alpha[2])), side = 1, line = 2)
    if(j==3) mtext(expression(paste(alpha[3])), side = 1, line = 2)
    # if(j==1) mtext(expression(paste("p(",alpha[1],"|Data)")), side = 2, line = 2)
    # if(j==2) mtext(expression(paste("p(",alpha[2],"|Data)")), side = 2, line = 2)
    # if(j==3) mtext(expression(paste("p(",alpha[3],"|Data)")), side = 2, line = 2)
    mtext("Density", side = 2, line = 2.5)
    lines(density(alphaU[,j]), lwd = 2, lty = 2)
  }
  

  J <- 2
  alphaU <- do.call(rbind, lapply(OUT[[J]],function(out) out[[1]]))
  for(j in 1:1)
  {
    Rx <- Rxx[[j]]
    plot(density(alphaB[,j], from = Rx[1], to =  Rx[2],  n = 2500), xlim = Rx,  
         ylab = "", xlab = "", main = "", lwd = 2, ylim = c(0, 2))
    title(main = expression(paste(epsilon[S] == 10, " and ", epsilon[n] == 0.1)))
    if(j==1) mtext(expression(paste(alpha[1])), side = 1, line = 2)
    if(j==2) mtext(expression(paste(alpha[2])), side = 1, line = 2)
    if(j==3) mtext(expression(paste(alpha[3])), side = 1, line = 2)
    # if(j==1) mtext(expression(paste("p(",alpha[1],"|Data)")), side = 2, line = 2)
    # if(j==2) mtext(expression(paste("p(",alpha[2],"|Data)")), side = 2, line = 2)
    # if(j==3) mtext(expression(paste("p(",alpha[3],"|Data)")), side = 2, line = 2)
    mtext("Density", side = 2, line = 2.5)
    lines(density(alphaU[,j]), lwd = 2, lty = 2)
  }
  
  J <- 3
  alphaU <- do.call(rbind, lapply(OUT[[J]],function(out) out[[1]]))
  for(j in 1:1)
  {
    Rx <- Rxx[[j]]
    plot(density(alphaB[,j], from = Rx[1], to =  Rx[2],  n = 2500), xlim = Rx,  
         ylab = "", xlab = "", main = "", lwd = 2, ylim = c(0, 2))
    title(main = expression(paste(epsilon[S] == 10, " and ", epsilon[n] == 1)))
    if(j==1) mtext(expression(paste(alpha[1])), side = 1, line = 2)
    if(j==2) mtext(expression(paste(alpha[2])), side = 1, line = 2)
    if(j==3) mtext(expression(paste(alpha[3])), side = 1, line = 2)
    # if(j==1) mtext(expression(paste("p(",alpha[1],"|Data)")), side = 2, line = 2)
    # if(j==2) mtext(expression(paste("p(",alpha[2],"|Data)")), side = 2, line = 2)
    # if(j==3) mtext(expression(paste("p(",alpha[3],"|Data)")), side = 2, line = 2)
    mtext("Density", side = 2, line = 2.5)
    lines(density(alphaU[,j]), lwd = 2, lty = 2)
    
    if(j==3)
      legend("topright", legend = c("Bounded", "Unbounded"), 
             lty = c(1, 2),  lwd = 2) 
  }
  
  dev.off()
}


