##################################################################
##################################################################
### This function computes the l1-sensitivy associated with AtA
### and output a noisy version of AtA using Laplace mechanism.
##################################################################
##################################################################
# INPUTS:
### Names:    list with variables names
### Type:     list with variables types (numeric or factor) 
### Bounds:   list with upper and bounds for numerical variables
###           and levels for categorical variable
### Ref:      list with reference levels for categorical predictors
### A:        dataframe with predictors and response
### epsilon:  privacy parameter
# OUTPUTS:
### AtA:      AtA matrix
### NAtA:     sanitized AtA matrix
### SensitivyAtA: sensitivy for the entries of AtA
### EpsilonSplitAtA: indicate how to split epsilon
### QuerySplitAtA: entries to be sanitized in AtA
### Error:  error added to AtA

NoisyAtAl1 = function(A, Names, Type, Bounds, Ref, epsilon)
{
  ep_N <- epsilon$ep_N
  epsilon <- epsilon$eps 
  
  ### Applying hard-thresholding to numerical varaibles
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "numeric")
    {
      A[,Names[j]] = as.numeric(A[,Names[j]])
      A[,Names[j]] = ifelse(A[,Names[j]] < Bounds[[Names[j]]][1], Bounds[[Names[j]]][1], A[,Names[j]])
      A[,Names[j]] = ifelse(A[,Names[j]] > Bounds[[Names[j]]][2], Bounds[[Names[j]]][2], A[,Names[j]])
    }
  }
  ### Removing and adding levels to categorical predictors
  addLevel = function(x, newlevel = NULL) {
    if(is.factor(x)) {
      if (is.na(match(newlevel, levels(x))))
        return(factor(x, levels = c(levels(x), newlevel)))
    }
    return(x)
  }
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "factor")
    {
      A[,Names[j]] = as.factor(A[,Names[j]])
      Ind = A[,Names[j]] %in% Bounds[[Names[j]]]
      A = A[Ind,]
      A[,Names[j]] = droplevels(A[,Names[j]])
      Newlevels = Bounds[[Names[j]]][!(Bounds[[Names[j]]]%in%levels(A[,Names[j]]))]
      if(length(Newlevels) != 0)
        A[,Names[j]] = addLevel(A[,Names[j]],Newlevels)
      A[,Names[j]] = relevel(A[,Names[j]], ref = as.character(Ref[Names[j]]))
    }
  }
  ### If the final A (after thresholding and removing and adding levels)
  ### is empty, the following code creates a fake dataset with single row
  if(nrow(A) == 0)
  {
    for(j in 1:length(Names))
    {
      if(Type[Names[j]] == "numeric")
        A[1,Names[j]] = runif(1,Bounds[[Names[j]]][1],Bounds[[Names[j]]][2])
      if(Type[Names[j]] == "factor")
        A[1,Names[j]] = sample(Bounds[[Names[j]]],1)
    }
  }
  ### Output
  # AtA matrix
  a = model.matrix(as.formula(paste("~",paste(Names,collapse = "+"))),data = A)
  AtA = matrix(NA, ncol = ncol(a), nrow = ncol(a))
  colnames(AtA) = rownames(AtA) = colnames(a)
  # Sensitivy for AtA entries
  SensitivyAtA = AtA
  # Epsilon used at each entry of AtA
  EpsilonSplitAtA = AtA
  count.epsilon = 0
  # Entries to be sanitized in AtA
  QuerySplitAtA = AtA
  count.queries = 0
  
  ### Computing sensitivy and identifying entries to be sanitized
  # Intercept
  j = 1
  a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
  ata = (t(a)%*%a)[1,1,drop = F]
  AtA[rownames(ata),colnames(ata)] = ata
  SensitivyAtA[rownames(ata),colnames(ata)] = 1
  count.queries = count.queries + 1
  QuerySplitAtA[rownames(ata),colnames(ata)] = count.queries
  count.epsilon = count.epsilon + 1
  EpsilonSplitAtA[rownames(ata),colnames(ata)] = count.epsilon
  
  # Diagonal and first row
  for(j in 1:length(Names))
  {
    if(Type[j] == "numeric")
    {
      a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
      ata = (t(a)%*%a)
      ata = diag(ata)[-1]
      AtA[names(ata),names(ata)] = ata
      
      # Before
      # L = optim(mean(Bounds[[Names[j]]] + c(0, 0.1)), function(x) x^2, 
      #           lower = Bounds[[Names[j]]][1], upper = Bounds[[Names[j]]][2],
      #           method = "L-BFGS-B")$value
      # U = -optim(mean(Bounds[[Names[j]]] + c(0, 0.1)), function(x) -x^2, 
      #            lower = Bounds[[Names[j]]][1], upper = Bounds[[Names[j]]][2],
      #            method = "L-BFGS-B")$value
      # We replace this computation with the one below because it can be done analytically.
      tmp <- Bounds[[Names[j]]]
      if(prod(tmp) < 0) tmp <- c(0,tmp)
      U <- max(tmp^2)
      L <- min(tmp^2)
      
      SensitivyAtA[names(ata),names(ata)] = max(abs(c(U,L))) # U-L
      count.queries = count.queries + 1
      QuerySplitAtA[names(ata),names(ata)] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[names(ata),names(ata)] = count.epsilon
      
      ata = (t(a)%*%a)[1,2,drop = F]
      AtA[rownames(ata),colnames(ata)] = ata
      SensitivyAtA[rownames(ata),colnames(ata)] = max(abs(Bounds[[Names[j]]])) # diff(Bounds[[Names[j]]])
      count.queries = count.queries + 1
      QuerySplitAtA[rownames(ata),colnames(ata)] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rownames(ata),colnames(ata)] = count.epsilon
    }
    if(Type[j] == "factor")
    {
      a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
      ata = (t(a)%*%a)
      ata = diag(ata)[-1]
      
      if(length(ata) == 1)
      {
        AtA[names(ata),names(ata)] = ata
        AtA[1,names(ata)] = ata
        
        SensitivyAtA[names(ata),names(ata)] = diag(rep(1,length(ata)))
        SensitivyAtA[1,names(ata)] = 1
        
        count.queries = count.queries + 1
        QuerySplitAtA[names(ata),names(ata)] = count.queries
        QuerySplitAtA[1,names(ata)] = count.queries
        
        count.epsilon = count.epsilon + 1
        EpsilonSplitAtA[names(ata),names(ata)] = count.epsilon
        EpsilonSplitAtA[1,names(ata)] = count.epsilon
      }
      if(length(ata)>1)
      {
        AtA[names(ata),names(ata)] = diag(ata)
        AtA[1,names(ata)] = ata
        
        SensitivyAtA[names(ata),names(ata)] = diag(rep(1,length(ata)))
        SensitivyAtA[1,names(ata)] = 1
        
        count.queries = count.queries + length(ata)
        QuerySplitAtA[names(ata),names(ata)] = diag(count.queries:(count.queries-length(ata)+1))
        QuerySplitAtA[1,names(ata)] = (count.queries:(count.queries-length(ata)+1))
        
        count.epsilon = count.epsilon + 1
        EpsilonSplitAtA[names(ata),names(ata)] = diag(rep(count.epsilon,length(ata)))
        EpsilonSplitAtA[1,names(ata)] = count.epsilon
      }
    }
  }
  
  # Off-diagonal and first row
  J = combn(1:ncol(A),2)
  
  for(j in 1:ncol(J))
  {
    j1 = J[1,j]; j2 = J[2,j]
    a = model.matrix(as.formula(paste("~",paste(Names[c(j1,j2)],collapse = "+"))),data = A)
    ata = (t(a)%*%a)
    rowNames = colnames(model.matrix(as.formula(paste("~",paste(Names[c(j1,j1)],collapse = "+"))),data = A))[-1]
    colNames = colnames(model.matrix(as.formula(paste("~",paste(Names[c(j2,j2)],collapse = "+"))),data = A))[-1]
    AtA[rowNames,colNames] = ata[rowNames,colNames]
    
    if(sum(Type[Names[c(j1,j2)]] == "numeric") == 2)
    {
      # Before
      # L = optim(c(mean(Bounds[[Names[j1]]] + c(0, 0.1)),mean(Bounds[[Names[j2]]] + c(0, 0.1))), function(x) prod(x), 
      #           lower = c(Bounds[[Names[j1]]][1],Bounds[[Names[j2]]][1]), 
      #           upper = c(Bounds[[Names[j1]]][2],Bounds[[Names[j2]]][2]),
      #           method = "L-BFGS-B")$value
      # U = -optim(c(mean(Bounds[[Names[j1]]] + c(0, 0.1)),mean(Bounds[[Names[j2]]] + c(0, 0.1))), function(x) -prod(x), 
      #            lower = c(Bounds[[Names[j1]]][1],Bounds[[Names[j2]]][1]), 
      #            upper = c(Bounds[[Names[j1]]][2],Bounds[[Names[j2]]][2]),
      #            method = "L-BFGS-B")$value
      # We replace this computation with the one below because it can be done analytically.
      tmp1 <- Bounds[[Names[j1]]]
      tmp2 <- Bounds[[Names[j2]]]
      if(prod(tmp1) < 0) tmp1 <- c(0,tmp1)
      if(prod(tmp2) < 0) tmp2 <- c(0,tmp2)
      tmp <- do.call(rbind,lapply(tmp1, function(x) cbind(x,tmp2)))
      products <- apply(tmp, 1, prod)
      L <- min(products)
      U <- max(products)
      
      SensitivyAtA[rowNames,colNames] = max(abs(c(U,L))) # U-L
      count.queries = count.queries + 1
      QuerySplitAtA[rowNames,colNames] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
    
    if(sum(Type[c(j1,j2)] == "numeric") == 1)
    {
      SensitivyAtA[rowNames,colNames] = diff(Bounds[[Names[c(j1,j2)][Type[Names[c(j1,j2)]] == "numeric"]]])
      count.queries = count.queries + length(ata[rowNames,colNames])
      QuerySplitAtA[rowNames,colNames] = (count.queries:(count.queries-length(ata[rowNames,colNames])+1))
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
    
    if(sum(Type[c(j1,j2)] == "numeric") == 0)
    {
      SensitivyAtA[rowNames,colNames] = 1
      count.queries = count.queries + length(ata[rowNames,colNames])
      QuerySplitAtA[rowNames,colNames] = (count.queries:(count.queries-length(ata[rowNames,colNames])+1))
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
  }
  
  AtA = ifelse(is.na(AtA),0,AtA)
  AtA = t(AtA)+AtA*upper.tri(diag(1,ncol(AtA)))
  
  SensitivyAtA = ifelse(is.na(SensitivyAtA),0,SensitivyAtA)
  SensitivyAtA = t(SensitivyAtA)+SensitivyAtA*upper.tri(diag(1,ncol(AtA)))
  
  EpsilonSplitAtA = ifelse(is.na(EpsilonSplitAtA),0,EpsilonSplitAtA)
  EpsilonSplitAtA = t(EpsilonSplitAtA)+EpsilonSplitAtA*upper.tri(diag(1,ncol(AtA)))
  
  QuerySplitAtA = ifelse(is.na(QuerySplitAtA),0,QuerySplitAtA)
  QuerySplitAtA = t(QuerySplitAtA)+QuerySplitAtA*upper.tri(diag(1,ncol(AtA)))
  
  # Spliting epsilon
  Epsilon = epsilon/max(EpsilonSplitAtA - 1) # n does not count
  
  # Generating error
  Error = AtA-AtA
  
  j <- 1
  error = rlaplace(n = 1, location = 0, scale = max(SensitivyAtA*(QuerySplitAtA == j))/ep_N)
  Error = ifelse(QuerySplitAtA == j,error,Error)
  
  for(j in 2:max(QuerySplitAtA))
  {
    error = rlaplace(n = 1, location = 0, scale = max(SensitivyAtA*(QuerySplitAtA == j))/Epsilon)
    Error = ifelse(QuerySplitAtA == j,error,Error)
  }
  
  round(Error,3)
  
  # Noisy version of AtA
  NAtA = AtA+Error
  
  list(AtA = AtA, SensitivyAtA = SensitivyAtA, 
       EpsilonSplitAtA = EpsilonSplitAtA,
       QuerySplitAtA = QuerySplitAtA,
       Error = Error,
       NAtA = NAtA,
       A = A)
}

##################################################################
##################################################################
### This function draw error matrices using Laplace mechanism.
##################################################################
##################################################################
# INPUTs:
### N:        number of draws
### Names:    list with variables names
### Type:     list with variables types (numeric or factor) 
### Bounds:   list with upper and bounds for numerical variables
###           and levels for categorical variable
### Ref:      list with reference levels for categorical predictors
### A:        dataframe with predictors and response
### epsilon:  privacy parameter
# OUTPUTS:
### N draws from Error (NAtA = AtA + Error)
rNoisyAtAl1 = function(N, A, Names, Type, Bounds, Ref, epsilon)
{
  ep_N <- epsilon$ep_N
  epsilon <- epsilon$eps 
  
  ### Applying hard-thresholding to numerical varaibles
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "numeric")
    {
      A[,Names[j]] = as.numeric(A[,Names[j]])
      A[,Names[j]] = ifelse(A[,Names[j]] < Bounds[[Names[j]]][1], Bounds[[Names[j]]][1], A[,Names[j]])
      A[,Names[j]] = ifelse(A[,Names[j]] > Bounds[[Names[j]]][2], Bounds[[Names[j]]][2], A[,Names[j]])
    }
  }
  ### Removing and adding levels to categorical predictors
  addLevel = function(x, newlevel = NULL) {
    if(is.factor(x)) {
      if (is.na(match(newlevel, levels(x))))
        return(factor(x, levels = c(levels(x), newlevel)))
    }
    return(x)
  }
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "factor")
    {
      A[,Names[j]] = as.factor(A[,Names[j]])
      Ind = A[,Names[j]] %in% Bounds[[Names[j]]]
      A = A[Ind,]
      A[,Names[j]] = droplevels(A[,Names[j]])
      Newlevels = Bounds[[Names[j]]][!(Bounds[[Names[j]]]%in%levels(A[,Names[j]]))]
      if(length(Newlevels) != 0)
        A[,Names[j]] = addLevel(A[,Names[j]],Newlevels)
      A[,Names[j]] = relevel(A[,Names[j]], ref = as.character(Ref[Names[j]]))
    }
  }
  ### If the final A (after thresholding and removing and adding levels)
  ### is empty, the following code creates a fake dataset with single row
  if(nrow(A) == 0)
  {
    for(j in 1:length(Names))
    {
      if(Type[Names[j]] == "numeric")
        A[1,Names[j]] = runif(1,Bounds[[Names[j]]][1],Bounds[[Names[j]]][2])
      if(Type[Names[j]] == "factor")
        A[1,Names[j]] = sample(Bounds[[Names[j]]],1)
    }
  }
  ### Output
  # AtA matrix
  a = model.matrix(as.formula(paste("~",paste(Names,collapse = "+"))),data = A)
  AtA = matrix(NA, ncol = ncol(a), nrow = ncol(a))
  colnames(AtA) = rownames(AtA) = colnames(a)
  # Sensitivy for AtA entries
  SensitivyAtA = AtA
  # Epsilon used at each entry of AtA
  EpsilonSplitAtA = AtA
  count.epsilon = 0
  # Entries to be sanitized in AtA
  QuerySplitAtA = AtA
  count.queries = 0
  
  ### Computing sensitivy and identifying entries to be sanitized
  # Intercept
  j = 1
  a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
  ata = (t(a)%*%a)[1,1,drop = F]
  AtA[rownames(ata),colnames(ata)] = ata
  SensitivyAtA[rownames(ata),colnames(ata)] = 1
  count.queries = count.queries + 1
  QuerySplitAtA[rownames(ata),colnames(ata)] = count.queries
  count.epsilon = count.epsilon + 1
  EpsilonSplitAtA[rownames(ata),colnames(ata)] = count.epsilon
  
  # Diagonal and first row
  for(j in 1:length(Names))
  {
    if(Type[j] == "numeric")
    {
      a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
      ata = (t(a)%*%a)
      ata = diag(ata)[-1]
      AtA[names(ata),names(ata)] = ata
      
      # Before
      # L = optim(mean(Bounds[[Names[j]]] + c(0, 0.1)), function(x) x^2, 
      #           lower = Bounds[[Names[j]]][1], upper = Bounds[[Names[j]]][2],
      #           method = "L-BFGS-B")$value
      # U = -optim(mean(Bounds[[Names[j]]] + c(0, 0.1)), function(x) -x^2, 
      #            lower = Bounds[[Names[j]]][1], upper = Bounds[[Names[j]]][2],
      #            method = "L-BFGS-B")$value
      # We replace this computation with the one below because it can be done analytically.
      tmp <- Bounds[[Names[j]]]
      if(prod(tmp) < 0) tmp <- c(0,tmp)
      U <- max(tmp^2)
      L <- min(tmp^2)
      
      SensitivyAtA[names(ata),names(ata)] = max(abs(c(U,L))) # U-L
      count.queries = count.queries + 1
      QuerySplitAtA[names(ata),names(ata)] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[names(ata),names(ata)] = count.epsilon
      
      ata = (t(a)%*%a)[1,2,drop = F]
      AtA[rownames(ata),colnames(ata)] = ata
      SensitivyAtA[rownames(ata),colnames(ata)] = max(abs(Bounds[[Names[j]]])) # diff(Bounds[[Names[j]]])
      count.queries = count.queries + 1
      QuerySplitAtA[rownames(ata),colnames(ata)] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rownames(ata),colnames(ata)] = count.epsilon
    }
    if(Type[j] == "factor")
    {
      a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
      ata = (t(a)%*%a)
      ata = diag(ata)[-1]
      
      if(length(ata) == 1)
      {
        AtA[names(ata),names(ata)] = ata
        AtA[1,names(ata)] = ata
        
        SensitivyAtA[names(ata),names(ata)] = diag(rep(1,length(ata)))
        SensitivyAtA[1,names(ata)] = 1
        
        count.queries = count.queries + 1
        QuerySplitAtA[names(ata),names(ata)] = count.queries
        QuerySplitAtA[1,names(ata)] = count.queries
        
        count.epsilon = count.epsilon + 1
        EpsilonSplitAtA[names(ata),names(ata)] = count.epsilon
        EpsilonSplitAtA[1,names(ata)] = count.epsilon
      }
      if(length(ata)>1)
      {
        AtA[names(ata),names(ata)] = diag(ata)
        AtA[1,names(ata)] = ata
        
        SensitivyAtA[names(ata),names(ata)] = diag(rep(1,length(ata)))
        SensitivyAtA[1,names(ata)] = 1
        
        count.queries = count.queries + length(ata)
        QuerySplitAtA[names(ata),names(ata)] = diag(count.queries:(count.queries-length(ata)+1))
        QuerySplitAtA[1,names(ata)] = (count.queries:(count.queries-length(ata)+1))
        
        count.epsilon = count.epsilon + 1
        EpsilonSplitAtA[names(ata),names(ata)] = diag(rep(count.epsilon,length(ata)))
        EpsilonSplitAtA[1,names(ata)] = count.epsilon
      }
    }
  }
  
  # Off-diagonal and first row
  J = combn(1:ncol(A),2)
  
  for(j in 1:ncol(J))
  {
    j1 = J[1,j]; j2 = J[2,j]
    a = model.matrix(as.formula(paste("~",paste(Names[c(j1,j2)],collapse = "+"))),data = A)
    ata = (t(a)%*%a)
    rowNames = colnames(model.matrix(as.formula(paste("~",paste(Names[c(j1,j1)],collapse = "+"))),data = A))[-1]
    colNames = colnames(model.matrix(as.formula(paste("~",paste(Names[c(j2,j2)],collapse = "+"))),data = A))[-1]
    AtA[rowNames,colNames] = ata[rowNames,colNames]
    
    if(sum(Type[Names[c(j1,j2)]] == "numeric") == 2)
    {
      # Before
      # L = optim(c(mean(Bounds[[Names[j1]]] + c(0, 0.1)),mean(Bounds[[Names[j2]]] + c(0, 0.1))), function(x) prod(x), 
      #           lower = c(Bounds[[Names[j1]]][1],Bounds[[Names[j2]]][1]), 
      #           upper = c(Bounds[[Names[j1]]][2],Bounds[[Names[j2]]][2]),
      #           method = "L-BFGS-B")$value
      # U = -optim(c(mean(Bounds[[Names[j1]]] + c(0, 0.1)),mean(Bounds[[Names[j2]]] + c(0, 0.1))), function(x) -prod(x), 
      #            lower = c(Bounds[[Names[j1]]][1],Bounds[[Names[j2]]][1]), 
      #            upper = c(Bounds[[Names[j1]]][2],Bounds[[Names[j2]]][2]),
      #            method = "L-BFGS-B")$value
      # We replace this computation with the one below because it can be done analytically.
      tmp1 <- Bounds[[Names[j1]]]
      tmp2 <- Bounds[[Names[j2]]]
      if(prod(tmp1) < 0) tmp1 <- c(0,tmp1)
      if(prod(tmp2) < 0) tmp2 <- c(0,tmp2)
      tmp <- do.call(rbind,lapply(tmp1, function(x) cbind(x,tmp2)))
      products <- apply(tmp, 1, prod)
      L <- min(products)
      U <- max(products)
      
      SensitivyAtA[rowNames,colNames] = max(abs(c(U,L))) # U-L
      count.queries = count.queries + 1
      QuerySplitAtA[rowNames,colNames] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
    
    if(sum(Type[c(j1,j2)] == "numeric") == 1)
    {
      SensitivyAtA[rowNames,colNames] = diff(Bounds[[Names[c(j1,j2)][Type[Names[c(j1,j2)]] == "numeric"]]])
      count.queries = count.queries + length(ata[rowNames,colNames])
      QuerySplitAtA[rowNames,colNames] = (count.queries:(count.queries-length(ata[rowNames,colNames])+1))
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
    
    if(sum(Type[c(j1,j2)] == "numeric") == 0)
    {
      SensitivyAtA[rowNames,colNames] = 1
      count.queries = count.queries + length(ata[rowNames,colNames])
      QuerySplitAtA[rowNames,colNames] = (count.queries:(count.queries-length(ata[rowNames,colNames])+1))
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
  }
  
  AtA = ifelse(is.na(AtA),0,AtA)
  AtA = t(AtA)+AtA*upper.tri(diag(1,ncol(AtA)))
  
  SensitivyAtA = ifelse(is.na(SensitivyAtA),0,SensitivyAtA)
  SensitivyAtA = t(SensitivyAtA)+SensitivyAtA*upper.tri(diag(1,ncol(AtA)))
  
  EpsilonSplitAtA = ifelse(is.na(EpsilonSplitAtA),0,EpsilonSplitAtA)
  EpsilonSplitAtA = t(EpsilonSplitAtA)+EpsilonSplitAtA*upper.tri(diag(1,ncol(AtA)))
  
  QuerySplitAtA = ifelse(is.na(QuerySplitAtA),0,QuerySplitAtA)
  QuerySplitAtA = t(QuerySplitAtA)+QuerySplitAtA*upper.tri(diag(1,ncol(AtA)))
  
  # Spliting epsilon
  Epsilon = epsilon/max(EpsilonSplitAtA - 1) # n does not count
  
  # Generating error
  ftmp = function(iter)
  {
    Error = AtA-AtA
    
    j <- 1
    error = rlaplace(n = 1, location = 0, scale = max(SensitivyAtA*(QuerySplitAtA == j))/ep_N)
    Error = ifelse(QuerySplitAtA == j,error,Error)
    
    for(j in 2:max(QuerySplitAtA))
    {
      error = rlaplace(n = 1, location = 0, scale = max(SensitivyAtA*(QuerySplitAtA == j))/Epsilon)
      Error = ifelse(QuerySplitAtA == j,error,Error)
    }
    
    Error
  }
  mclapply(1:N, ftmp, mc.cores = 1)
}


##################################################################
##################################################################
### This function computes differentially private regression 
### coefficients, confidence intervals, and standard errors using 
### the Laplace mechanism.
### This implementation is an adaptation of the approach proposed
### by Ferrando, Cecilia, Shufan Wang, and Daniel Sheldon. 
### "General-Purpose Differentially-Private Confidence Intervals." 
### arXiv preprint arXiv:2006.07749 (2020).
##################################################################
##################################################################
# INPUTS:
### Names:    list with variables names (last variable should be
###           be the response)
### Type:     list with variables types (numeric or factor) 
### Bounds:   list with upper and bounds for numerical variables
###           and levels for categorical variable
### Ref:      list with reference levels for categorical predictors
### Dataset:  confidential data (data.frame)
### epsilon:  privacy parameter
### alpha:    significance level
# OUTPUTS:
### asymptotic: DP estimates using asymptotic+plug-in argument 
### bootstrap:  DP estimates using bootstrap
### bootsBeta:  draws for the bootstrap
### epsilon:  privacy parameter
### NAtA:     sanitized AtA matrix
### sensitivity: sensitivy for the entries of AtA
### epsilonSplit: indicate how to split epsilon

DpBoots.lm = function(Dataset, Names, Type, Bounds, Ref, epsilon, alpha)
{
  ### Data
  A = cbind(Dataset[,Names])
  
  ### Compute the l1-sensitivy and noisy version of AtA.
  out = NoisyAtAl1(A, Names, Type, Bounds, Ref, epsilon)
  p = nrow(out$NAtA)-1 # number of predictors
  
  ##### Summary statistics
  AtA  = out$AtA
  
  ##### Sensitivity (it is over estimated)
  SensitivyAtA = out$SensitivyAtA
  
  ##### Noisy statistics
  NAtA = out$NAtA
  p = ncol(NAtA)-1
  
  tilde.tXX = NAtA[1:p,1:p]
  hat.n = NAtA[1,1]
  hat.Q = tilde.tXX/hat.n
  
  ##### Regularization
  R = diag(0,p+1)
  Noise = rNoisyAtAl1(10000, A, Names, Type, Bounds, Ref, epsilon)
  minEV = sapply(Noise,function(x)min(eigen(x)$val))
  
  # regularization proportional to the sensitivity and quantile of the
  # min eigen value
  if(!is.positive.definite(NAtA+R))
  {
    # R = diag(diag(out$SensitivyAtA))*abs(quantile(minEV,0.01))/sum(diag(out$SensitivyAtA))
    # R = diag(abs(quantile(minEV, 0.01)), nrow(NAtA))
    R = diag(2*abs(min(eigen(NAtA)$val)), nrow(NAtA))
    # R = diag(R)
  }
  
  # if(length(R)==1) R = diag(R,p+1)
  # if(length(R)>1 & !is.matrix(R)) R = diag(R)
  # regularization proportonal equal to 3 times the eigen value of NAtA
  # if(!is.positive.definite(NAtA+R))
    # R = 3*abs(min(eigen(NAtA)$val))
  
  ##### Noisy statistics Beta
  # if(length(R)==1) R = diag(R,p+1)
  # if(length(R)>1 & !is.matrix(R)) R = diag(R)
  RNAtA = NAtA+R
  Rtilde.tXX = RNAtA[1:p,1:p]
  Rhat.n = RNAtA[1,1]
  # Rhat.n = max(NAtA[1,1],p+1)
  Rhat.n.sigma <- RNAtA[1,1]
  Rhat.Q = Rtilde.tXX/Rhat.n
  
  ##### Asymptotic Beta and sigma2 (Plug-in)
  hat.Beta = solve(Rtilde.tXX)%*%RNAtA[1:p,p+1]
  
  hat.sigma2 = ((RNAtA[p+1,p+1]-RNAtA[1:p,p+1]%*%
                   solve(Rtilde.tXX)%*%
                   RNAtA[1:p,p+1])/((Rhat.n.sigma)-p))[1,1]
  
  hat.l = hat.Beta-qnorm(1-alpha/2)*sqrt(hat.sigma2*diag(solve(Rtilde.tXX)))
  hat.u = hat.Beta+qnorm(1-alpha/2)*sqrt(hat.sigma2*diag(solve(Rtilde.tXX)))
  
  hat.se = sqrt(hat.sigma2*diag(solve(Rtilde.tXX)))
  
  ##### Bootstrap
  fboots = function(iter,Noise)
  {
    # N1 = rlaplace(length(diagAtA),0,DeltaAtA/epsilon)
    # N2 = rlaplace(length(unique_offdiagAtA),0,DeltaAtA/epsilon)
    # N = diag(N1)
    # for(j in 1:length(unique_offdiagAtA)) 
    # N  = (AtA==unique_offdiagAtA[j])*N2[j] + N
    N = Noise[[iter]]
    
    tilde.inv_n.V = N[1:p,1:p] # Noise for tXX
    tilde.inv_n.V = tilde.inv_n.V/Rhat.n
    tilde.inv_n.W = N[1:p,p+1]/Rhat.n
    tilde.inv_n.tXu = mvrnorm(1,mu = rep(0,p), Sigma = hat.sigma2*Rhat.Q/Rhat.n)
    
    tilde.Beta = solve(Rhat.Q)%*%(Rhat.Q - tilde.inv_n.V)%*%hat.Beta + 
      solve(Rhat.Q)%*%(tilde.inv_n.tXu+tilde.inv_n.W)
    
    tilde.Beta
  }
  
  # tilde.BETA = sapply(1:10000,function(iter) try(fboots(iter,Noise)))
  tilde.BETA = mclapply(1:10000,function(iter) try(fboots(iter,Noise)), 
                        mc.cores = 1)
  
  if(is.list(tilde.BETA))
  {
    Ind = which(simplify2array(lapply(tilde.BETA,is.numeric)))
    tilde.BETA = sapply(Ind,function(iter) tilde.BETA[[iter]])
  }
  
  tilde.Beta = rowMeans(tilde.BETA)
  tilde.l = apply(tilde.BETA,1,quantile,prob = (alpha/2))
  tilde.u = apply(tilde.BETA,1,quantile,prob = (1-alpha/2))
  tilde.se = apply(tilde.BETA,1,sd)
  
  # Outputs
  asymptotic = data.frame(Beta = hat.Beta, CI.l = hat.l, CI.u = hat.u, se = hat.se)
  boots = data.frame(Beta = tilde.Beta, CI.l = tilde.l, CI.u = tilde.u, se = tilde.se)
  rownames(boots) = rownames(asymptotic)
  bootsBeta = t(tilde.BETA)
  colnames(bootsBeta) = row.names(asymptotic)
  
  list(asymptotic = asymptotic,
       bootstrap = boots,
       bootsBeta = bootsBeta,
       NAtA = NAtA,
       RNAtA = RNAtA,
       epsilon = epsilon,
       sensitivity = out$SensitivyAtA,
       epsilonSplit = out$EpsilonSplitAtA,
       hat.sigma2 = hat.sigma2,
       Rndp = Rhat.n.sigma)
}