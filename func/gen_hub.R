edge.gen = function(p, sparsity, hubnumber, hubsparsity){
  G = diag(p)
  
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if(runif(1)>sparsity) G[i,j] = 1
    }
  }
  
  hubcol = sample(1:p, hubnumber, replace = FALSE)
  G[, hubcol] = rbinom(hubnumber * p, 1, 1 - hubsparsity) 
  diag(G) = 0
  G = (G + t(G))
  G[G==2] = 1
  return(list(G = G, hubindex = hubcol))
}

# Function to generate Omega1 to OmegaM from E1,...,EM, given G.list[[1]] to[[M]]
Omega.gen = function(G, M){
  p = nrow(G)
  Omega = matrix(0, nrow = p*M, ncol = p*M)
  for(l in 1:M){
    Omega.l = diag(p)
    for(j in 1:(p-1)){
      for(i in (j+1):p){
        if(G[i,j]==1){
          u.D = runif(1, 0.25, 0.75)
          if(runif(1)<0.5) u.D = -u.D
          Omega.l[i,j] = u.D
        }
      }
    }
    Omega.l = Omega.l + t(Omega.l)
    diag(Omega.l) = 0
    ee = eigen(Omega.l, only.values = T)$values[p]
    diag(Omega.l) = ifelse(ee < 0, -ee + 0.1, 0.1)
    Omega[(l-1)*p + 1:p, (l-1)*p + 1:p] = Omega.l
  }
  return(Omega)
}

# Function to generate Omega1 to OmegaM from E1,...,EM, given G.list[[1]] to[[M]]
Omega.gen.new = function(G, m, delta = 0.2){
  p = nrow(G)
  Omega = matrix(0, nrow = p*m, ncol = p*m)
  for(i in 1:p){
    index.i = (((i-1)*m+1):(i*m))
    for(j in 1:p){
      index.j = (((j-1)*m+1):(j*m))
      Omega[index.i, index.j] = 0.3*diag(m)
    }
  }
  
  ee = eigen(Omega, only.values = T)$values[p]
  diag(Omega) = ifelse(ee < 0, -ee + 1, 1)
  
  return(Omega)
}


hub.trans = function(Omega, p){
  M = dim(Omega)[1]/p
  Omega.new = matrix(0, M*p, M*p)
  for(i in 1:p){
    for(j in 1:p){
      for(s in 1:M){
        for(t in 1:M){
          Omega.new[(i-1)*M+s, (j-1)*M+t] = Omega[(s-1)*p+i,(t-1)*p+j]
        }
      }
    }
  }
  return(Omega.new)
}

