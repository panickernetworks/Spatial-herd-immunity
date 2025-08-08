#CODE FOR SPATIAL HERD IMMUNITY PAPER

#numerical solution of AI model for fixed k
kSIR = function(N, b, k,f) {
  l = 1
  tmax = 120
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)   #Density       
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta)**(I[h]*(rho/N)*pi*b*(b/(f*f))*(
      1+(f*f-1)*(1-k))))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}


disk_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}


copy <- function(x) {
  y <- c()
  for (i in 1:length(x)){
    y <- c(y, x[i])
  }
  return(y)
}

#AI model, fixed k, Fully mixed
kSIRSIM <- function(N, r1,k,f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initially_infected <- sample(1:N, 1)
  node_types[initially_infected] <- "I"
  prevalence <- numeric(num_steps)
  infected_nodes <- which(node_types == "I")
  node_r = rep(r1,N)
  r2=r1/f
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2 
  }
  x=runif(N,0,1)
  y=runif(N,0,1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N #append infected
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- which((prev_node_types == "S" & dist < r1))
        }else{
          neighbors <- which((prev_node_types == "S" & dist < r2))
        }
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    node_r =rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes <- infected_nodes[runif(length(infected_nodes)) < k]  
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
    x=runif(N,0,1)
    y=runif(N,0,1)
  }
  return(prevalence)  
}

#AI model, fixed k, static
kSIRSTI <- function(N, r1,k,f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initially_infected <- sample(1:N, 1)
  node_types[initially_infected] <- "I"
  prevalence <- numeric(num_steps)
  infected_nodes <- which(node_types == "I")
  node_r = rep(r1,N)
  r2=r1/f
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2  
  }
  x=runif(N,0,1)
  y=runif(N,0,1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- which((prev_node_types == "S" & dist < r1))
        }else{
          neighbors <- which((prev_node_types == "S" & dist < r2))
        }
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    node_r =rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes <- infected_nodes[runif(length(infected_nodes)) < k]  
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)
}


#AS model, fixed k, fully mixed
kSIRSIM2 <- function(N, r1,k,f) {
  beta=0.5
  num_steps <- 120
  node_types <- rep("S", N)
  initially_infected <- sample(1:N, 1)
  node_types[initially_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  node_r= rep(r1,N)
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N #append infected
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "S") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- length(which(prev_node_types == "I" & dist < r1))
        }else{
          neighbors <- length(which(prev_node_types == "I" & dist < r2))
        }
        if(neighbors==0){
          node_types[i] =="S"
        }
        else{
          prob=1-(1-beta)**(neighbors)
          if (runif(1) < prob) {
            node_types[i] <- "I"
          }
        }
      }
      if(prev_node_types[i] == "I"){
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    node_r= rep(r1,N)
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
    x <- runif(N, 0, 1)
    y <- runif(N, 0, 1)
  }
  return(prevalence)
}

#AS model, fixed k, static
kSIRSTI2 <- function(N, r1,k,f) {
  beta=0.5
  num_steps <- 120
  node_types <- rep("S", N)
  initially_infected <- sample(1:N, 1)
  node_types[initially_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  node_r= rep(r1,N)
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I") 
    prevalence[time_step] <- num_infected / N #append infected
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "S") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- length(which(prev_node_types == "I" & dist < r1))
        }else{
          neighbors <- length(which(prev_node_types == "I" & dist < r2))
        }
        if(neighbors==0){
          node_types[i] =="S"
        }
        else{
          prob=1-(1-beta)**(neighbors)
          if (runif(1) < prob) {
            node_types[i] <- "I"
          }
        }
      }
      if(prev_node_types[i] == "I"){
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    node_r= rep(r1,N)
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)
}

#Numerical solution of AS model for fixed k
kSIR2 = function(N, b, k,f) {
  l = 1
  tmax = 120
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density 
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] *(1-k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b))+
      S[h] *(k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b/(f*f)))    
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}

#AIS model, fixed k, fully mixed
kSIRSIM3 <- function(N, r1,k,f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  node_r= rep(r1,N)
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2
  }
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N #append infected
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        # Calculate distances with periodic boundary
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r2){
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- c()
        }else{
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- which(prev_node_types == "S" & dist < r1 & dist > r2 & node_r==r1)
        }
        neighbors=c(neighbors1,neighbors2)
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    node_r= rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes1 <-  infected_nodes[runif(length(infected_nodes)) < k] 
      node_r[adapted_nodes1] <- r2  
    }else{
      node_r = node_r
    }
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
    x <- runif(N, 0, 1)
    y <- runif(N, 0, 1)
  }
  return(prevalence)
}

#AIS model, fixed k, static
kSIRSTI3 <- function(N, r1,k,f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  node_r= rep(r1,N)
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2
  }
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N #append infected
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        # Calculate distances with periodic boundary
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r2){
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- c()
        }else{
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- which(prev_node_types == "S" & dist < r1 & dist > r2 & node_r==r1)
        }
        neighbors=c(neighbors1,neighbors2)
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    node_r= rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes1 <-  infected_nodes[runif(length(infected_nodes)) < k] 
      node_r[adapted_nodes1] <- r2  
    }else{
      node_r = node_r
    }
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)
}

#Numerical solution of AIS model for fixed k
kSIR3 = function(N, b, k,f) {
  l = 1
  tmax = 120
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density 
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] *(1-k)*(1 - (1 - beta)**((I[h]*(rho/N)*pi*(b*b)/(f*f))+
                                          (I[h]*(rho/N)*pi*((b*b) - (b*b)/(f*f))*(1-k))))+
      S[h] *(k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b/(f*f)))   
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}
#############################################################################
#Numerical solution: AI mode;power law function
SIRm = function(N, b,m,f,tmax) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density 
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    k=(I[h]/N)**m
    dS <- S[h] * (1 - (1 - beta)**(I[h]*(rho/N)*pi*b*(b/(f*f))*(
      1+(f*f-1)*(1-k))))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}

#Numerical solution: AS mode;power law function
SIRm2 = function(N, b,m,f,tmax) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density 
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    k=(I[h]/N)**m
    dS <- S[h] *(1-k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b))+
      S[h] *(k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b/(f*f)))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}

#Numerical solution: AIS mode;power law function
SIRm3 = function(N, b,m,f,tmax) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density 
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    k=(I[h]/N)**m
    dS <- S[h] *(1-k)*(1 - (1 - beta)**((I[h]*(rho/N)*pi*(b*b)/(f*f))+
                                          (I[h]*(rho/N)*pi*((b*b) - (b*b)/(f*f))*(1-k))))+
      S[h] *(k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b/(f*f)))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}
################################################
#AI mode, power law, static
SIRmSTI <- function(N, r1, m, f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  infected_nodes <- which(node_types == "I")
  node_r = rep(r1,N)
  r2=r1/f
  k = (1/N)**m
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2  # Apply adaptation with probability p
  }
  x=runif(N,0,1)
  y=runif(N,0,1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N #append infected
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- which((prev_node_types == "S" & dist < r1))
        }else{
          neighbors <- which((prev_node_types == "S" & dist < r2))
        }
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    s = sum(node_types == "I") / N
    k = s**m
    node_r =rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes <- infected_nodes[runif(length(infected_nodes)) < k]  
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)  # Return the maximum prevalence during the simulation
}

##################################################
#AS mode, power law, static
SIRmSTI2 <- function(N, r1,m,f) {
  beta=0.5
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  node_r= rep(r1,N)
  k = (1/N)**m
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N#append infected
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "S") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- length(which(prev_node_types == "I" & dist < r1))
        }else{
          neighbors <- length(which(prev_node_types == "I" & dist < r2))
        }
        if(neighbors==0){
          node_types[i] =="S"
        }
        else{
          prob=1-(1-beta)**(neighbors)
          if (runif(1) < prob) {
            node_types[i] <- "I"
          }
        }
      }
      if(prev_node_types[i] == "I"){
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    s = sum(node_types == "I") / N
    k = s**m
    node_r= rep(r1,N)
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)
}
#################################################
#AIS mode power law static
SIRmSTI3 <- function(N, r1, m, f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  k = (1/N)**m
  node_r= rep(r1,N)
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2
  }
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N #append infected
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        # Calculate distances with periodic boundary
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r2){
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- c()
        }else{
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- which(prev_node_types == "S" & dist < r1 & dist > r2 & node_r==r1)
        }
        neighbors=c(neighbors1,neighbors2)
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    s = sum(node_types == "I") / N
    k = s**m
    node_r= rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes1 <-  infected_nodes[runif(length(infected_nodes)) < k] 
      node_r[adapted_nodes1] <- r2  
    }else{
      node_r = node_r
    }
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)
}
######################################################
#Numerical solution of AI model: sigmoid law
SIRq1 = function(N, b,q,f,tmax) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  kupd=c()
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    k = 1/(1+exp(-((I[h]/N)-0.25)/q))
    dS <- S[h] * (1 - (1 - beta)**(I[h]*(rho/N)*pi*b*(b/(f*f))*(
      1+(f*f-1)*(1-k))))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
    kupd=c(kupd,k)
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}

######################################################
#Numerical solution of AS model: sigmoid law
SIRq2 = function(N, b,q,f,tmax) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density 
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  kupd=c()
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    k = 1/(1+exp(-((I[h]/N)-0.25)/q))
    dS <- S[h] *(1-k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b))+
      S[h] *(k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b/(f*f)))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
    kupd=c(kupd,k)
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}

#################################################
#Numerical solution of AIS model: sigmoid law
SIRq3 = function(N, b,q,f,tmax) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)
  beta = 0.5
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  rho = N / (l ** 2)           # Density
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  kupd=c()
  gamma = 0.05      # Recovery rate (gamma) # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    k = 1/(1+exp(-((I[h]/N)-0.25)/q))
    dS <- S[h] *(1-k)*(1 - (1 - beta)**((I[h]*(rho/N)*pi*(b*b)/(f*f))+
                                          (I[h]*(rho/N)*pi*((b*b) - (b*b)/(f*f))*(1-k))))+
      S[h] *(k)*(1 - (1 - beta)**(I[h]*(rho/N)*pi*b*b/(f*f)))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
    kupd=c(kupd,k)
  }
  # Return the list of daily infected counts (prevalence)
  return(c(1/N,inf/N)[1:tmax])
}

#AI model, sigmoid, Static
SIRSTIq <- function(N, r1, q, f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  infected_nodes <- which(node_types == "I")
  node_r = rep(r1,N)
  r2=r1/f
  k = 1/(1+exp(-((1/N)-0.25)/q))
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2  # Apply adaptation with probability p
  }
  x=runif(N,0,1)
  y=runif(N,0,1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- which((prev_node_types == "S" & dist < r1))
        }else{
          neighbors <- which((prev_node_types == "S" & dist < r2))
        }
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    s = sum(node_types == "I") / N
    k = 1/(1+exp(-(s-0.25)/q))
    node_r =rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes <- infected_nodes[runif(length(infected_nodes)) < k]  
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)   # Return the prevalence during the simulation
}

#AS model, sigmoid, Static
SIRSTIq2 <- function(N, r1,q,f) {
  beta=0.5
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  node_r= rep(r1,N)
  k = 1/(1+exp(-((1/N)-0.25)/q))
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "S") {
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r1){
          neighbors <- length(which(prev_node_types == "I" & dist < r1))
        }else{
          neighbors <- length(which(prev_node_types == "I" & dist < r2))
        }
        if(neighbors==0){
          node_types[i] =="S"
        }
        else{
          prob=1-(1-beta)**(neighbors)
          if (runif(1) < prob) {
            node_types[i] <- "I"
          }
        }
      }
      if(prev_node_types[i] == "I"){
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    s = sum(node_types == "I") / N
    k = 1/(1+exp(-(s-0.25)/q))
    node_r= rep(r1,N)
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)   # Return the  prevalence during the simulation
}


#AIS model, sigmoid, Static
SIRSTIq3 <- function(N, r1, q, f) {
  num_steps <- 120
  node_types <- rep("S", N)
  initial_infected <- sample(1:N, 1)
  node_types[initial_infected] <- "I"
  prevalence <- numeric(num_steps)
  r2=r1/f
  k = 1/(1+exp(-((1/N)-0.25)/q))
  node_r= rep(r1,N)
  if (runif(1) < k) {  
    node_r[initial_infected] <- r2
  }
  susceptible_nodes <- which(node_types == "S")  
  if(length(susceptible_nodes)>0){
    adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
    node_r[adapted_nodes] <- r2  
  }else{
    node_r = node_r
  }
  x <- runif(N, 0, 1)
  y <- runif(N, 0, 1)
  for (time_step in 1:num_steps) {
    prev_node_types <- copy(node_types)
    num_infected <- sum(node_types == "I")
    prevalence[time_step] <- num_infected / N
    neighbors <- integer(0)
    for (i in 1:N) {
      if (prev_node_types[i] == "I") {
        # Calculate distances with periodic boundary
        dx <- abs(x - x[i])
        dy <- abs(y - y[i])
        dx <- ifelse(dx > 0.5, 1 - dx, dx)
        dy <- ifelse(dy > 0.5, 1 - dy, dy)
        dist <- sqrt(dx^2 + dy^2)
        if(node_r[i] == r2){
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- c()
        }else{
          neighbors1 <- which(prev_node_types == "S" & dist < r2)
          neighbors2 <- which(prev_node_types == "S" & dist < r1 & dist > r2 & node_r==r1)
        }
        neighbors=c(neighbors1,neighbors2)
        for (neighbor in neighbors) {
          if (runif(1) < 0.5) {
            node_types[neighbor] <- "I"
          }
        }
        if (runif(1) < 0.05) {
          node_types[i] <- "R"
        }
      }
    }
    s = sum(node_types == "I") / N
    k = 1/(1+exp(-(s-0.25)/q))
    node_r= rep(r1,N)
    infected_nodes <- which(node_types == "I")  
    if(length(infected_nodes)>0){
      adapted_nodes1 <-  infected_nodes[runif(length(infected_nodes)) < k] 
      node_r[adapted_nodes1] <- r2  
    }else{
      node_r = node_r
    }
    susceptible_nodes <- which(node_types == "S")  
    if(length(susceptible_nodes)>0){
      adapted_nodes <-  susceptible_nodes[runif(length(susceptible_nodes)) < k] 
      node_r[adapted_nodes] <- r2  
    }else{
      node_r = node_r
    }
  }
  return(prevalence)   # Return the  prevalence during the simulation
}
