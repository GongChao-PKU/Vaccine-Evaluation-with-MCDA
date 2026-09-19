
file <- readxl::read_excel("param/vaccine_1.xlsx")
file_sen <- readxl::read_excel("param/vaccine_sensitivity analysis_1.xlsx")

file_D <- readxl::read_excel("param/vaccine_2.xlsx")
file_sen_D <- readxl::read_excel("param/vaccine_sensitivity analysis_2.xlsx")
# weight extraction
weight <- file[1,-1]

# cost_effectiveness calculation and point value of parameters display
param_add_cost_effect <- function(x=file){
  table_pre <- x[-c(1,2),-1]
  table_pre[table_pre$persistence==2,]$cost_effectiveness <- round(table_pre[table_pre$persistence==2,]$price/
                                                                     (table_pre[table_pre$persistence==2,]$mortality*table_pre[table_pre$persistence==2,]$effectiveness)/78,0)
  table_pre[table_pre$persistence==1,]$cost_effectiveness <- round(table_pre[table_pre$persistence==1,]$price/
                                                                     (table_pre[table_pre$persistence==1,]$mortality*table_pre[table_pre$persistence==1,]$effectiveness),0)
  table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness <- round(table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness/2,0)
  return(table_pre)
}

# CI of parameters display
param_CI <- function(x=file_sen){
  table_pre <- x[,-1]
  table_pre <- table_pre[,unlist(lapply(table_pre,function(x){x[1]!="NA"}))]
  return(table_pre)
}

# point results of AHP
compare_result_AHP <- function(point_data=file){
  table_weight <- point_data[1,-1]
  table_pre <- point_data[-c(1,2),-1]
  table_pre[table_pre$persistence==2,]$cost_effectiveness <- table_pre[table_pre$persistence==2,]$price/
    (table_pre[table_pre$persistence==2,]$mortality*table_pre[table_pre$persistence==2,]$effectiveness)/78     
  table_pre[table_pre$persistence==1,]$cost_effectiveness <- table_pre[table_pre$persistence==1,]$price/       
    (table_pre[table_pre$persistence==1,]$mortality*table_pre[table_pre$persistence==1,]$effectiveness)
  table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness <- table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness/2
  table_rank <- data.frame(
    lapply(table_pre,rank),
    row.names=rownames(table_pre)
  )
  table_pre_direction <- point_data[2,-1]
  table_rank[,which(table_pre_direction==1)] <- lapply(-table_pre[,which(table_pre_direction==1)],rank)
  
  result_AHP <- data.frame(
    score <- round(as.matrix(table_rank)%*%(t(as.matrix(table_weight))),3),
    rank <- rank(-score)
  )
  colnames(result_AHP)=c("score","rank")
  return(result_AHP)
}

# point parameter extraction (n_times times)
param_value <- function(param_type = NA,
                        param_point = 0,
                        param_low = -1,
                        param_high = 1,
                        lower_bound = -Inf,
                        upper_bound = Inf) {
  library(mc2d)
  if (param_type == "NA") {
    param_value <- rep(param_point, n_times)
  } else if (param_type == "normal") {
    sd_value <- (param_high - param_low) / (2 * 1.96)
    param_value <- rnorm(n_times,param_point,sd_value)
    # 超出理论边界的值重新抽样
    while (any(param_value < lower_bound | param_value > upper_bound)) {
      index <- which(param_value < lower_bound | param_value > upper_bound)
      param_value[index] <- rnorm(length(index),param_point,sd_value)
    }
  } else if (param_type == "lognormal") {
    param_value <- exp(rnorm(n_times,log(param_point, 10),log(param_high / param_low, 10) / (2 * 1.96)))
  } else if (param_type == "triangle") {
    param_value <- rtri(n_times,param_low,param_high,param_point)
  } else if (param_type == "uniform") {
    param_value <- runif(n_times,param_low,param_high)
  }
  return(param_value)
}

# PSA results of AHP
compare_result_AHP_graph <- function(point_data=file,interval_data=file_sen){
  table_weight <- point_data[1,-1]
  table_pre_point <- point_data[-c(1,2),-1]
  
  table_pre_type <- interval_data[1,-1]
  table_pre_ci <- interval_data[-c(1,2),-1]
  
  table_pre_all <- array(0,dim=c(nrow(table_pre_point),ncol(table_pre_point),n_times))
  for(i in 1:nrow(table_pre_point)){
    for(j in 1:ncol(table_pre_point)){
      table_pre_all[i,j,] <- param_value(param_type=table_pre_type[j],param_point=as.numeric(table_pre_point[i,j]),
                                         param_low=as.numeric(table_pre_ci[i,j]),
                                         param_high=as.numeric(table_pre_ci[i,(j+ncol(table_pre_point))]),
                                      
                                         lower_bound=as.numeric(interval_data[2,j+1]),
                                         upper_bound=as.numeric(interval_data[2,(j+1+ncol(table_pre_point))]))
    }
  }
  
  table_score <- array(0,dim=c(nrow(table_pre_point),n_times))
  table_sort <- array(0,dim=c(nrow(table_pre_point),n_times))
  rownames(table_score) <- rownames(table_sort) <- unlist(point_data[-c(1,2),1])
  for(k in 1:n_times){
    table_pre <- table_pre_point
    table_pre[1:nrow(table_pre),1:ncol(table_pre)] <- table_pre_all[,,k]
    table_pre[table_pre$persistence==2,]$cost_effectiveness <- table_pre[table_pre$persistence==2,]$price/
      (table_pre[table_pre$persistence==2,]$mortality*table_pre[table_pre$persistence==2,]$effectiveness)/78     
    table_pre[table_pre$persistence==1,]$cost_effectiveness <- table_pre[table_pre$persistence==1,]$price/       
      (table_pre[table_pre$persistence==1,]$mortality*table_pre[table_pre$persistence==1,]$effectiveness)
    table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness <- table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness/2
    table_rank <- data.frame(
      lapply(table_pre,rank),
      row.names=rownames(table_pre)
    )
    
    table_pre_direction <- point_data[2,-1]
    table_rank[,which(table_pre_direction==1)] <- lapply(-table_pre[,which(table_pre_direction==1)],rank)
    
    result_AHP <- data.frame(
      score <- as.matrix(table_rank)%*%(t(as.matrix(table_weight))),
      rank <- rank(-score)
    )
    colnames(result_AHP)=c("score","rank")
    table_score[,k] <- result_AHP$score
    table_sort[,k] <- result_AHP$rank
  }
  table_score <- reshape2::melt(table_score, id = rownames(table_score))
  library(ggplot2)
  p1 <- ggplot(table_score, aes(`Var1`,`value`,group=`Var1`)) + geom_boxplot() + 
    labs(x="vaccine",y="MCDA score") + 
    theme_classic()
    
  table_sort_pre <- matrix(0,nrow(table_pre_point),nrow(table_pre_point))
  for(i in 1:ncol(table_sort_pre)){
    for(j in 1:nrow(table_sort_pre)){
      table_sort_pre[j,i] <- sum(table_sort[j,]==i)
    }
  }
  table_sort_pre <- data.frame(table_sort_pre)
  colnames(table_sort_pre)[1:6] <- 1:6
  table_sort_pre[,'vaccine'] <- point_data[-c(1,2),1]
  table_sort <- tidyr::gather(table_sort_pre,key="rank",value="proportion",-"vaccine")
  table_sort$proportion <- as.numeric(table_sort$proportion)/n_times
  library(ggplot2)
  p2 <- ggplot(table_sort,aes(x=`vaccine`,y=`proportion`,fill=`rank`)) + 
    scale_fill_brewer(palette = "Blues")+
    geom_bar(stat = 'identity',width = 0.5,colour="black") + 
    theme_classic()
  
  return(list(p1,p2,table_score,table_sort))
}

# point results of TOPSIS
compare_result_TOPSIS <- function(point_data=file){
  table_pre <- point_data[-c(1,2),-1]
  table_pre[table_pre$persistence==2,]$cost_effectiveness <- table_pre[table_pre$persistence==2,]$price/
    (table_pre[table_pre$persistence==2,]$mortality*table_pre[table_pre$persistence==2,]$effectiveness)/78
  table_pre[table_pre$persistence==1,]$cost_effectiveness <- table_pre[table_pre$persistence==1,]$price/
    (table_pre[table_pre$persistence==1,]$mortality*table_pre[table_pre$persistence==1,]$effectiveness)
  table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness <- table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness/2
  
  for(i in 1:ncol(table_pre)){
    table_pre[(table_pre[,i]==max(table_pre[,i])),i][1,] <-
      table_pre[(table_pre[,i]==max(table_pre[,i])),i][1,]*1.001
    table_pre[(table_pre[,i]==min(table_pre[,i])),i][1,] <-
      table_pre[(table_pre[,i]==min(table_pre[,i])),i][1,]*0.999
  }
  
  # entropy of positive indicator
  entropy_positive <- function(x){
    # Indicator normalization
    y = (x -min(x)) / (max(x) - min(x))
    # Indicator weight
    p = y / sum(y)
    # Indicator entropy
    entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
  }
  
  # entropy of negative indicator
  entropy_negative <- function(x){
    # Indicator normalization
    y = (max(x) -x) / (max(x) - min(x))
    # Indicator weight
    p = y/ sum(y)
    # Indicator entropy
    entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
  }
  
  table_pre_direction <- point_data[2,-1]
  table_weight_pre <- data.frame(t(apply(table_pre,2,entropy_positive)))
  table_weight_pre[,which(table_pre_direction==1)] <- 
    apply(table_pre[,which(table_pre_direction==1)],2,entropy_negative)
  
  table_weight <- (1-table_weight_pre)/sum(1-table_weight_pre)
  
  # Vector normalization
  vector_normalize <- function(x){
    x / sqrt(sum(x^2))
  }
  
  table_normal <- data.frame(apply(table_pre,2,vector_normalize))
  
  # Matrix weighted normalization
  table_normal_weight <-  data.frame(t(apply(table_normal,1,
                                             function(x){x*unlist(table_weight)})))
  # Determine positive ideal solutions and negative ideal solutions
  best_case <- apply(table_normal_weight,2,max)
  worst_case <- apply(table_normal_weight,2,min)
  
  distance_best <- apply(table_normal_weight,1,
                         function(x){sqrt(sum((x -best_case)^2))})
  distance_worst <- apply(table_normal_weight,1,
                          function(x){sqrt(sum((x -worst_case)^2))})
  proximity_data <- data.frame(distance_best  = round(distance_best,3),
                               distance_worst = round(distance_worst,3))
  # Calculate the distance between each solution and the positive and negative ideal solutions
  proximity_data[,'proximity'] = round(distance_worst / (distance_best + distance_worst),3)
  proximity_data[,'rank'] = rank(-proximity_data[,'proximity'])
  result <- proximity_data
  result_1 <- round(table_weight,3)
  rownames(result_1)[1] <- "weight"
  result_2 <- round(table_normal_weight,3)
  return(list(result,result_1,result_2))
}

# PSA results of TOPSIS
compare_result_TOPSIS_graph <- function(point_data=file,interval_data=file_sen){
  table_pre_point <- point_data[-c(1,2),-1]
  table_pre_type <- interval_data[1,-1]
  table_pre_ci <- interval_data[-c(1,2),-1]
  
  table_pre_all <- array(0,dim=c(nrow(table_pre_point),ncol(table_pre_point),n_times))
  for(i in 1:nrow(table_pre_point)){
    for(j in 1:ncol(table_pre_point)){
      table_pre_all[i,j,] <- param_value(param_type=table_pre_type[j],param_point=as.numeric(table_pre_point[i,j]),
                                         param_low=as.numeric(table_pre_ci[i,j]),
                                         param_high=as.numeric(table_pre_ci[i,(j+ncol(table_pre_point))]),
                                         
                                         lower_bound=as.numeric(interval_data[2,j+1]),
                                         upper_bound=as.numeric(interval_data[2,(j+1+ncol(table_pre_point))]))
      
      
      }
  }
  
  
  table_score <- array(0,dim=c(nrow(table_pre_point),n_times))
  table_sort <- array(0,dim=c(nrow(table_pre_point),n_times))
  rownames(table_score) <- rownames(table_sort) <- unlist(point_data[-c(1,2),1])
  
  for(k in 1:n_times){
    table_pre <- table_pre_point
    table_pre[1:nrow(table_pre),1:ncol(table_pre)] <- table_pre_all[,,k]
    
    table_pre[table_pre$persistence==2,]$cost_effectiveness <- table_pre[table_pre$persistence==2,]$price/
      (table_pre[table_pre$persistence==2,]$mortality*table_pre[table_pre$persistence==2,]$effectiveness)/78
    table_pre[table_pre$persistence==1,]$cost_effectiveness <- table_pre[table_pre$persistence==1,]$price/
      (table_pre[table_pre$persistence==1,]$mortality*table_pre[table_pre$persistence==1,]$effectiveness)
    table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness <- table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness/2
    
    for(i in 1:ncol(table_pre)){
      table_pre[(table_pre[,i]==max(table_pre[,i])),i][1,] <-
        table_pre[(table_pre[,i]==max(table_pre[,i])),i][1,]*1.001
      table_pre[(table_pre[,i]==min(table_pre[,i])),i][1,] <-
        table_pre[(table_pre[,i]==min(table_pre[,i])),i][1,]*0.999
    }
    
    # entropy of positive indicator
    entropy_positive <- function(x){
      # Indicator normalization
      y = (x -min(x)) / (max(x) - min(x))
      # Indicator weight
      p = y / sum(y)
      # Indicator entropy
      entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
    }
    
    # entropy of negative indicator
    entropy_negative <- function(x){
      # Indicator normalization
      y = (max(x) -x) / (max(x) - min(x))
      # Indicator weight
      p = y/ sum(y)
      # Indicator entropy
      entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
    }
    
    table_pre_direction <- point_data[2,-1]
    table_weight_pre <- data.frame(t(apply(table_pre,2,entropy_positive)))
    table_weight_pre[,which(table_pre_direction==1)] <- 
      apply(table_pre[,which(table_pre_direction==1)],2,entropy_negative)
    
    table_weight <- (1-table_weight_pre)/sum(1-table_weight_pre)
    
    # Vector normalization
    vector_normalize <- function(x){
      x / sqrt(sum(x^2))
    }
    
    table_normal <- data.frame(apply(table_pre,2,vector_normalize))
    
    # Matrix weighted normalization
    table_normal_weight <-  data.frame(t(apply(table_normal,1,
                                               function(x){x*unlist(table_weight)})))
    # Determine positive ideal solutions and negative ideal solutions
    best_case <- apply(table_normal_weight,2,max)
    worst_case <- apply(table_normal_weight,2,min)
    # Determine positive ideal solutions and negative ideal solutions
    distance_best <- apply(table_normal_weight,1,
                           function(x){sqrt(sum((x -best_case)^2))})
    distance_worst <- apply(table_normal_weight,1,
                            function(x){sqrt(sum((x -worst_case)^2))})
    proximity_data <- data.frame(distance_best  = round(distance_best,3),
                                 distance_worst = round(distance_worst,3))
    # Calculate the distance between each solution and the positive and negative ideal solutions
    proximity_data[,'proximity'] = round(distance_worst / (distance_best + distance_worst),5)
    proximity_data[,'rank'] = rank(-proximity_data[,'proximity'])
    result <- proximity_data
    
    table_score[,k] <- result$proximity
    table_sort[,k] <- result$rank
  }
  
  table_score <- reshape2::melt(table_score, id = rownames(table_score))
  library(ggplot2)
  p1 <- ggplot(table_score, aes(`Var1`,`value`,group=`Var1`)) + geom_boxplot() + 
    labs(x="vaccine",y="MCDA score") + 
    theme_classic()
  
  table_sort_pre <- matrix(0,nrow(table_pre_point),nrow(table_pre_point))
  for(i in 1:ncol(table_sort_pre)){
    for(j in 1:nrow(table_sort_pre)){
      table_sort_pre[j,i] <- sum(table_sort[j,]==i)
    }
  }
  
  table_sort_pre <- data.frame(table_sort_pre)
  colnames(table_sort_pre)[1:6] <- 1:6
  table_sort_pre[,'vaccine'] <- point_data[-c(1,2),1]
  table_sort <- tidyr::gather(table_sort_pre,key="rank",value="proportion",-"vaccine")
  table_sort$proportion <- as.numeric(table_sort$proportion)/n_times
  library(ggplot2)
  p2 <- ggplot(table_sort,aes(x=`vaccine`,y=`proportion`,fill=`rank`)) + 
    scale_fill_brewer(palette = "Blues")+
    geom_bar(stat = 'identity',width = 0.5,colour="black") + 
    theme_classic()
  
  return(list(p1,p2,table_score,table_sort))
}

# point results of TOPSIS (Rank)
compare_result_TOPSIS_Rank <- function(point_data=file){
  table_pre <- point_data[-c(1,2),-1]
  table_pre[table_pre$persistence==2,]$cost_effectiveness <- table_pre[table_pre$persistence==2,]$price/
    (table_pre[table_pre$persistence==2,]$mortality*table_pre[table_pre$persistence==2,]$effectiveness)/78
  table_pre[table_pre$persistence==1,]$cost_effectiveness <- table_pre[table_pre$persistence==1,]$price/
    (table_pre[table_pre$persistence==1,]$mortality*table_pre[table_pre$persistence==1,]$effectiveness)
  table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness <- table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness/2
  
  # entropy of positive indicator
  entropy_positive <- function(x){
    # Indicator normalization
    y = (x -min(x)) / (max(x) - min(x))
    # Indicator weight
    p = y / sum(y)
    # Indicator entropy
    entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
  }
  
  # entropy of negative indicator
  entropy_negative <- function(x){
    # Indicator normalization
    y = (max(x) -x) / (max(x) - min(x))
    # Indicator weight
    p = y/ sum(y)
    # Indicator entropy
    entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
  }
  
  table_pre <- apply(table_pre,2,rank)
  
  for(i in 1:ncol(table_pre)){
    table_pre[(table_pre[,i]==max(table_pre[,i])),i][1] <-
      table_pre[(table_pre[,i]==max(table_pre[,i])),i][1]*1.001
    table_pre[(table_pre[,i]==min(table_pre[,i])),i][1] <-
      table_pre[(table_pre[,i]==min(table_pre[,i])),i][1]*0.999
  }
  
  table_pre_direction <- point_data[2,-1]
  table_weight_pre <- data.frame(t(apply(table_pre,2,entropy_positive)))
  table_weight_pre[,which(table_pre_direction==1)] <- 
    apply(table_pre[,which(table_pre_direction==1)],2,entropy_negative)
  
  table_weight <- (1-table_weight_pre)/sum(1-table_weight_pre)
  
  # Vector normalization
  vector_normalize <- function(x){
    x / sqrt(sum(x^2))
  }
  # Matrix weighted normalization
  table_normal <- data.frame(apply(table_pre,2,vector_normalize))
  
  table_normal_weight <-  data.frame(t(apply(table_normal,1,
                                             function(x){x*unlist(table_weight)})))
  # Determine positive ideal solutions and negative ideal solutions
  best_case <- apply(table_normal_weight,2,max)
  worst_case <- apply(table_normal_weight,2,min)
  # Calculate the distance between each solution and the positive and negative ideal solutions
  distance_best <- apply(table_normal_weight,1,
                         function(x){sqrt(sum((x -best_case)^2))})
  distance_worst <- apply(table_normal_weight,1,
                          function(x){sqrt(sum((x -worst_case)^2))})
  proximity_data <- data.frame(distance_best  = round(distance_best,3),
                               distance_worst = round(distance_worst,3))
  #计算相对接近度并进行排序
  proximity_data[,'proximity'] = round(distance_worst / (distance_best + distance_worst),3)
  proximity_data[,'rank'] = rank(-proximity_data[,'proximity'])
  result <- proximity_data
  result_1 <- round(table_weight,3)
  rownames(result_1)[1] <- "weight"
  result_2 <- round(table_normal_weight,3)
  return(list(result,result_1,result_2))
}

# PSA results of TOPSIS (Rank)
compare_result_TOPSIS_Rank_graph <- function(point_data=file,interval_data=file_sen){
  table_pre_point <- point_data[-c(1,2),-1]
  table_pre_type <- interval_data[1,-1]
  table_pre_ci <- interval_data[-c(1,2),-1]
  
  table_pre_all <- array(0,dim=c(nrow(table_pre_point),ncol(table_pre_point),n_times))
  for(i in 1:nrow(table_pre_point)){
    for(j in 1:ncol(table_pre_point)){
      table_pre_all[i,j,] <- param_value(param_type=table_pre_type[j],param_point=as.numeric(table_pre_point[i,j]),
                                         param_low=as.numeric(table_pre_ci[i,j]),
                                         param_high=as.numeric(table_pre_ci[i,(j+ncol(table_pre_point))]),
                                         
                                         lower_bound=as.numeric(interval_data[2,j+1]),
                                         upper_bound=as.numeric(interval_data[2,(j+1+ncol(table_pre_point))]))
    }
  }
  
  table_score <- array(0,dim=c(nrow(table_pre_point),n_times))
  table_sort <- array(0,dim=c(nrow(table_pre_point),n_times))
  rownames(table_score) <- rownames(table_sort) <- unlist(point_data[-c(1,2),1])
  
  for(k in 1:n_times){
    table_pre <- table_pre_point
    table_pre[1:nrow(table_pre),1:ncol(table_pre)] <- table_pre_all[,,k]
    
    table_pre[table_pre$persistence==2,]$cost_effectiveness <- table_pre[table_pre$persistence==2,]$price/
      (table_pre[table_pre$persistence==2,]$mortality*table_pre[table_pre$persistence==2,]$effectiveness)/78
    table_pre[table_pre$persistence==1,]$cost_effectiveness <- table_pre[table_pre$persistence==1,]$price/
      (table_pre[table_pre$persistence==1,]$mortality*table_pre[table_pre$persistence==1,]$effectiveness)
    table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness <- table_pre[rownames(table_pre)=="HPV",]$cost_effectiveness/2
    
    # entropy of positive indicator
    entropy_positive <- function(x){
      # Indicator normalization
      y = (x -min(x)) / (max(x) - min(x))
      # Indicator weight
      p = y / sum(y)
      # Indicator entropy
      entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
    }
    
    # entropy of negative indicator
    entropy_negative <- function(x){
      # Indicator normalization
      y = (max(x) -x) / (max(x) - min(x))
      # Indicator weight
      p = y/ sum(y)
      # Indicator entropy
      entropy = -1/log(length(x)) * sum(ifelse(p == 0,0,p *log(p)))
    }
    
    table_pre <- apply(table_pre,2,rank)
    for(i in 1:ncol(table_pre)){
      table_pre[(table_pre[,i]==max(table_pre[,i])),i][1] <-
        table_pre[(table_pre[,i]==max(table_pre[,i])),i][1]*1.001
      table_pre[(table_pre[,i]==min(table_pre[,i])),i][1] <-
        table_pre[(table_pre[,i]==min(table_pre[,i])),i][1]*0.999
    }
    
    table_pre_direction <- point_data[2,-1]
    table_weight_pre <- data.frame(t(apply(table_pre,2,entropy_positive)))
    table_weight_pre[,which(table_pre_direction==1)] <- 
      apply(table_pre[,which(table_pre_direction==1)],2,entropy_negative)
    
    table_weight <- (1-table_weight_pre)/sum(1-table_weight_pre)
    
    # Vector normalization
    vector_normalize <- function(x){
      x / sqrt(sum(x^2))
    }
    # Matrix weighted normalization
    table_normal <- data.frame(apply(table_pre,2,vector_normalize))
    
    table_normal_weight <-  data.frame(t(apply(table_normal,1,
                                               function(x){x*unlist(table_weight)})))
    # Determine positive ideal solutions and negative ideal solutions
    best_case <- apply(table_normal_weight,2,max)
    worst_case <- apply(table_normal_weight,2,min)
    # Calculate the distance between each solution and the positive and negative ideal solutions
    distance_best <- apply(table_normal_weight,1,
                           function(x){sqrt(sum((x -best_case)^2))})
    distance_worst <- apply(table_normal_weight,1,
                            function(x){sqrt(sum((x -worst_case)^2))})
    proximity_data <- data.frame(distance_best  = round(distance_best,3),
                                 distance_worst = round(distance_worst,3))
    
    proximity_data[,'proximity'] = round(distance_worst / (distance_best + distance_worst),5)
    proximity_data[,'rank'] = rank(-proximity_data[,'proximity'])
    result <- proximity_data
    
    table_score[,k] <- result$proximity
    table_sort[,k] <- result$rank
  }
  
  table_score <- reshape2::melt(table_score, id = rownames(table_score))
  library(ggplot2)
  p1 <- ggplot(table_score, aes(`Var1`,`value`,group=`Var1`)) + geom_boxplot() + 
    labs(x="vaccine",y="MCDA score") + 
    theme_classic()
  
  table_sort_pre <- matrix(0,nrow(table_pre_point),nrow(table_pre_point))
  for(i in 1:ncol(table_sort_pre)){
    for(j in 1:nrow(table_sort_pre)){
      table_sort_pre[j,i] <- sum(table_sort[j,]==i)
    }
  }
  table_sort_pre <- data.frame(table_sort_pre)
  colnames(table_sort_pre)[1:6] <- 1:6
  table_sort_pre[,'vaccine'] <- point_data[-c(1,2),1]
  table_sort <- tidyr::gather(table_sort_pre,key="rank",value="proportion",-"vaccine")
  table_sort$proportion <- as.numeric(table_sort$proportion)/n_times
  library(ggplot2)
  p2 <- ggplot(table_sort,aes(x=`vaccine`,y=`proportion`,fill=`rank`)) + 
    scale_fill_brewer(palette = "Blues")+
    geom_bar(stat = 'identity',width = 0.5,colour="black") + 
    theme_classic()
  
  return(list(p1,p2,table_score,table_sort))
}


n_times=1000

# #point value and CI of parameters
# output$table_baseline <- param_add_cost_effect()
# output$table_baseline_CI <- param_CI()
output <- list()
#results of AHP
output$result_AHP <- compare_result_AHP(point_data = file_D)
output$result_AHP_graph1 <- compare_result_AHP_graph(point_data = file_D,interval_data = file_sen_D)[[1]]
output$result_AHP_graph2 <- compare_result_AHP_graph(point_data = file_D,interval_data = file_sen_D)[[2]]
compare_result_TOPSIS_graph()[[1]]
#results of TOPSIS
output$result_TOPSIS <- compare_result_TOPSIS()[[1]]
output$result_TOPSIS_mid1 <- compare_result_TOPSIS()[[2]]
output$result_TOPSIS_mid2 <- compare_result_TOPSIS()[[3]]
output$result_TOPSIS_graph1 <- compare_result_TOPSIS_graph()[[1]]
output$result_TOPSIS_graph2 <- compare_result_TOPSIS_graph()[[2]]

#results of TOPSIS (Rank)
output$result_TOPSIS_Rank <- compare_result_TOPSIS_Rank()[[1]]
output$result_TOPSIS_Rank_mid1 <- compare_result_TOPSIS_Rank()[[2]]
output$result_TOPSIS_Rank_mid2 <- compare_result_TOPSIS_Rank()[[3]]
output$result_TOPSIS_Rank_graph1 <- compare_result_TOPSIS_Rank_graph()[[1]]
output$result_TOPSIS_Rank_graph2 <- compare_result_TOPSIS_Rank_graph()[[2]]

#results of Delphi
output$result_Delphi <- compare_result_AHP()
output$result_Delphi_graph1 <- compare_result_AHP_graph()[[1]]
output$result_Delphi_graph2 <- compare_result_AHP_graph()[[2]]

library(patchwork)
jpeg("Figure 1-Probabilistic sensitivity analysis results according to AHP method.jpg",width=5000,height=5000,res=1000)
output$result_AHP_graph1/output$result_AHP_graph2
dev.off()
jpeg("Figure 2-Probabilistic sensitivity analysis results according to Delphi method.jpg",width=5000,height=5000,res=1000)
output$result_Delphi_graph1/output$result_Delphi_graph2
dev.off()

jpeg("Figure 3-Probabilistic sensitivity analysis results according to TOPSIS method-1.jpg",width=5000,height=5000,res=1000)
output$result_TOPSIS_Rank_graph1/output$result_TOPSIS_Rank_graph2
dev.off()
jpeg("Figure 4-Probabilistic sensitivity analysis results according to TOPSIS method-2.jpg",width=5000,height=5000,res=1000)
output$result_TOPSIS_graph1/output$result_TOPSIS_graph2
dev.off()

###########################V4_R1
library(tidyverse)

n_times_pre = c(1:10,(2:10)*10)*100
result_for_covergence <- data.frame(
  n_times = rep(n_times_pre,each=24),
  Method = rep(c("AHP","TOPSIS","TOPSIS_Rank","Delphi"),each=6,
               length.out=24*length(n_times_pre)),
  Var1 = rep(c("HPV","RV","varicella","PCV","Hib","influenza"),
             length.out=24*length(n_times_pre)),
  mean = NA,
  lower_95 = NA,
  upper_95 = NA
)

for(i in 1:length(n_times_pre)){
  n_times = n_times_pre[i]
  output <- list()
  output$result_AHP_all <- compare_result_AHP_graph(point_data = file_D,
                                                    interval_data = file_sen_D)[[3]]  %>%
    group_by(Var1) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      lower_95 = quantile(value, 0.025, na.rm = TRUE),
      upper_95 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  output$result_TOPSIS_all <- compare_result_TOPSIS_graph()[[3]]  %>%
    group_by(Var1) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      lower_95 = quantile(value, 0.025, na.rm = TRUE),
      upper_95 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  output$result_TOPSIS_Rank_all <- compare_result_TOPSIS_Rank_graph()[[3]]  %>%
    group_by(Var1) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      lower_95 = quantile(value, 0.025, na.rm = TRUE),
      upper_95 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  output$result_Delphi_all <- compare_result_AHP_graph()[[3]]  %>%
    group_by(Var1) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      lower_95 = quantile(value, 0.025, na.rm = TRUE),
      upper_95 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  result_pre <- rbind(output$result_AHP_all,
                      output$result_TOPSIS_all,
                      output$result_TOPSIS_Rank_all,
                      output$result_Delphi_all)
  
  result_for_covergence[(24*(i-1)+1):(24*i),4:6] <- result_pre[,2:4]
  print(i)
}

write.csv(result_for_covergence,"result_for_covergence.csv")

result_for_covergence <- read.csv("result_for_covergence.csv",row.names = 1)
library(ggplot2)
library(dplyr)

result_for_covergence <- result_for_covergence %>%
  mutate(
    Var1 = as.factor(Var1),
    n_times = as.numeric(n_times)
  )

png("covergence.png",width = 3000,height=2000,res = 300)
ggplot(
  result_for_covergence,
  # result_for_covergence[result_for_covergence$n_times>2000,],
  aes(
    x = n_times,
    y = mean,
    group = Var1,
    color = Var1
  )
) +
  # geom_ribbon(
  #   aes(
  #     ymin = lower_95,
  #     ymax = upper_95,
  #     fill = Var1
  #   ),
  #   alpha = 0.15,
  #   colour = NA
  # ) +
geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ Method,
    scales = "free_y"
  ) +
  labs(
    x = "Number of Monte Carlo simulations",
    y = "Mean",
    color = "Vaccine",
    fill = "Vaccine"
  ) +
  scale_x_continuous(
    breaks = c(1:10)*1000
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
dev.off()


