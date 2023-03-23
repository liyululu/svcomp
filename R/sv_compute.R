### svcomp
# This package computes the allocation of cost according to Shapley Value.
# The input can be matrix or tibble form contains the name of composition and volume.
# This package must use tidyverse and gtools, you would better install them first.

#' Shapley Value Computation of the transaction cost
#'
#' @description
#' This package computes the allocation of cost according to Shapley Value. The input can be matrix or tibble form contains the name of composition and volume. This package must use tidyverse and gtools, you would better install them first.
#'
#' @param CS_table A tibble, represents the cost structure.
#' @param table A matrix or a tibble.
#' @returns A list with 4 elements.
#'
#' @examples
#' library('tidyverse')
#' CS <- tibble('lower_v'=c(1,300001,3000001,20000001,100000001),'upper_v'=c(300000,3000000,20000000,100000000,Inf),'Commission_rate'=c(0.0035,0.0020,0.0015,0.0010,0.0005))
#' example1<-tibble('Composition'=c('1','2'),'v'=c(200000,250000))
#' sv_compute(CS,example1)


## Shapley Value Compute
sv_compute<-function(CS_table,table){
  ## Get the rate
  check_rate<-function(table,v){
    return (filter(table, v<= upper_v,v > lower_v)$Commission_rate)
  }
  ## Generate "1 & 2" name
  paste_pair<-function(pair){
    n<-length(pair)
    tmp<-as.character(pair[1])
    for(i in 2:n){
      tmp<-paste(tmp,"&",as.character(pair[i]))
    }
    tmp
  }
  table<-tibble(table)
  colnames(table)[1] <- 'Composition'
  colnames(table)[2] <- 'v'
  n<-nrow(table)
  ## Get the rate and commission as initial table
  table<-table %>% mutate('Commission_rate'=rep(0,n),'Commission'=rep(0,n))
  for(i in 1:n){
    table$Commission_rate[i]<-check_rate(CS_table,table$v[i])
    table$Commission[i]<-table$v[i]*table$Commission_rate[i]
  }
  initial_table<-table
  ## Get all combination
  for(j in 2:n){
    CB<-combn(n,j)
    for (i in 1:ncol(CB)){
      tmp<-CB[,i]
      v<-sum(table$v[tmp])
      cr<-check_rate(CS_table,v)
      table <- table %>%
        add_row('Composition'=paste_pair(tmp),
                'v'=v,
                'Commission_rate'=cr,
                'Commission'=cr*v)
    }
  }
  combination<-table
  ## Get Shapley value table
  ordername<-c('1st','2nd','3rd',paste(4:20, 'th', sep = ""))
  table.initial<-matrix(0,nrow=factorial(n),ncol=2*n+1)
  table.initial[,1]<-rep(1/factorial(n),factorial(n))
  table.initial[,2:(n+1)]<-permutations(n,n)
  table.initial<-as.data.frame(table.initial)
  names(table.initial)<-c('Probability',ordername[1:n],paste(1:n, 'marginal', sep = " "))
  for(i in 1:factorial(n)){
    tmp<-as.numeric(table.initial[i,])
    pr.tmp<-tmp[1:n+1]
    tmpx<-rep(0,n)
    tmpx[1]<-filter(table,Composition==as.character(pr.tmp[1]))$Commission
    table.initial[i,pr.tmp[1]+n+1]<-tmpx[1]
    for(j in 2:length(pr.tmp)){
      x.tmp<-pr.tmp[j]
      all.tmp<-pr.tmp[1:j]
      tmpx[j]<-filter(table,Composition==paste_pair(all.tmp[order(all.tmp)]))$Commission
      table.initial[i,x.tmp+n+1]<-tmpx[j]-tmpx[j-1]
    }
  }
  SV_table<-tibble(table.initial)
  ## Get final allocation according to Shapley value
  weight<-t(SV_table$Probability)%*%as.matrix(select(SV_table,`1 marginal`:paste(n, 'marginal')))
  ratio<-rep(0,n)
  ratio<-weight/sum(weight)
  allocation<-rep(0,n)
  total.amount<-combination$Commission[nrow(combination)]
  allocation<-as.data.frame(total.amount*ratio)
  names(allocation)<-initial_table$Composition[1:n]
  ## Get all results
  z<-list(initial_table=initial_table,combination=combination,SV_table=SV_table,allocation=allocation)
  return(z)
}

