##########################
#### Score functions #####
##########################

#F1 micro
f1_micro <- function(truth,pred){
  
  #Aggregate precision and recall
  prec <- sum(pred==truth) / sum(pred > 0)
  rec <- sum(pred==truth) / length(truth > 0) #Ignore the unclassifiable/negative cells in the ground truth data set
  
  #F1 micro
  score <- 2 * prec * rec / (prec + rec)
  return(score)
}

#F1 macro
f1_macro <- function(truth,pred,bars){
  
  #Sample-specific precision and recall
  precs <- c()
  recs <- c()
  for(i in 1:length(bars)){
    bar_i <- bars[i]
    precs[i] <- sum(pred==bar_i & truth==bar_i) / sum(pred==bar_i)
    recs[i] <- sum(pred==bar_i & truth==bar_i) / sum(truth==bar_i)
  }
  
  #F1 macro
  scores <- 2 * precs * recs / (precs + recs)
  scores[is.na(scores)] <- 0 #If no predictions are made, penalize by setting F1 to 0
  score <- mean(scores)
  return(score)
}


#F1 weighted
f1_weighted <- function(truth,pred,bars,weights){
  
  #Sample-specific precision and recall
  precs <- c()
  recs <- c()
  for(i in 1:length(bars)){
    bar_i <- bars[i]
    precs[i] <- sum(pred==bar_i & truth==bar_i) / sum(pred==bar_i)
    recs[i] <- sum(pred==bar_i & truth==bar_i) / sum(truth==bar_i)
  }
  
  #F1 weighted
  scores <- 2 * precs * recs / (precs + recs)
  scores[is.na(scores)] <- 0 #If no predictions are made, penalize by setting F1 to 0
  score <- sum(weights * scores)
  return(score)
}

#Matthews Correlation Coefficient
mcc <- function(truth,pred,bars){
  
  #Ensure equal labels between truth and pred
  classes <- unique(c(0,bars))
  truth <- factor(truth, levels = classes)
  pred  <- factor(pred, levels = classes)
  
  #Confusion matrix
  conf <- table(truth, pred)
  
  #Number of calls in ground truth
  t_k <- rowSums(conf)
  
  #Number of calls by method
  p_k <- colSums(conf)
  
  #Number of correctly assigned cells
  c <- sum(diag(conf))
  
  #Total number of cells
  s <- sum(conf)
  
  #Normalization factor
  denom <- sqrt(
    (s^2 - sum(p_k^2)) *
      (s^2 - sum(t_k^2))
  )
  
  if (denom == 0) {
    return(NA)
  }
  
  #MCC
  mcc <- (c * s - sum(t_k * p_k)) / denom
  
  return(mcc)
  
}

#Concordance
concordance <- function(truth,pred){
  con <- sum(truth==pred) / length(truth)
  return(con)
  
}


#F1 scores
f1_scores <- function(truth,pred,bars){
  
  #Sample-specific precision and recall
  precs <- c()
  recs <- c()
  for(i in 1:length(bars)){
    bar_i <- bars[i]
    precs[i] <- sum(pred==bar_i & truth==bar_i) / sum(pred==bar_i)
    recs[i] <- sum(pred==bar_i & truth==bar_i) / sum(truth==bar_i)
  }
  
  #F1 macro
  scores <- 2 * precs * recs / (precs + recs)
  scores[is.na(scores)] <- 0 #If no predictions are made, penalize by setting F1 to 0
  return(scores)
}

