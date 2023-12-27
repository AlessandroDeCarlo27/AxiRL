get_epsilon<-function(episode){
  max_epsilon=1
  min_epsilon=0.3
  n_episode_saturation=25000
  
  
  
  
  lambda=-log(min_epsilon/max_epsilon)/(n_episode_saturation-1) #-1 perchè episode parte da 1 e non 0
  epsilon<-max(min_epsilon,max_epsilon*exp(-lambda*(episode-1)))
  
  return(epsilon)
}