library(ggplot2)
library(resample)
library(Rmisc)
library(dplyr)

ratings <- read.csv("dados/ml-latest-small/ratings.csv")
movies <- read.csv("dados/ml-latest-small/movies.csv")

ratings <- read.csv("dados/ml-latest/ratings.csv")
movies <- read.csv("dados/ml-latest/movies.csv")


ratings = full_join(ratings, movies)
rf = filter(ratings, grepl("Horror", ratings$genres))

ratings$i = 1


ratings2 = ratings %>% 
  group_by(movieId, title, genres) %>% 
  summarise(rating = mean(rating), n = sum(i)) %>% 
  filter(n >= 100)

set.seed(123)

sample_size = 20
get_confidence = function(sample_size, df){
  real.mean = mean(df$rating)
  runs = 1000
  contained = 0
  contained_clt = 0
  for(i in seq(1, runs)){
    #message(i)
    a_sample <- sample(rf$rating, sample_size)
    
    b = bootstrap(a_sample, mean, R = 1000)
    mean.news = CI.bca(b, probs = c(.025, .975))
    if(mean.news[1] <= real.mean & mean.news[2] >= real.mean){
      contained = contained+1
    }
    
    cltci = CI(a_sample)
    if(cltci[3] <= real.mean & cltci[1] >= real.mean){
      contained_clt = contained_clt+1
    }
  }
  return(c(bootstrap = contained / runs, clt = contained_clt / runs))
}

#get_confidence(100, rf)
get_confidence(20, rf)
get_confidence(10, rf)

ratings$length = floor(nchar(as.character(ratings$genres)) / 10)
ratings$length = as.factor(as.integer(ratings$length))

ggplot(ratings %>% filter(length %in% 1:5), aes(x = length, y = rating)) + 
  stat_summary(fun.y = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2)
