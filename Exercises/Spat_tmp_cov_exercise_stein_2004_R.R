---
title: "Exercise_12_15"
author: "Joonwon Lee"
date: "2023-12-16"
output:
  html_document:
    df_print: paged
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r libraries, include=FALSE}
library(fields)
library(GpGp)
library(MASS)
library(ggplot2)
library(gridExtra)
library(Matrix)         # is.positive.definite(cov_list[[2]])
library(pracma)
```


```{r}
matern_cov = function(v,r){
  out = 2^{1-v}/gamma(r) * r^v * besselK(r,v)
  return(out)
}
```


```{r, cache=TRUE}

my_sim1 = function(x_min=0,y_min=0,x_max=100,y_max=100, n=500, theta=3,range=5,cov_index=1){
  
  x <- runif(n, min = x_min, max = x_max)
  y <- runif(n, min = y_min, max = y_max) 
  
  s <- (outer(x,x,"-")) 
  t= (outer(y,y,"-"))
  
  # 1) separable model
  cov_mat = theta * exp(-abs(t)/range) * theta * exp(- (abs(s)/range))
  
  eps = 1e-8         
  abs_s = abs(s)  

  a = 0.5  
  c=0.3
  eps = 1e-8
  
  # 2) symmetric, smoother along the axes  MA 2003 (5.3)
  tmp1 = exp(-a*t) * erfc( sqrt(abs_s+ eps )- a*t/(2*sqrt(abs_s+ eps )    )  )
  tmp2 = exp( a*t) * erfc( sqrt(abs_s+ eps ) + a*t/(2*sqrt(abs_s+ eps )   )  ) # without eps, singular
  cov_mat2 = tmp1 + tmp2
  d= dim(cov_mat2)[1]
  
  eigen(cov_mat2)$values
  
  # symmetric, smoother along the axes  (6)
  a = 0.5  
  c=0.3
  eps = 1e-8
  
  r = abs(s) +eps

  tmp1 = pi^2/(16*c^6)*exp(a*r ) * (erfc( c*a* abs(t)^{1/2} + r/(2*c*abs(t)^{1/2}))  )
  tmp1_2 = 1/a^3 - r/a^2 + (4*c^4*t^2)/r
  tmp2 = pi^2/(16*c^6)*exp(-a*r) * (erfc( c*a* abs(t)^{1/2} - r/(2*c*abs(t)^{1/2}))  )
  tmp2_2 = 1/a^3 + r/a^2 - (4*c^4*t^2)/r
  tmp3 = pi^{3/2}*abs(t)^{1/2}/ (4*c^5*a^2)* exp(  -c^2*a^2*abs(t)- r^2/(4*c^2*abs(t)))
  cov_mat2 = tmp1*tmp1_2 + tmp2*tmp2_2 + tmp3

  eigen(cov_mat2)$values
  
  # lacking full symmetry
  
  a= 10
  beta_1 = 3
  beta_2 = 4; c1 = a/beta_1^2; c2= a/beta_2^2
  b=1
  var = 1/12* (x_max-x_min)^2
  
  eps = 1e-8
  abs_s = abs(s) +eps  # k(x,t) depends on(x,t) through (|x|,t)
  t = t + eps 
  v= 1.5; d=1; z=1  # v=1/2 expo v=3/2 (1+d/a)*exp(-d/a)
  yy= (beta_1^2*abs_s^2+beta_2^2*t^2)^{1/2}  # defined below (7) in page 6.
  tau = b* beta_1* beta_2/(2*a)       
  
  yy= yy/0.5
  tmp3 = var*(1+(yy))*exp(-(yy))
  # tmp3[n] = tmp3[n-1]
  tmp4 = var*exp(-(yy))
  # tmp4[n] = tmp4[n-1]
  cov_mat3 = a*(
            (2*v+d+1)*tmp3 - 2*tau*(beta_1*abs_s*z)*beta_2*t*tmp4
  )      
  eigen(cov_mat3)$values
  
  cov_list = list(cov_mat,cov_mat2,cov_mat3)

  ev = eigen(cov_mat)
  ev2 = eigen(cov_mat2)
  ev3 = eigen(cov_mat3)
  
  sim_data <- rmvnorm(n=1, mean = rep(0, length(n)), sigma= cov_list[[cov_index]])
  sim_data = t(sim_data)
  
  tmp1 = rep(0,n)
  tmp2 = rep(0,n)
  tmp3 = rep(0,n)
  
  for (i in 1:n){
  tmp1[i] = t(ev$vectors[,i]) %*%(sim_data) / ev$values[i]^{1/2} 
  tmp2[i] = t(ev2$vectors[,i]) %*%(sim_data) / ev2$values[i]^{1/2}
  tmp3[i] = t(ev3$vectors[,i]) %*%(sim_data) / ev3$values[i]^{1/2}
  tmp1[i] = log(tmp1[i]^2)
  tmp2[i] = log(tmp2[i]^2)
  tmp3[i] = log(tmp3[i]^2)
  }
  
  data = data.frame(x= c(1:n), y1 = tmp1, y2= tmp2, y3 = tmp3 )
  
  par(mfrow=c(1,3))
  ggplot(data, aes(x,y1)) +
    geom_point()+
    geom_smooth(method = "loess", span = 0.5)  +
    labs(
      x = "Index j",
      y = "Model: exp"
    ) +
    ylim(c(-10, 10))
  
    ggplot(data, aes(x,y2)) +
    geom_point()+
    geom_smooth(method = "loess", span = 0.5)  +
    labs(
      x = "Index j",
      y = "Model: exp sq"
    ) +
    ylim(c(-10, 10))
    
    ggplot(data, aes(x,y3)) +
    geom_point()+
    geom_smooth(method = "loess", span = 0.5)  +
    labs(
      x = "Index j",
      y = "Model: matern"
    ) +
    ylim(c(-10, 10))
    
    ggplot(data, aes(x = x)) +
    geom_point(aes(y = y1), color = "red", size=0.5, alpha = 0.5) +
    geom_smooth(aes(y = y1), method = "lm", color = "red", se = FALSE) +
  
    geom_point(aes(y = y2), color = "green", size=0.5, alpha = 0.5) +
    geom_smooth(aes(y = y2), method = "lm", color = "green", se = FALSE) +
    geom_point(aes(y = y3), color = "blue", size=0.5, alpha = 0.5) +
    geom_smooth(aes(y = y3), method = "lm", color = "blue", se = FALSE) +
    ylim(c(-10, 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Smoothed Lines for y1, y2, and y3")
  
}  
```

```{r, warnings=FALSE}
my_sim1(x_min=00,y_min=00,x_max=100,y_max=100, n=500, theta=3,range=5,cov_index=1)

```


# 2.
Observe index j vs $\log(\hat{\lambda}_j)$

```{r}

my_sim2 = function(x_min=0,y_min=0,x_max=100,y_max=100, n=500, theta=3,range=5){
  
   
  x <- runif(n, min = x_min, max = x_max)
  y <- runif(n, min = y_min, max = y_max) 
  
  s <- (outer(x,x,"-")) 
  t= (outer(y,y,"-"))
  
  # 1) separable model
  cov_mat = theta * exp(-abs(t)/range) * theta * exp(- (abs(s)/range))
  
  eps = 1e-8
    
  abs_s = abs(s)  
  a = 0.5  
  c=0.3
  eps = 1e-8

  # 2) symmetric, smoother along the axes  MA 2003 (5.3)
  tmp1 = exp(-a*t) * erfc( sqrt(abs_s+ eps )- a*t/(2*sqrt(abs_s+ eps )          )  )
  tmp2 = exp( a*t) * erfc( sqrt(abs_s+ eps ) + a*t/(2*sqrt(abs_s+ eps )         )  ) # without eps, singular
  cov_mat2 = tmp1 + tmp2
  d= dim(cov_mat2)[1]
  
  eigen(cov_mat2)$values
  
  # symmetric, smoother along the axes  (6)
  a = 0.5  
  c=0.3
  eps = 1e-8
  
  r = abs(s) +eps

  tmp1 = pi^2/(16*c^6)*exp(a*r ) * (erfc( c*a* abs(t)^{1/2} + r/(2*c*abs(t)^{1/2}))  )
  tmp1_2 = 1/a^3 - r/a^2 + (4*c^4*t^2)/r
  tmp2 = pi^2/(16*c^6)*exp(-a*r) * (erfc( c*a* abs(t)^{1/2} - r/(2*c*abs(t)^{1/2}))  )
  tmp2_2 = 1/a^3 + r/a^2 - (4*c^4*t^2)/r
  tmp3 = pi^{3/2}*abs(t)^{1/2}/ (4*c^5*a^2)* exp(  -c^2*a^2*abs(t)- r^2/(4*c^2*abs(t)))
  cov_mat2 = tmp1*tmp1_2 + tmp2*tmp2_2 + tmp3

  eigen(cov_mat2)$values
  
  # lacking full symmetry
  
  a = 10
  beta_1 = 3
  beta_2 = 4; c1 = a/beta_1^2; c2= a/beta_2^2
  b=1
  var = 1/12* (x_max-x_min)^2
  
  eps = 1e-8
  abs_s = abs(s) +eps  # k(x,t) depends on(x,t) through (|x|,t)
  t = t+ eps 
  v= 1.5; d=1; z=1  # v=1/2 expo v=3/2 (1+d/a)*exp(-d/a)
  yy= (beta_1^2*abs_s^2+beta_2^2*t^2)^{1/2}  # defined below (7) in page 6.
  tau = b* beta_1* beta_2/(2*a)       
  
  yy= yy/0.5
  tmp3 = var*(1+(yy))*exp(-(yy))
  # tmp3[n] = tmp3[n-1]
  tmp4 = var*exp(-(yy))
  # tmp4[n] = tmp4[n-1]
  cov_mat3 = a*(
            (2*v+d+1)*tmp3 - 2*tau*(beta_1*abs_s*z)*beta_2*t*tmp4
  )      
  eigen(cov_mat3)$values
  
  cov_list = list(cov_mat,cov_mat2,cov_mat3)
  
  ev = eigen(cov_mat)
  ev2 = eigen(cov_mat2)
  ev3 = eigen(cov_mat3)  
  
  data = data.frame(x=c(1:n), y1 =log(ev$values), y2 = log(ev2$values), y3 = log(ev3$values) )
  
  ggplot(data, aes(x = x)) +
  geom_point(aes(y = y1), color = "red", size=0.5, alpha = 0.5) +
  # geom_smooth(aes(y = y1), method = "lm", color = "red", se = FALSE) +

  geom_point(aes(y = y2), color = "green", size=0.5, alpha = 0.3) +
  # geom_smooth(aes(y = y2), method = "lm", color = "green", se = FALSE) +
  geom_point(aes(y = y3), color = "blue", size=0.5, alpha = 0.5) +
  # geom_smooth(aes(y = y3), method = "lm", color = "blue", se = FALSE) +
  ylim(c(-10, 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("Smoothed Lines for y1, y2, and y3")
}

```



```{r}
my_sim2(x_min=000,y_min=000,x_max=20,y_max=20, n=500, theta=3,range=6)
```



