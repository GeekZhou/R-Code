#A comparison of R effciency
#给出了三个解答，slowwalk是运行最慢的，
#在我的电脑上大约要运算8秒，
#这段代码慢的原因在于用c来补进新的向量元素。
#而在R中，对一个向量进行新增的操作实际上是进行了复制重建，所以速度是很慢的。
#为了改进这一点，speedwalk1函数中，先用一个空向量占好了位置，再向这个向量中放东西，
#这样不需要反复新建，速度会增加10倍，只花了0.7秒左右。但R最重要的运算思路是向量化，
#也就是同时对向量各元素进行计算，而不是把向量元素一个个放到循环中。
#speedwalk2函数就是向量化计算的实现，这样下来代码速度又加快了十几倍，
#大约只需0.05秒，而且还没有用到多核心计算。
slowwalk <- function(n) {
  dir <- c(TRUE, FALSE) ## N/S or W/E
  step <- c(1, -1) ## N or S | E or W
  ## x is W/E and y is N/S
  x <- y <- 0 ## results
  .x <- .y <- 0 ## current position
  for (i in 1:n) {
    .dir <- sample(dir, 1)
    .step <- sample(step, 1)
    if (.dir) { ## N/S
      x <- c(x, .x) ## W/E unchanged
      .y <- .y + .step 
      y <- c(y, .y)
    } else { ## W/E
      y <- c(y, .y) ## N/S unchanged
      .x <- .x + .step 
      x <- c(x, .x)      
    }
  }
  return(list(x = x, y = y))
}
 
speedwalk1 <- function(n) {
  res <- matrix(rep(0,2*n),ncol=2)
  colnames(res) <- c('x','y')
  dir <- c(TRUE, FALSE) ## N/S or W/E
  step <- c(1, -1) ## N or S | E or W
  for (i in 2:n) {
    .dir <- sample(dir, 1)
    .step <- sample(step, 1)
    if (.dir) { ## N/S
      res[i,'x'] <- res[i-1,'x'] ## W/E unchanged
      res[i,'y'] <- res[i-1,'y'] + .step 
    } else { ## W/E
      res[i,'y'] <- res[i-1,'y'] ## N/S unchanged
      res[i,'x'] <- res[i-1,'x'] + .step 
    }
  }
  return(res)
}
 
speedwalk2 <- function(n) {
  NS <- sample(c(0,1),size=n,replace=T)
  WE <- 1 - NS
  .x <- sample(c(1,-1),size=n,replace=T)
  .x <- .x * NS
  .y <- sample(c(1,-1),size=n,replace=T)
  .y <- .y * WE
  x <- cumsum(.x)
  y <- cumsum(.y)
  res <- cbind(x,y)
  return(res)
}
 
system.time({
  trajectory <- slowwalk(2^15)
})
 
system.time({
  trajectory <- speedwalk1(2^15)
})
 
system.time({
  trajectory <- speedwalk2(2^15)
})
 
plot(trajectory,type='l')
