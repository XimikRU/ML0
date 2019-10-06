## Евклидово расстояние
euclideanDistance <- function(u, v)
{
  sqrt(sum((u - v)^2))
}

## Сортируем объекты согласно расстояния до объекта z
sortObjectsByDist <- function(xl, z, metricFunction = euclideanDistance)
{
  l <- dim(xl)[1]
  n <- dim(xl)[2] - 1
  ## Создаём матрицу расстояний
  distances <- matrix(NA, l, 2)
  for (i in 1:l)
  {
    distances[i, ] <- c(i, metricFunction(xl[i, 1:n], z))
  }
  ## Сортируем
  orderedXl <- xl[order(distances[, 2]), ]
  return (orderedXl);
}

## Применяем метод kNN
kNN <- function(xl, z, k = c(6))
{
  ## Сортируем выборку согласно классифицируемого	объекта
  orderedXl <- sortObjectsByDist(xl, z)
  n <- dim(orderedXl)[2] - 1
  ## Получаем классы первых k соседей
  classes <- orderedXl[1:k[length(k)], n + 1]
  
  class <- c(seq(length(k)))
  
  for (i in seq(length(k)))
  {
    ## Составляем таблицу(+считает общее число элементов каждого "типа") встречаемости каждого класса
    counts <- table(orderedXl[1:i,n+1])
    ## Находим класс, который чаще всего встречается среди первых k соседей
    class[i] <- names(which.max(counts))
  }
  
  return (class)
}

looCV_kNN <- function(dat, algo) {
  l <- length(dat[[1]])
  correct <- rep(0, length(l))
  
  for (i in seq(l)) {
    print(i)
    trainData <- dat[-i, ]
    control <- dat[i, ]
    
    res <- algo(trainData, control[1:2], seq(l))
    correct <- correct + (res != control$Species)
  }
  
  res <- (correct/l)
  return (res)
}

demo <- function()
{
  dat <- iris[,3:5]
  res <- looCV_kNN(dat,kNN)
  
  dataRes <- as.data.frame(res)
  
  plot(data.frame("k"=seq(length(dat[[1]])), "LOO(k)"=res), type="l")
  
  minK <- which.min(dataRes$res)
  points(minK,dataRes$res[minK],type="p",pch = 16, col = "green",bg = "green")
}
demo()