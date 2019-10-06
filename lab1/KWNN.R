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

## Метод взвешенных k-ближайших соседей
kwNN <- function(xl, z, k, q)
{
  ## Сортируем выборку согласно классифицируемого	объекта
  orderedXl <- sortObjectsByDist(xl, z)
  n <- dim(orderedXl)[2]
  ## Получаем классы первых k соседей
  classes <- orderedXl[1:k, n]
  ## Получаем матрицу с весом каждого ближайшего элемента
  weighted_classes <- cbind(as.matrix(classes),q^as.numeric(1:k))
  
  lvls <- levels(classes)
  ## Суммируем все одинаковые классы и добавляем значения в матрицу
  result <- as.matrix(sapply(1:length(lvls), function(x) sum(as.numeric(weighted_classes[which(weighted_classes[,1] == lvls[x]),2]))))
  
  ## Находим класс, который доминирует среди первых k соседей
  class <- lvls[which.max(result)]	
  return (class)
}


colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
plot(iris[, 3:4], pch = 21, bg = colors[iris$Species], col = colors[iris$Species], asp = 1)

##set the 
z <- c(2.7,1)
xl <- iris[,3:5]

##kwNN
class <- kwNN(xl, z, k=6, q=0.5)
points(z[1], z[2], pch = 22, bg = colors[class], asp = 1)