chr19_20 <- read.delim("data/allhic/chr19_chr20.txt")
og <- read.delim("allhic/chr19_chr20.txt")

chr19_20[is.na(chr19_20)] = 0
colnames(chr19_20)[1] = "ROW"
colnames(chr19_20)[2] = "COLUMN"
colnames(chr19_20)[3] = "VALUE"
C <- function(x) as.integer((x/250000))+1
A <- function(x) log(x+1)

chr19_20[1] <- lapply(chr19_20[1], C)
chr19_20[2] <- lapply(chr19_20[2], C)
chr19_20[3] <- lapply(chr19_20[3], A)

library(plotly)
test1 <- max(chr19_20$ROW)
new_mat <- matrix(0, ncol = max(chr19_20$COLUMN), nrow = max(chr19_20$ROW))

for(i in 1:nrow(chr19_20)){

  new_mat[chr19_20[i,]$ROW, chr19_20[i,]$COLUMN] <- chr19_20[i,]$VALUE
}

plot_ly(z = new_mat, type = "heatmap",yaxis = list(title = 'Sepal Width (cm)'), xaxis = list(title = 'Sepal Width (cm)'))

fig2 <-  plot_ly(z= new_mat ,type = 'heatmap', mode = 'markers')%>%
  layout(title = 'HeatMap', plot_bgcolor = "#e5ecf6", xaxis = list(title = 'Chr 19'), 
         yaxis = list(title = 'Chr 20'), legend = list(title=list(text='<b> Species of Iris </b>')))
fig2
