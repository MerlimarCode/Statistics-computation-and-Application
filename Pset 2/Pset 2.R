
interSum <-0
tot2 <-0
totSq  <-0
A <- function(x) log(x+1)
B <- function(x) x*x

for (x in 1:21) {
  for (y in 1:21) {
    if (x != y & x<y){
      fileName = sprintf("allhic/chr%d_chr%d.txt",x,y);
      data <- read.delim(fileName)
      
      colnames(data)[1] = "X"
      colnames(data)[2] = "Y"
      colnames(data)[3] = "interFreq"
      maxX = 1+ (max(data$X)/250000)
      maxY = 1+ (max(data$Y)/250000)
      est = maxX*maxY

      inter = data.frame(lapply(data[3],A))
      interSq = data.frame(lapply(inter,B))
      total  = sum(inter$interFreq,na.rm=TRUE)
      tot2  = tot2 + est
      
      totSq = totSq + sum(interSq,na.rm=TRUE)
      interSum = interSum + total
    }
  }}
print(interSum/tot2)
print(totSq/tot2)
varience = (totSq/tot2 - (interSum/tot2)*(interSum/tot2))
sd = varience^.5
test <- read.csv("allhic/chr1_chr1.txt")

chr19 <- read.csv("allhic/chr1_chr1.txt")



