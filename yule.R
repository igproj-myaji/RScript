library(languageR)
path <- '~/R/data/textmining/discussion/20100421'
filename <- 'hato-w.txt'
setwd(path)
temp.raw <- readLines(filename)
temp.list <- strsplit(temp.raw, split = "[[:blank:]]|[[:punct:]]")
temp.vec <- unlist(temp.list)
word.vec <- temp.vec[temp.vec != ""]
plot(growth.fnc(word.vec, size=200, nchunks=length(word.vec)%/%200))
word.spc <- spectrum.fnc(word.vec)
yule.fnc(word.spc)
