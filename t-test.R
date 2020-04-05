require(data.table)

#Import the survival data
#Percentages of the colonies formed 

survival = fread("survival.txt")

#Dataframe consists of 3 columns, cell line, percentages of colonies formed and dose rates

#   head(survival, n=5)
#           values cline rate
#   1: 100.0000000 lemur  0.0
#   2:  56.5292096 lemur  2.5
#   3:  15.4639175 lemur  5.0
#   4:   4.4673540 lemur  7.5
#   5:   0.6872852 lemur 10.0

#Perform t-test
for (i in c(0.0,2.5,5.0,7.5,10.0)) {
    
    x = survival[which(rate == i & cline == "human"),]$values
    y = survival[which(rate == i & cline == "lemur"),]$values
    print(t.test(x,y))
    
}

#Import Immunoslot blot repair assay data for (6-4)PP
data_64 = fread("data_64.txt")
#Similarly, Immunoslot blot has the same columns.
#Perform t-test
for (i in c(0.0,2.5,5.0,7.5,10.0)) {
    
    x = data_64[which(rate == i & cline == "human"),]$values
    y = data_64[which(rate == i & cline == "lemur"),]$values
    print(t.test(x,y))
    
}

#Import Immunoslot blot repair assay data for CPD
data_cpd = fread("data_cpd.txt")
#Perform t-test
for (i in c(0.0,2.5,5.0,7.5,10.0)) {
    
    x = data_cpd[which(rate == i & cline == "human"),]$values
    y = data_cpd[which(rate == i & cline == "lemur"),]$values
    print(t.test(x,y))
    
}






