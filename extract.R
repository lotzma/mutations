## R functions for preparing the data set for analysis

#library(deconstructSigs)
library(robust)

# Extract exposure of mutation catalog M to signatures P using non-negative least squares

# Function to extract signatures from mutation matrix

# Load data
load('data/tcga_samples.Rdata')

sigdata.tcga <- data.tcga.new

# Compute averages in separate data frame
ages <- sort(unique(sigdata.tcga$Age))
averages <- data.frame(Age=ages)
averages$Sig1Mean <- NA
averages$Sig1Med <- NA
averages$Sig7Mean <- NA
averages$Sig7Med <- NA
for (a in ages) {
    i <- match(a,averages$Age)
    averages[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a),"Sig1Total"])
    averages[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a),"Sig1Total"])
    averages[i,c("Sig7Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a),"Sig7Total"])
    averages[i,c("Sig7Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a),"Sig7Total"])
    if (nrow(subset(sigdata.tcga, Age==a))>0) {
      averages[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a), family=poisson))["(Intercept)"])
      averages[i,c("Sig7RobMean")] <- exp(coef(glmRob(as.integer(Sig7Total)~1, data=subset(sigdata.tcga, Age==a), family=poisson))["(Intercept)"])
    }
}
   
# Compute averages in separate data frame (BRAF)
ages <- sort(unique(subset(sigdata.tcga, Cohort=="BRAF")$Age))
averages.braf <- data.frame(Age=ages)
averages.braf$Sig1Mean <- NA
averages.braf$Sig1Med <- NA
averages.braf$Sig1RobMean <- NA
averages.braf$Sig7Mean <- NA
averages.braf$Sig7RobMean <- NA
averages.braf$Sig7Med <- NA

for (a in ages) {
    i <- match(a,averages.braf$Age)
    averages.braf[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="BRAF"),"Sig1Total"])
    averages.braf[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="BRAF"),"Sig1Total"])
    averages.braf[i,c("Sig7Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="BRAF"),"Sig7Total"])
    averages.braf[i,c("Sig7Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="BRAF"),"Sig7Total"])
    if (nrow(subset(sigdata.tcga, Age==a & Cohort=="BRAF"))>0) {
      averages.braf[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="BRAF"), family=poisson))["(Intercept)"])
      averages.braf[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig7Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="BRAF"), family=poisson))["(Intercept)"])
    }
}

# Compute averages in separate data frame (NRAS)
ages <- sort(unique(subset(sigdata.tcga, Cohort=="NRAS")$Age))
averages.nras <- data.frame(Age=ages)
averages.nras$Sig1Mean <- NA
averages.nras$Sig1Med <- NA
averages.nras$Sig1RobMean <- NA
averages.nras$Sig7Mean <- NA
averages.nras$Sig7RobMean <- NA
averages.nras$Sig7Med <- NA
for (a in ages) {
    i <- match(a,averages.nras$Age)
    averages.nras[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NRAS"),"Sig1Total"])
    averages.nras[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NRAS"),"Sig1Total"])
    averages.nras[i,c("Sig7Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NRAS"),"Sig7Total"])
    averages.nras[i,c("Sig7Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NRAS"),"Sig7Total"])
    if (nrow(subset(sigdata.tcga, Age==a & Cohort=="NRAS"))>0) {
      averages.nras[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="NRAS"), family=poisson))["(Intercept)"])
      averages.nras[i,c("Sig7RobMean")] <- exp(coef(glmRob(as.integer(Sig7Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="NRAS"), family=poisson))["(Intercept)"])
    }
}

# Compute averages in separate data frame (NF1)
ages <- sort(unique(subset(sigdata.tcga, Cohort=="NF1")$Age))
averages.nf1 <- data.frame(Age=ages)
averages.nf1$Sig1Mean <- NA
averages.nf1$Sig1Med <- NA
averages.nf1$Sig1RobMean <- NA
averages.nf1$Sig7Mean <- NA
averages.nf1$Sig7RobMean <- NA
averages.nf1$Sig7Med <- NA
for (a in ages) {
  i <- match(a,averages.nf1$Age)
  averages.nf1[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NF1"),"Sig1Total"])
  averages.nf1[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NF1"),"Sig1Total"])
  averages.nf1[i,c("Sig7Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NF1"),"Sig7Total"])
  averages.nf1[i,c("Sig7Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NF1"),"Sig7Total"])
  if (nrow(subset(sigdata.tcga, Age==a & Cohort=="NF1"))>0) {
    averages.nf1[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="NF1"), family=poisson))["(Intercept)"])
    averages.nf1[i,c("Sig7RobMean")] <- exp(coef(glmRob(as.integer(Sig7Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="NF1"), family=poisson))["(Intercept)"])
  }
}

# Compute averages in separate data frame (W3)
ages <- sort(unique(subset(sigdata.tcga, Cohort=="W3")$Age))
averages.w3 <- data.frame(Age=ages)
averages.w3$Sig1Mean <- NA
averages.w3$Sig1Med <- NA
averages.w3$Sig1RobMean <- NA
averages.w3$Sig7Mean <- NA
averages.w3$Sig7RobMean <- NA
averages.w3$Sig7Med <- NA
for (a in ages) {
  i <- match(a,averages.w3$Age)
  averages.w3[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="W3"),"Sig1Total"])
  averages.w3[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="W3"),"Sig1Total"])
  averages.w3[i,c("Sig7Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="W3"),"Sig7Total"])
  averages.w3[i,c("Sig7Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="W3"),"Sig7Total"])
  if (nrow(subset(sigdata.tcga, Age==a & Cohort=="W3"))>0) {
    averages.w3[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="W3"), family=poisson))["(Intercept)"])
    averages.w3[i,c("Sig7RobMean")] <- exp(coef(glmRob(as.integer(Sig7Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="W3"), family=poisson))["(Intercept)"])
  }
}

# Compute averages in separate data frame (MALE)
ages <- sort(unique(subset(sigdata.tcga, Gender=="MALE")$Age))
averages.male <- data.frame(Age=ages)
averages.male$Sig1Mean <- NA
averages.male$Sig1Med <- NA
averages.male$Sig1RobMean <- NA
averages.male$Sig7Mean <- NA
averages.male$Sig7RobMean <- NA
averages.male$Sig7Med <- NA
for (a in ages) {
  i <- match(a,averages.male$Age)
  averages.male[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="MALE"),"Sig1Total"])
  averages.male[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="MALE"),"Sig1Total"])
  averages.male[i,c("Sig7Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="MALE"),"Sig7Total"])
  averages.male[i,c("Sig7Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="MALE"),"Sig7Total"])
  if (nrow(subset(sigdata.tcga, Age==a & Gender=="MALE"))>0) {
    averages.male[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Gender=="MALE"), family=poisson))["(Intercept)"])
    averages.male[i,c("Sig7RobMean")] <- exp(coef(glmRob(as.integer(Sig7Total)~1, data=subset(sigdata.tcga, Age==a & Gender=="MALE"), family=poisson))["(Intercept)"])
  }
}

# Compute averages in separate data frame (FEMALE)
ages <- sort(unique(subset(sigdata.tcga, Gender=="FEMALE")$Age))
averages.female <- data.frame(Age=ages)
averages.female$Sig1Mean <- NA
averages.female$Sig1Med <- NA
averages.female$Sig1RobMean <- NA
averages.female$Sig7Mean <- NA
averages.female$Sig7RobMean <- NA
averages.female$Sig7Med <- NA
for (a in ages) {
  i <- match(a,averages.female$Age)
  averages.female[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="FEMALE"),"Sig1Total"])
  averages.female[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="FEMALE"),"Sig1Total"])
  averages.female[i,c("Sig7Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="FEMALE"),"Sig7Total"])
  averages.female[i,c("Sig7Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Gender=="FEMALE"),"Sig7Total"])
  if (nrow(subset(sigdata.tcga, Age==a & Gender=="FEMALE"))>0) {
    averages.female[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Gender=="FEMALE"), family=poisson))["(Intercept)"])
    averages.female[i,c("Sig7RobMean")] <- exp(coef(glmRob(as.integer(Sig7Total)~1, data=subset(sigdata.tcga, Age==a & Gender=="FEMALE"), family=poisson))["(Intercept)"])
  }
}

averages.braf$Cohort <- "BRAF"
averages.nras$Cohort <- "NRAS"
averages.nf1$Cohort <- "NF1"
averages.w3$Cohort <- "W3"
averages.brafnras <- rbind(averages.braf, averages.nras)
averages.all <- rbind(averages.braf, averages.nras, averages.nf1, averages.w3)
averages.male$Gender <- "MALE"
averages.female$Gender <- "FEMALE"
averages.gender <- rbind(averages.male, averages.female)
    
# Save into file
write.csv(sigdata.tcga, 'data/tcga_4class.csv', sep=',')
save(sigdata.tcga, file = 'data/tcga_4class.Rdata')
write.table(averages, 'data/averages.csv', sep=',')
write.table(averages.braf, 'data/averages_braf.csv', sep=',')
write.table(averages.nras, 'data/averages_nras.csv', sep=',')
write.table(averages.nf1, 'data/averages_nf1.csv', sep=',')
write.table(averages.w3, 'data/averages_w3.csv', sep=',')
write.table(averages.brafnras, 'data/averages_brafnras.csv', sep=',')
write.table(averages.all, 'data/averages_all.csv', sep=',')
write.table(averages.male, 'data/averages_male.csv', sep=',')
write.table(averages.female, 'data/averages_female.csv', sep=',')
write.table(averages.gender, 'data/averages_gender.csv', sep=',')