# Are healthcare policies, factors, and outcomes related to COVID-19 deaths?
# Author: Chris Walter - chriswalter.info
# Data sources: County health data - Robert Wood Johnson Foundation and 
#                                    the University of Wisconsin Population Health Institute
#                                    Retrieved from https://www.arcgis.com/home/item.html?id=4dc9f89d893b40cf85e4fd59e5a444f0
#                                    on July 8, 2020
#                                    Metadata: https://www.countyhealthrankings.org/explore-health-rankings/our-methods
#               Covid 19 data      - Aggregated from state-level data by USAfacts.org
#                                    Retrieved from https://usafacts.org/visualizations/coronavirus-covid-19-spread-map/
#                                    on July 8, 2020
#                                    Metadata: https://usafacts.org/articles/detailed-methodology-covid-19-data/


# load data, calculate covid death sum by county, clean and merge data
hlth <- read.csv("countyhealth.csv", header = TRUE)     # load county health data  
ronad <- read.csv("coviddeaths.csv", header = TRUE)     # load covid death data
ronac <- read.csv("covidcases.csv", header = TRUE)      # load covid case data
ronad$covdeaths <- ronad[, 172]                         # create total deaths column
ronac$covcases <- ronac[, 172]                          # total cases column
ronad <- ronad[, c(1, 173)]                             # remove un-needed columns
ronac <- ronac[, c(1, 173)]                             # remove un-needed columns
colnames(ronad)[1] <- "FIPS"                            # change colname for easy merging
colnames(ronac)[1] <- "FIPS"                            # change colname for easy merging
hlth <- merge(hlth, ronad, by = "FIPS", all.x = TRUE)   # merge ronad with hlth
hlth <- merge(hlth, ronac, by = "FIPS", all.x = TRUE)   # merge ronac with hlth
hlth <- hlth[!is.na(hlth$covdeaths), ]                  # remove NAs


# calculate covid deaths as percentage of cases
hlth$dthpcas <- 0                                       # initialize column       
for (i in 1:length(hlth[, 1])) {                        # death per case calculation
  if (hlth$covdeaths[i] > 0) {
    hlth$dthpcas[i] <- (hlth$covdeaths[i] / hlth$covcases[i]) * 100
  } else {
    hlth$dthpcas[i] <- 0
    }
}

# clean data and select variables for analysis
library(expss)                                          # load expss
hlthcas <- subset(hlth, hlth$covcases > 0)              # subset counties with cases
countieswithcases <- length(hlthcas[,1])                # num. of counties with at least 1 case
counts <- data.frame(                                   # initialize dataframe
                     "nam" = as.character(),
                     "count" = as.double()
)
for (i in 1: length(colnames(hlthcas))) {               # get count for each variable
  c <- data.frame(
                  "nam" = colnames(hlthcas)[i],
                  "count" = length(hlthcas[,1]) - count_if(NA, hlthcas[, i])
  )
  counts <- rbind(counts, c)
}
counts <- counts[-c(grep("numerator", counts$nam)), ]    # remove numerator variables
counts <- counts[-c(grep("denominator", counts$nam)), ]  # remove denominator variables
counts <- counts[-c(grep("rank", counts$nam)), ]         # remove county ranking variables
counts <- counts[-c(grep("Rank", counts$nam)), ]         # remove county Ranking variables
counts <- counts[-c(grep("Shape", counts$nam)), ]        # remove shape variables
counts <- counts[-c(grep("state", counts$nam)), ]        # remove extra state variable
counts <- counts[-c(grep("OBJECT", counts$nam)), ]       # remove OBJECTID variables
counts <- counts[-c(grep("Code", counts$nam)), ]         # remove fips code variables
counts <- counts[-which(counts$count < 3033), ]          # remove columns with < 3033 instances 
hlthcas <- hlthcas[, which(colnames(hlthcas) %in% counts$nam)]  # apply removal to dataframe
hlthcas <- hlthcas[complete.cases(hlthcas), ]            # remove incomplete cases (row with at least 1 NA)
countiesusedinanalysis <- length(hlthcas[,1])            # num. of counties used in analysis
(countiesusedinanalysis / countieswithcases) * 100       # ~96.4 % of counties with cases used in analysis
hlthcas <- hlthcas[order(hlthcas$FIPS), ]                # order data by fips
  # regarding which variables were selected for use: 
  # there is a natural break at 3033 - where state rankings are calculated
  # this appears to be the complete case data cutoff for the original county
  # health ranking analysis. cutting the complete cases here preserves ~96.4% of
  # the counties with covid cases for the analysis. Additionally, numerators
  # and denominators were removed, their dividend was included if it occurred in
  # 3033 or more cases. Finally, rankings were removed, as they are done by state,
  # and are not comparable across states. Furthermore, they are not very comparable
  # within states.

# remove population variables and additional numerators for nmds analysis,
# will save them to add them back to the final data frame
counts <- counts[-c(grep("Population", counts$nam)), ]   # remove Population variables
counts <- counts[-c(grep("population", counts$nam)), ]   # remove population variables
hlthnmds <- hlthcas[, which(colnames(hlthcas) %in% counts$nam)] # apply removal to dataframe
hlthnmds <- hlthnmds[, -c(1:4, 55:57)]                  # remove fips, area, id, and covid vars.

# calculate nmds for health variables
lapply(c("vegan", "scales"), library, character.only = TRUE)  # load vegan and scales
bray.hlth <- vegdist(hlthnmds, method = "bray")          # calculate bray-curtis distances
nmds.hlth <- metaMDS(bray.hlth)                          # calculate nmds scores
stressplot(nmds.hlth)                                    # look at shepard plot
plot(nmds.hlth, type = "t")                              # look at generic nmds plot
hlthcas <-cbind(hlthcas, scores(nmds.hlth ))             # add nmds scores to data
hlthcas$hlthidx <- 1 - rescale(hlthcas$NMDS1, to = c(0, 1))  # rescale NMDS to 0-1 index and flip scale
write.csv(hlthcas, "healthdata.csv", row.names = FALSE) # write data to file

