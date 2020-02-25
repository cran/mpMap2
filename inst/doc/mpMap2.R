### R code from vignette source 'mpMap2.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: loadPackages
###################################################
library(mpMap2)
library(qtl)


###################################################
### code chunk number 2: pedigreeGraph
###################################################
p <- pedigree(lineNames = c("A", "B", "C", "F1-1", "F1-2", "F1-3", "F1-4", 
    "F2-1", "F2-2", "F3"), mother = c(0, 0, 0, 1, 1, 1, 2, 4, 6, 8), 
    father = c(0, 0, 0, 2, 2, 3, 3, 5, 7, 9), selfing = "finite")
plot(pedigreeToGraph(p))


###################################################
### code chunk number 3: simulation1
###################################################
#Generate map
map <- qtl::sim.map(len = rep(300, 2), n.mar = 301, anchor.tel = TRUE, 
    include.x = FALSE, eq.spacing = TRUE)
#Generate random funnels pedigree
pedigreeRF <- fourParentPedigreeRandomFunnels(initialPopulationSize = 1000, 
    nSeeds = 1, intercrossingGenerations = 1, selfingGenerations = 2)
#Analysis pedigreeRF will assume finite generations of selfing (two)
selfing(pedigreeRF) <- "finite"
#Prefix line names with RF
lineNames(pedigreeRF) <- paste0("RF", lineNames(pedigreeRF))
#Generate single funnel pedigree
pedigreeSF <- fourParentPedigreeSingleFunnel(initialPopulationSize = 1000, 
    nSeeds = 1, intercrossingGenerations = 1, selfingGenerations = 2)
#Analysis pedigreeSF will assume finite generations of selfing (two)
selfing(pedigreeSF) <- "finite"
#Prefix line names with SF
lineNames(pedigreeSF) <- paste0("SF", lineNames(pedigreeSF))
crossSingleFunnel <- simulateMPCross(map = map, pedigree = pedigreeSF, 
    mapFunction = haldane, seed = 1)
crossRandomFunnels <- simulateMPCross(map = map, pedigree = pedigreeRF,
    mapFunction = haldane, seed = 1)


###################################################
### code chunk number 4: combinedExperiments
###################################################
length(crossSingleFunnel@geneticData)
length(crossRandomFunnels@geneticData)
combined <- crossSingleFunnel + crossRandomFunnels
length(combined@geneticData)


###################################################
### code chunk number 5: generics1
###################################################
nMarkers(crossSingleFunnel)
nFounders(crossSingleFunnel)
nFounders(combined)
nLines(crossSingleFunnel)
nLines(combined)


###################################################
### code chunk number 6: summaryExample
###################################################
print(crossSingleFunnel)


###################################################
### code chunk number 7: fullyInformativeExample
###################################################
#Equivalent to crossSingleFunnel@geneticData[[1]]@founders[,1:5]
founders(crossSingleFunnel)[,1:5]


###################################################
### code chunk number 8: markerEncodingExample
###################################################
data <- rbind(rep(1, 5), rep(10, 5), rep(100, 5), rep(200, 5))
colnames(data) <- c("D1M1", "D1M2", "D1M3", "D1M4", "D1M5")
rownames(data) <- c("SFL1", "SLF2", "SLF3", "SLF4")
print(data)


###################################################
### code chunk number 9: hetDataExample
###################################################
#Equivalent to crossSingleFunnel@geneticData[[1]]@hetData[["D1M1"]]
hetData(crossSingleFunnel, "D1M1")


###################################################
### code chunk number 10: allelesDistribution
###################################################
table(finals(crossSingleFunnel)[,1])


###################################################
### code chunk number 11: lessInformative1
###################################################
combinedSNP <- combined + multiparentSNP(keepHets = TRUE)


###################################################
### code chunk number 12: lessInformative2
###################################################
combinedSNP <- combined
combinedSNP@geneticData[[1]] <- combinedSNP@geneticData[[1]] +
    multiparentSNP(keepHets = TRUE)
combinedSNP@geneticData[[2]] <- combinedSNP@geneticData[[2]] +
    multiparentSNP(keepHets = FALSE)
founders(combinedSNP@geneticData[[1]])[, 1:5]
hetData(combinedSNP, "D1M1") 


###################################################
### code chunk number 13: hetDataConstruction1
###################################################
nMarkers <- 10
hetData <- replicate(nMarkers, rbind(rep(0, 3), rep(1, 3)), simplify=FALSE)
names(hetData) <- paste0("M", 1:10)
hetData <- new("hetData", hetData)
hetData[[1]]


###################################################
### code chunk number 14: hetDataConstruction2
###################################################
nMarkers <- 10
hetData <- replicate(nMarkers, rbind(rep(0, 3), rep(1, 3), c(0, 1, 2), 
    c(1, 0, 2)), simplify=FALSE)
names(hetData) <- paste0("M", 1:10)
hetData <- new("hetData", hetData)
hetData[[1]]


###################################################
### code chunk number 15: construction1
###################################################
founders <- founders(crossSingleFunnel)
finals <- finals(crossSingleFunnel)
hetData <- hetData(crossSingleFunnel)
crossSingleFunnel <- mpcross(founders = founders, finals = finals, 
    pedigree = pedigreeSF, hetData = hetData)


###################################################
### code chunk number 16: errors1
###################################################
#Put in two errors for the founders
founders[3,3] <- 10
founders[2,2] <- NA
#Put in an error for the finals
finals[1, 1] <- 100
#Put in two errors for the hetData
hetData[4] <- list(rbind(rep(0, 3), rep(1, 3)))
hetData[[5]][1,1] <- NA

error <- try(crossSingleFunnel <- mpcross(founders = founders, finals = finals, 
    pedigree = pedigreeSF, hetData = hetData))
cat(error)


###################################################
### code chunk number 17: errors2
###################################################
errors <- listCodingErrors(founders = founders, finals = finals, hetData = hetData)
errors$invalidHetData
errors$null
head(errors$finals)


###################################################
### code chunk number 18: simulatedMap1
###################################################
simulatedMap <- qtl::sim.map(len = rep(100, 2), n.mar = 11, anchor.tel = TRUE, 
    include.x = FALSE, eq.spacing = FALSE)
#map object has class "map"
class(simulatedMap)
#Names of entries are chromosomoe names
names(simulatedMap)
#Markers are in increasing order. 
simulatedMap[["1"]]


###################################################
### code chunk number 19: flatExample1
###################################################
cbind(M1 = c("Founder 1" = 1, "Founder 2" = 0, "Founder 3" = 0, "Founder 4" = 1), M2 = c("Founder 1" = 0, "Founder 2" = 1, "Founder 3" = 0, "Founder 4" = 1))


###################################################
### code chunk number 20: approxFlat1
###################################################
cbind(M1 = c("Founder 1" = 1, "Founder 2" = 0, "Founder 3" = 0, "Founder 4" = 1, "Founder 5" = 1, "Founder 6" = 0, "Founder 7" = 1, "Founder 8" = 1), M2 = c("Founder 1" = 0, "Founder 2" = 0, "Founder 3" = 1, "Founder 4" = 0, "Founder 5" = 0, "Founder 6" = 0, "Founder 7" = 1, "Founder 8" = 1))


###################################################
### code chunk number 21: notFlat1
###################################################
cbind(M1 = c("Founder 1" = 1, "Founder 2" = 1, "Founder 3" = 0, "Founder 4" = 1, "Founder 5" = 1, "Founder 6" = 0, "Founder 7" = 1, "Founder 8" = 1), M2 = c("Founder 1" = 0, "Founder 2" = 1, "Founder 3" = 0, "Founder 4" = 0, "Founder 5" = 0, "Founder 6" = 0, "Founder 7" = 1, "Founder 8" = 1))


###################################################
### code chunk number 22: symmetricExample
###################################################
cbind(M1 = c("Founder 1" = 1, "Founder 2" = 1, "Founder 3" = 0, "Founder 4" = 1, "Founder 5" = 1, "Founder 6" = 0, "Founder 7" = 1, "Founder 8" = 1), M2 = c("Founder 1" = 0, "Founder 2" = 1, "Founder 3" = 1, "Founder 4" = 0, "Founder 5" = 0, "Founder 6" = 0, "Founder 7" = 1, "Founder 8" = 1))


###################################################
### code chunk number 23: rfExample
###################################################
rf <- estimateRF(object = combinedSNP, verbose = list(progressStyle = 1))


###################################################
### code chunk number 24: imputeExample1
###################################################
mappedSNP <- new("mpcrossMapped", combinedSNP, map = map)
imputed <- imputeFounders(mappedSNP, extraPositions = 
    list("2" = c("a" = 3.14, "b" = 66)))


###################################################
### code chunk number 25: imputeExample2
###################################################
imputed <- imputeFounders(mappedSNP, extraPositions = generateGridPositions(10))


###################################################
### code chunk number 26: imputeExample3
###################################################
imputed <- imputeFounders(mappedSNP)
imputationKey(imputed, experiment = 1)
table(imputationData(imputed, experiment = 1), finals(combined)[[1]])


###################################################
### code chunk number 27: imputeExample4
###################################################
table(imputationData(imputed, experiment = 2), finals(combined)[[2]])
imputed <- imputeFounders(mappedSNP, heterozygoteMissingProb = 1, 
    homozygoteMissingProb = 0.05)
table(imputationData(imputed, experiment = 2), finals(combined)[[2]])


###################################################
### code chunk number 28: fourParentExample1
###################################################
pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 800, 
    intercrossingGenerations = 0, selfingGenerations = 0, nSeeds = 1)
selfing(pedigree) <- "finite"
map <- qtl::sim.map(len = rep(300, 3), n.mar = 101, anchor.tel = TRUE, 
    include.x = FALSE, eq.spacing = FALSE)
cross <- simulateMPCross(pedigree = pedigree, map = map, 
    mapFunction = haldane, seed = 1)
crossSNP <- cross + multiparentSNP(keepHets=TRUE)
table(finals(cross))


###################################################
### code chunk number 29: fourParentExample2
###################################################
#Randomly rearrange markers
crossSNP <- subset(crossSNP, markers = sample(markers(cross)))
rf <- estimateRF(crossSNP, verbose = list(progressStyle = 1))


###################################################
### code chunk number 30: fourParentExample3
###################################################
grouped <- formGroups(rf, groups = 3, method = "average", clusterBy="theta")
try(omp_set_num_threads(1), silent = TRUE)
ordered <- orderCross(grouped, effortMultiplier = 2)
imputedTheta <- impute(ordered, verbose = list(progressStyle = 1)) 


###################################################
### code chunk number 31: fourParentExample4
###################################################
estimatedMap <- estimateMap(imputedTheta, maxOffset = 10)
estimatedMap <- jitterMap(estimatedMap)
#match up estimated chromosomes with original chromosomes
estChrFunc <- function(x) which.max(unlist(lapply(estimatedMap, 
    function(y) length(intersect(names(y), names(map[[x]]))))))
estimatedChromosomes <- sapply(1:3, estChrFunc)
tail(estimatedMap[[estimatedChromosomes[[1]]]])
tail(map[[1]])


###################################################
### code chunk number 32: fourParentExample5
###################################################
mappedObject <- new("mpcrossMapped", imputedTheta, map = estimatedMap)
imputedFounders <- imputeFounders(mappedObject)
summary <- table(imputedFounders@geneticData[[1]]@imputed@data[,markers(cross)], 
    finals(cross))
sum(diag(summary))/sum(summary)


