###########################################################
### PMM algorithm

PMMsearchMet <- function(yHatMis, yHatObs) {
	posNear <- which(abs(yHatMis-yHatObs) == min(abs(yHatMis-yHatObs))) 
	posNear <- ifelse(length(posNear) > 1, sample(posNear,1), posNear)
	return(posNear)
}
