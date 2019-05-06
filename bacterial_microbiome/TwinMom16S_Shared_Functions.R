# Shared Functions for Twin Mom 16S Analyses

#---------------------------------#

# S3 function for extracting eigenvalues from ordination objects
#   Taken directly from 
#   https://github.com/joey711/phyloseq/blob/0a338569162922fb1841c1203a94384ce21af9b0/R/plot-methods.R
extract_eigenvalue = function(ordination) UseMethod("extract_eigenvalue", ordination)
# Default is to return NULL (e.g. for NMDS, or non-supported ordinations/classes).
extract_eigenvalue.default = function(ordination) NULL
# for pcoa objects
extract_eigenvalue.pcoa = function(ordination) ordination$values$Relative_eig
# for CCA objects
extract_eigenvalue.cca = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for RDA objects
extract_eigenvalue.rda = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for dpcoa objects
extract_eigenvalue.dpcoa = function(ordination) ordination$eig
# for decorana (dca) objects
extract_eigenvalue.decorana = function(ordination) ordination$evals

#---------------------------------#

GetOrdEigenLabs <- function(ordinationObj) {
  
  # Code taken directly from
  # https://github.com/joey711/phyloseq/blob/0a338569162922fb1841c1203a94384ce21af9b0/R/plot-methods.R
  
  eigvec = extract_eigenvalue(ordination = ordinationObj)
  # Fraction variability, fracvar
  fracvar = eigvec[axes = c(1,2)] / sum(eigvec)
  # Percent variability, percvar
  percvar = round(100*fracvar, 1)
  # The string to add to each axis label, strivar
  # Start with the curent axis labels in the plot
  strivar = c("Axis.1", "Axis.2")
  # paste the percent variability string at the end
  strivar = paste0(strivar, "   [", percvar, "%]")
  
  return(strivar)
  
}

#---------------------------------#
