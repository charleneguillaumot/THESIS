# jaccard.R
# Written in 2012 by Joona Lehtomäki <joona.lehtomaki@gmail.com>  

# To the extent possible under law, the author(s) have dedicated all 
# copyright and related and neighboring rights to this software to 
# the public domain worldwide. This software is distributed without any warranty. 

# You should have received a copy of the CC0 Public Domain Dedication along with 
# this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 

#' Calculate the Jaccard index between two Raster objects.
#'
#' The Jaccard coefficient measures similarity between sample sets, and is 
#' defined as the size of the intersection divided by the size of the union of 
#' the sample sets. Function assumes that values in rasters being compared
#' have values [0, 1] (other behaviour is not defined), such as rank rasters
#' produced by Zonation. 
#'
#' Specific parts of rasters can compared by using the threshhold argument.
#' For example, if threshhold is set to 0.9, the overlaps of pixels with 
#' values >= 0.9 are compared. In Zonation terms this translates into comparing
#' the best 10% of the landscape in both solutions.
#'
#' @param raster1 first raster to be compared
#' @param raster2 second raster to be compared
#' @param threshhold double value ([0, 1]) defining which fraction of the landscape is compared
#' @param warn.uneven boolean conrolling whether a waringn is raised if the fractions of landscape in different rasters are very uneven
#'
#' @return a double value defining the Jaccard index
#' @export
#' @references
#' @author Joona Lehtomaki \email{joona.lehtomaki@gmail.com}

jaccard <- function(raster1, raster2, threshhold, warn.uneven=TRUE) {
  
  # Get the values above the threshhold
  raster1.bin <- raster1 >= threshhold
  raster2.bin <- raster2 >= threshhold
  
  if (warn.uneven) {
    raster1.size <- raster::count(raster1.bin, 1)
    raster2.size <- raster::count(raster2.bin, 1)
    # Sort from smaller to larger
    sizes <- sort(c(raster1.size, raster2.size))
    if (sizes[2] / sizes[1] > 20) {
      warning("The extents of raster values above the threshhold differ more than 20-fold: Jaccard coefficient may not be informative.")
    }
  }
  
  # Calculate the intersection of the two rasters, this is given by adding 
  # the binary rasters together -> 2 indicates intersection
  combination <- raster1.bin + raster2.bin
  intersection <- combination == 2
  
  # Union is all the area covered by the both rasters
  union <- combination >= 1
  
  return(count(intersection, 1) / count(union, 1))
}
