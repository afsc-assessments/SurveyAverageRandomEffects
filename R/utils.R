#' Read output file from the random effects (RE) model into a
#' data.frame
#'
#' @param file A character string of the file to be read in,
#'   defaulting to 'rwout.rep' in the current working directory.
#' @param skip_data Whether to skip returning columns
#'   representing the data inputs. If TRUE these are included and
#'   will have NA values for any missing years
#' @param use_names Whether to use the names written by the RE
#'   model or to use more meaningful ones (default)
#' @return A dataframe with columns named based on the tags in
#'   the output file
#' @example
#' x <- read_re_output()
#' str(x)
#' library(ggplot)
#' ggplot(x, aes(year, biomA))
read_re_output <- function(file='rwout.rep', skip_data=FALSE,
                           use_names=FALSE){
  if(!file.exists(file))
    stop(paste(file, "not found, check path and file name"))
  out <- readLines(file)
  colnames <- gsub(" ","", out[-grep(' ', out)])
  out <- out[grep(' ', out)]
  y <- lapply(out, function(x)  as.numeric(strsplit(x,' ')[[1]][-1]))
  ## Since the data may have missing values split data and
  ## estimates into two then recombine
  d <- data.frame(do.call(cbind, y[1:3]))
  names(d) <- colnames[1:3]
  o <- data.frame(do.call(cbind, y[4:11]))
  names(o) <- colnames[4:11]
  if(!skip_data)
    o <- merge(o,d, by.x='yrs', by.y='yrs_srv', all.x=TRUE)
  if(!use_names){
    x <- c('year', 'LCI','biomass', 'UCI', 'low90th', 'upp90th',
           'log_biomass', 'SE_log', 'srv_est', 'srv_sd')
    names(o) <- x[1:ncol(o)]
  }
  return(o)
}

#' Write an input file for the RE model
#'
#' @param file A path to the output file, usually ending in .dat
#' @param years A vector of years for the observations
#' @param biomass A vector of biomasses for the observations
#' @param CV A vector of CV values for the observations
#' @param tag Optional character string to print to file for
#'   bookkeeping (e.g., "EBS shelf flathead sole")
#' @param yr_start Optional start year for predictions
#' @param yr_end Optional end year for predictions
#'
#' @return Nothing, a data file is written as \param{file}
#'
#' @seealso read_re_output
write_re_input <- function(file, years, biomass, CV, tag=NULL,
                           yr_start=NULL, yr_end=NULL) {
  ## input checks
  stopifnot(all.equal(length(years), length(biomass), length(CV)))
  if(!all.equal(years, sort(years)))
    stop("years is not sorted")
  stopifnot(all(is.finite(CV)))
  stopifnot(all(is.finite(years)))
  stopifnot(all(is.finite(biomass)))
  stopifnot(all(biomass>=0))
  if(is.null(yr_start)) yr_start <- min(years)
  if(is.null(yr_end)) yr_end <- max(years)
  stopifnot(yr_start<= min(years))
  stopifnot(yr_end>=max(years))

  out <- file
  write(paste("#", survey), out,ncolumns =  1 )

  ## Start new file
  N <- length(years)
  write(noquote(paste("# RE input file written on", date(),
                      "by write_re_input function")), out, append=FALSE)
  write(noquote(paste("# For",sp, survey)) , out, append=TRUE)
  ## syr and nyr
  write(noquote("# Start year and end year"), out, ncolumns=N, append=TRUE)
  write(paste(yr_start, yr_end), out, ncolumns=N, append=TRUE)
  ##nobs
  write(noquote("# Nobs survey"), out, append=TRUE)
  write(length(years), out, ncolumns=N, append=TRUE)
  ##Years
  write(noquote("# Years"),out, append=TRUE)
  write(years, out, ncolumns=N, append=TRUE)
  ##biomass
  write(noquote("# Survey biomass"), out, append=TRUE)
  write(biomass, out, ncolumns=N, append=TRUE)
  ##CV
  write(noquote("# Survey CV"), out, append=TRUE)
  write(CV, out, ncolumns=N, append=TRUE)
  message(paste("Finished writing RE input for", survey, sp))
}
