myecdf.ksCI <- function(x)
{
  n <- length(x)
  ec <- ecdf(x)
  xx <- get("x", envir=environment(ec))# = sort(x)
  yy <- get("y", envir=environment(ec))
  D <- approx.ksD(n)
  yyu <- pmin(yy+D, 1)
  yyl <- pmax(yy-D, 0)
  ecu <- stepfun(xx, c(yyu, 1) )
  ecl <- stepfun(xx, c(yyl, yyl[n]) )
  return(list(ecu = ecu, ecl = ecl, ec = ec))
}


approx.ksD <- function(n)
{
  ## approximations for the critical level for Kolmogorov-Smirnov
  ## statistic D,
  ## for confidence level 0.95. Taken from Bickel & Doksum, table IX,
  ## p.483
  ## and Lienert G.A.(1975) who attributes to Miller,L.H.(1956), JASA
  ifelse(n > 80,
         1.358 /( sqrt(n) + .12 + .11/sqrt(n)),##Bickel&Doksum, table
         ##IX,p.483
         
         splinefun(c(1:9, 10, 15, 10 * 2:8),# from Lienert
                   c(.975,   .84189, .70760, .62394, .56328,# 1:5
                     .51926, .48342, .45427, .43001, .40925,# 6:10
                     .33760, .29408, .24170, .21012,# 15,20,30,40
                     .18841, .17231, .15975, .14960)) (n))
}

