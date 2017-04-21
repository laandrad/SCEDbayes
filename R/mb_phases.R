mb.create.phases <- function(P, s) {
    ## create dummy ABAB variables and Time according to Moeyaert et al. (2014)
    nPoints = length(P)
    x = P[-1] - P[-nPoints]
    nPhase = which(!x==0)
    P1 = P
    T1 = s
    T2 = T1 - (nPhase[1]+1)

    list(P1, T1, T2)
}
