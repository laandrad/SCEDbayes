create.phases <- function(P, s) {
    ## create dummy ABAB variables and Time according to Moeyaert et al. (2014)
    nPoints = length(P)
    x = P[-1] - P[-nPoints]
    nPhase = which(!x==0)
    P1 = c(rep(0, nPhase[1]), rep(1, nPoints - nPhase[1]))
    P2 = c(rep(0, nPhase[2]), rep(1, nPoints - nPhase[2]))
    P3 = c(rep(0, nPhase[3]), rep(1, nPoints - nPhase[3]))
    T1 = s
    T2 = T1 - (nPhase[1]+1)
    T3 = T1 - (nPhase[2]+1)
    T4 = T1 - (nPhase[3]+1)

    list(P1, P2, P3, T1, T2, T3, T4)
}
