edgecluster <- function(data,
                        h1n, h2n,
                        maxval,
                        bw         = max(h1n,h2n)/qnorm(0.975),
                        asteps     = 4,
                        estimator  = "M_median",
                        kernel     = "gauss",
                        score      = "gauss",
                        sigma      = 1,
                        kernelfunc = NULL){

  if (estimator=="test_mean" || estimator=="test_median")
    test=TRUE
  else
    test=FALSE
  
  ep <- eplist(edgepoints(data,h1n,h2n,
                          asteps     = asteps,
                          estimator  = estimator,
                          kernel     = kernel,
                          score      = score,
                          sigma      = sigma,
                          kernelfunc = kernelfunc),
               maxval,test=test)

  list(oregMclust(ep,bw=bw),ep)
}
