edgepoints = function(data, h1n, h2n, asteps = 4, estimator  = "kernel", kernel = "mean",
  score = "gauss", sigma = 1, kernelfunc = NULL, margin = FALSE) {

  epDelta = function(x) {
    if (x < 0)
      -1
    else
      1
  }

  epAt = function(x, y) {
    if (x == 0) {
      if (y >= 0)
        pi/2
      else
        pi/2
    } else {
      atan(y/x)
    }
  }

  epR1 = function(theta, x, y) {
    sqrt(x^2 + y^2) * epDelta(x) * cos(epAt(x, y) - theta)
  }

  epR2 = function(theta, x, y) {
    sqrt(x^2 + y^2) * epDelta(x) * sin(epAt(x, y) - theta)
  }

  angle = matrix(double(length(data)), nrow=nrow(data))
  value = matrix(double(length(data)), nrow=nrow(data))

  es = NULL
  sc = NULL
  ms = sigma * max(data)

  if (estimator == "kernel")
    es = 0
  else if (estimator == "median")
    es = 1
  else if (estimator == "M_mean")
    es = 2
  else if (estimator == "M_median")
    es = 3
  else if (estimator == "test_mean")
    es = 5
  else if (estimator == "test_median")
    es = 6
  else
    stop("estimator \"", estimator, "\" unknown.")

  if (es==2 || es ==3) {
    if (score == "gauss") {
      sc = 0
    }
    if (score == "huber") {
      sc = 1
      #ms = sigma/2
    }
    if (score == "mean") {
      sc = 9
    }
  }

  env = ceiling(sqrt((h1n * nrow(data))^2 + (h2n * ncol(data))^2))

  kernmat = NULL
  if (kernel == "mean" || es >= 5) {
    kern = 0
  } else if (kernel == "linear") {
    kern = 1
  } else if (kernel == "linear2") {
    kern = 2
  } else if (kernel == "gauss") {
    kern = 3
  } else if (kernel == "func") {
    kern = 4
    kernmat = double(asteps * (2 * env + 1)^2)
    for (i in ((-env):env)) {           
      for(j in ((-env):env)) {        
        for (k in (0:(asteps - 1))) {   
          theta = -pi/2 + (k * pi/asteps)
          x = epR1(theta, i/nrow(data), j/ncol(data))/h1n
          y = epR2(theta, i/nrow(data), j/ncol(data))/h2n
          kernmat[k * (2 * env + 1)^2 + (i + env) * (2 * env + 1) + (j + env) + 1] =
            kernelfunc(2 * x, y)
        }
      }
    }
  } else {
    stop("kernel \"",kernel,"\" unknown.")
  }

  if (es == 1)
    kern = 0

  if (!is.null(es)) {
    result = .C("c_edgepoints",
      as.double(data),
      nrow(data),          
      ncol(data),
      as.integer(kern),     # kernel
      as.double(h1n),
      as.double(h2n),
      as.integer(es),
      as.integer(sc),       # Typ der Scorefunktion
      as.double(sigma),     # Sigma
      as.double(kernmat),   # Gewichtsmatrix
      as.double(ms),        # Max_Schritt
      as.integer(asteps),
      angle = angle,
      value = value,
      PACKAGE = "edci")
  }

  value = result$value
  angle = result$angle

  if (es == 5 || es == 6)
    value = -value

  if (margin == FALSE) {
    if (es == 5 || es == 6)
      v = 1
    else
      v = 0
    value[c(1:env,(nrow(value) - env + 1):nrow(value)), ] = v
    value[, c(1:env, (ncol(value) - env+1):ncol(value))] = v
  } else if (margin == "cut") {
    value = value[(env + 1):(nrow(value) - env), (env + 1):(ncol(value) - env)]
    angle = angle[(env + 1):(nrow(angle) - env), (env + 1):(ncol(angle) - env)]
  }
  list(value = value, angle = angle)
}

eplist = function(data, maxval, test = FALSE, xc = NULL, yc = NULL) {
  if (test == TRUE) {
    data[[1]] = -data[[1]]
    maxval = -maxval
  }

  n = sum(data[[1]] > maxval)
  if (is.null(xc))
    xc = seq(1/nrow(data[[1]]), 1, 1/nrow(data[[1]]))
  if (is.null(yc))
    yc = seq(1/ncol(data[[1]]), 1, 1/ncol(data[[1]]))

  o = order(data[[1]], decreasing = TRUE)[1:n]
  result = cbind(xc[(o - 1) %% nrow(data[[1]]) + 1], yc[(o - 1) %/% nrow(data[[1]]) +1 ],
    data[[2]][o])
  colnames(result) = c("x", "y", "angle")

  result
}
