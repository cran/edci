# delete duplicate rows from a matrix
deldupMclust = function(clust,  prec = NULL,  ncol = NULL,  dz = TRUE) {

  c = class(clust)

  if (is.matrix(clust) && dz && any(colnames(clust) == "value")) {
    if (sum(clust[1:nrow(clust),  "value"] != 0) == 0)
      clust = clust[1,  ]
    else
      clust = clust[1:nrow(clust) * (clust[1:nrow(clust),  "value"] != 0),  ]
  }

  if (is.matrix(clust)) {

    if (all(is.null(dimnames(clust))) |
        all(dimnames(clust)[[2]] != "count")) {
      count = rep(1,  nrow(clust))
      clust = cbind(clust,  count)
    }

    if (is.null(ncol)) {
      if (!is.null(c) && c == "oregMclust")
        ncol = 2
      else if (!is.null(c) && c == "circMclust")
        ncol = 3
      else
        ncol = ncol(clust) - sum(colnames(clust) %in% c("value",  "count",  "proj",  "env"))
    }
    
    if (ncol > ncol(clust))
      ncol = ncol(clust)
    
    if (!is.null(prec)) {
      clust[,  1:ncol] = round(clust[,  1:ncol],  prec)
      if (!is.null(c) && c == "oregMclust")
        clust[, 1] = clust[, 1] %% round(2 * pi, prec)
    }

    result = clust
    count = 1
    
    for (i in 2:nrow(clust)) {
      j = 1
      ins = TRUE
      while (j <= count && ins) {
        if (is.na(all(clust[i, 1:ncol] == result[j, 1:ncol])) |
            all(clust[i, 1:ncol] == result[j, 1:ncol])) {
          ins = FALSE
          if (!is.na(all(clust[i, 1:ncol] == result[j, 1:ncol])))
            result[j, "count"] = result[j, "count"] + clust[i, "count"]
        }
        j = j + 1
      }
      if (ins) {
        count = count + 1
        result[count, ] = clust[i, ]
      }
    }
    result = result[1:count, ]
    if (count > 1)
      class(result) = c
  }
  else
    result = clust
  
  result
}

# calculation of regression clusters based on Orthogonal regression
oregMclust = function(datax,  datay,  bw,  method = "const", xrange = range(datax),
  yrange = range(datay),  prec = 4, na = 1,  sa = NULL,  nl = 10,  nc = NULL, brmaxit = 1000) {

  cat("Break with <CTRL>-C (linux) or <ESC> (windows)\n")

  if (is.matrix(datax)) {
    if (ncol(datax) >= 3)
      sa = datax[, 3]
    
    datay = datax[, 2]
    datax = datax[, 1]
  }
 
  n = min(length(datax), length(datay))
  
  # factors for normalizing the data
  cf = (xrange[2] - xrange[1] + yrange[2] - yrange[1])/2/pi 

  count = integer(1)

  if (method != "prob") {
    m = 0

    # normalize starting angles
    if (!is.null(sa)) {
      sa = as.double(sa + pi/2)
      m = 1
    }

    # return values for C-procedure
    if (method == "all") {
      alpha = double(n * (n + 1)/2)
      beta = double(n * (n + 1)/2)
      value = double(n * (n + 1)/2)
      m = 3
    }
    else if (m == 1) {
      alpha = double(n)
      beta = double(n)
      value = double(n)
    } else {
      alpha = double(n * na)
      beta = double(n * na)
      value = double(n * na)
    }

    #if ((m == 1) | ((m == 0) && (na == 1))) {
    #  xp = double(n)
    #  yp = double(n)
    #}
    #else {
    #  xp = NULL
    #  yp = NULL
    #}

    # calculate regression lines
    result = .C("c_oregMclust", 
      as.double(datax/cf), 
      as.double(datay/cf), 
      as.integer(n), 
      as.double(bw/cf), 
      as.integer(m), 
      as.integer(na), 
      sa, 
      0.0, 
      0.0, 
      as.integer(brmaxit), 
      alpha = alpha, 
      beta = beta, 
      value = value, 
      count = count, 
      #xp, 
      #yp, 
      PACKAGE = "edci")

    # rescale return values
    beta = result$beta * cf
    #xp = xp * cf
    #yp = yp * cf

    alpha = result$alpha
    value = result$value
    count = result$count

    # normalize parameteres
    for (i in 1:count) {
      if (beta[i] < 0) {
        beta[i] = beta[i] * (-1)
        alpha[i] = alpha[i] + pi
      }
    }
    alpha = alpha %% (2 * pi)
  } else { # method == "prob"
    alpha = double(n * (n + 1)/2)
    beta = double(n * (n + 1)/2)
    value = double(n * (n + 1)/2)

    m = 2
    
    noregs = -1
    count = 0

    anew = double(1)
    bnew = double(1)
    vnew = double(1)
    repeat {
      for (i in 1:nl) {
        repeat {
          i1 = round(runif(1, min = 1, max = n))
          i2 = round(runif(1, min = 1, max = n))
          if (i1 != i2)
            break
        }
        if (datay[i1] == datay[i2])
          a = pi/2
        else
          a = atan(-(datax[i2] - datax[i1])/(datay[i2] - datay[i1]))
        b = (cos(a) * datax[i1] + sin(a) * datay[i1])/cf

        count = count + 1

        result = .C("c_oregMclust",
          as.double(datax/cf),
          as.double(datay/cf),
          n,
          as.double(bw/cf),
          as.integer(m),
          0L,
          0.0,
          a,
          b,
          as.integer(brmaxit),
          anew = anew,
          bnew = bnew,
          vnew = vnew,
          count = count,
          #xp,
          #yp,
          PACKAGE = "edci")

        # rescale
        bnew = result$bnew * cf

        anew = result$anew
        vnew = result$vnew
        count = result$count

        # normalize
        if (bnew < 0) {
          bnew = -bnew
          anew = anew + pi
        }
        anew = anew %% (2 * pi)
        if (anew > pi) {
          anew = anew - pi
          bnew = -bnew
        }
        alpha[count] = anew
        beta[count] = bnew
        value[count] = vnew
      }

      # count different
      noregsnew = nrow(deldupMclust(cbind(round(alpha, prec),
        round(beta, prec))[1:count, ], prec=prec))
      if (is.null(noregsnew))
        noregsnew = 1
      cat("Found clusters: ", noregsnew, "\n")

      if (is.null(nc) && noregs == noregsnew)
        break
      if (!is.null(nc))
        if (noregsnew >= nc)
          break

      noregs = noregsnew
    }
  }

  value = -value
  reg = cbind(alpha, beta, value)[1:count, ]
  class(reg) = "oregMclust"

  # delete equivalant regression lines
  reg = deldupMclust(reg, prec = prec, ncol = 2)

  reg
}

plot.oregMclust = function(x, datax, datay, prec = 3, rcol = "black", rlty = 1, rlwd = 3, ...) {

  myline = function (alpha, beta) {
    if (alpha != 0)
      abline(beta/sin(alpha), -cos(alpha)/sin(alpha), col = rcol, lwd = rlwd, lty = rlty)
    else
      abline(v = beta, col = rcol, lwd = rlwd, lty = rlty)
  }

  if (is.matrix(datax)) {
    datay = datax[, 2]
    datax = datax[, 1]
  }

  plot(datax, datay, ...)
  if (is.matrix(x)) {
    for (i in 1:nrow(x))
      myline(x[i, 1], x[i, 2])
  } else {
    myline(x[1], x[2])
  }
}

print.oregMclust = function(x, ...) {
  print(x[], ...)
}

regparm = function(reg) {
  result = reg
  if (is.matrix(reg)) {
    result[, 1] = reg[, 2]/sin(reg[, 1])
    result[, 2] = -cos(reg[, 1])/sin(reg[, 1])
    colnames(result)[1] = "intersept"
    colnames(result)[2] = "slope"  
  } else {
    result[1:2] = c(reg[2]/sin(reg[1]), -cos(reg[1])/sin(reg[1]))
    names(result)[1] = "intersept"
    names(result)[2] = "slope"  
  }
  result
}

# calculation of circle clusters
circMclust = function(datax, datay, bw, method = "const", prec = 4,
  minsx = min(datax), maxsx = max(datax), nx = 10, minsy = min(datay), maxsy = max(datay),
  ny = 10, minsr = 0.01 * max(datax, datay), maxsr = (max(datax, datay) - min(datax, datay)),
  nr = 10, nsc = 5, nc = NULL, minsd = NULL, maxsd = NULL, brminx = minsx, brmaxx = maxsx, 
  brminy = minsy, brmaxy = maxsy, brminr = minsr, brmaxr = maxsr, brmaxit = 1000) {

  cat ("Break with <CTRL>-C (linux) or <ESC> (windows)\n")

  n = min(length(datax), length(datay))

  count = integer(1)

  if (method=="all") {
    m = 0

    # return values for C-procedure
    nmax = n * (n - 1)/2
    nmax = nmax*(nmax - 1)/2
    cx = double(nmax)
    cy = double(nmax)
    r = double(nmax)
    value = double(nmax)

    result = .C("c_oregMcirc", 
      as.double(datax), 
      as.double(datay), 
      n, 
      as.double(bw), 
      as.integer(m), 
      0L, 
      0L, 
      0L, 
      as.double(maxsr), 
      0.0, 
      0.0, 
      0.0, 
      0L, 
      as.double(brminx), as.double(brmaxx),
      as.double(brminy), as.double(brmaxy),
      as.double(brminr), as.double(brmaxr),
      as.integer(brmaxit),
      cx = cx,
      cy = cy,
      r = r,
      value = value,
      count = count,
      PACKAGE = "edci")

    cx = result$cx
    cy = result$cy
    r = result$r
    value = result$value
    count = result$count
  } else if (method == "prob") {
    m = 1

    nmax = 100
    cx = double(nmax)
    cy = double(nmax)
    r = double(nmax)
    value = double(nmax)

    nocirc = -1
    count = 0

    cxnew = double(1)
    cynew = double(1)
    rnew = double(1)
    vnew = double(1)
    c = integer(1)

    if (!is.null(minsd))
      minsdq = minsd^2
    if (!is.null(maxsd)) 
      maxsdq = maxsd^2

    repeat {
      i = 1
      while (i <= nsc) {
        repeat {
          i1 = round(runif(1, min = 1, max = n))
          i2 = round(runif(1, min = 1, max = n))
          i3 = round(runif(1, min = 1, max = n))
          if (i1 != i2 && i1 != i3 && i2 != i3)
            break
        }

        if ((is.null(minsd) ||
          (((datax[i1] - datax[i2])^2 + (datay[i1] - datay[i2])^2) > minsdq &&
            ((datax[i1] - datax[i3])^2 + (datay[i1] - datay[i3])^2) > minsdq &&
            ((datax[i3] - datax[i2])^2 + (datay[i3] - datay[i2])^2) > minsdq)) &&
            (is.null(maxsd) ||
          (((datax[i1] - datax[i2])^2 + (datay[i1] - datay[i2])^2) < maxsdq &&
            ((datax[i1] - datax[i3])^2 + (datay[i1] - datay[i3])^2) < maxsdq &&
            ((datax[i3] - datax[i2])^2 + (datay[i3] - datay[i2])^2) < maxsdq))) {

          result = .C("c_oregMcirc", 
            as.double(datax),
            as.double(datay),
            n,
            as.double(bw),
            as.integer(m),
            as.integer(i1 - 1),
            as.integer(i2 - 1),
            as.integer(i3 - 1),
            as.double(maxsr),
            0.0,
            0.0,
            0.0,
            0L,
            as.double(brminx), as.double(brmaxx),
            as.double(brminy), as.double(brmaxy),
            as.double(brminr), as.double(brmaxr),
            as.integer(brmaxit),
            cxnew = cxnew,
            cynew = cynew,
            rnew = rnew,
            vnew = vnew,
            c = c,
            PACKAGE = "edci")

          cxnew = result$cxnew
          cynew = result$cynew
          rnew = result$rnew
          vnew = result$vnew
          c = result$c

          if (c > 0) {
            count = count + 1

            if (length(cx) < count) {
              cx = c(cx, double(nmax))
              cy = c(cy, double(nmax))
              r = c(r, double(nmax))
              value = c(value, double(nmax))
            }

            cx[count] = cxnew
            cy[count] = cynew
            r[count] = rnew
            value[count] = vnew
          }
          i = i + 1
        }
      }

      # count different
      nocircnew = nrow(deldupMclust(cbind(cx, cy, r)[1:count, ], prec = prec))

      if (is.null(nocircnew))
        nocircnew = 1
      cat("Found clusters: ", nocircnew, "\n")
      
      if (is.null(nc) && nocirc == nocircnew)
        break
      if (!is.null(nc))
        if (nocircnew >= nc)
          break

      nocirc = nocircnew
    }
  } else if (method == "const") {
    m = 2

    # return values for C-procedure
    #nmax = ceiling((maxsx-minsx)/xstep)*
    #  ceiling((maxsy-minsy)/ystep)*
    #    ceiling((maxsr-minsr)/rstep)
    nmax = 100
    cx = double(nmax)
    cy = double(nmax)
    r = double(nmax)
    value = double(nmax)
    count = integer(1)
    circ = NULL

    xseq = seq(minsx, maxsx, length = nx)
    yseq = seq(minsy, maxsy, length = ny)
    rseq = seq(minsr, maxsr, length = nr)
    startpoints = matrix(c(rep(xseq, each = ny * nr), rep(rep(yseq, each = nr), nx), 
      rep(rseq, nx * ny)), ncol = 3)

    for (i in 1:(ceiling(nrow(startpoints)/nmax))) {
      s = seq(((i - 1) * nmax + 1), (min(i * nmax, nrow(startpoints))))

      cat("calculating clusters for startingvalues ",
        s[1], "-", s[length(s)], " of ", nrow(startpoints), "...\n")

      result = .C("c_oregMcirc",
        as.double(datax),
        as.double(datay),
        n,
        as.double(bw),
        as.integer(m),
        0L,
        0L,
        0L,
        as.double(maxsr),
        as.double(startpoints[s, 1]),
        as.double(startpoints[s, 2]),
        as.double(startpoints[s, 3]),
        as.integer(length(s)),
        as.double(brminx), as.double(brmaxx),
        as.double(brminy), as.double(brmaxy),
        as.double(brminr), as.double(brmaxr),
        as.integer(brmaxit),
        cx = cx,
        cy = cy,
        r = r,
        value = value,
        count = count,
        PACKAGE = "edci")

      cx = result$cx
      cy = result$cy
      r = result$r
      value = result$value
      count = result$count

      if (count > 0) {
        if (count == 1)
          circnew = c(cx[1], cy[1], r[1], value[1], count)
        else
          circnew = deldupMclust(cbind(cx,  cy,  r,  value)[1:count, ], prec = prec, ncol = 3)
        if (i == 1)
          circ = circnew
        else 
          circ = deldupMclust(rbind(circ, circnew), prec = prec, ncol = 3)
      }
    }
    if (is.null(circ)) {
      cat ("no cluster found!\n")
    } else {
      cat ("finished\n")
      if (is.matrix(circ))
        circ[, 4] = -circ[, 4]
      else
        circ[4] = -circ[4]
    }
  } else {
    cat ("unknown method\n")
  }

  if (method != "const") {
    value = -value
    circ = cbind(cx, cy, r, value)[1:count, ]

    # delete equivalant regression lines  
    circ = deldupMclust(circ, prec = prec, ncol = 3)
  }

  if (!is.null(circ)) {
    class(circ) = "circMclust"
    rownames(circ) = NULL
  }
  circ
}

plot.circMclust = function(x, datax, datay, ccol = "black", clty = 1, clwd = 3, ...) {
  circle = function(cx,  cy,  r,  ...) {
    z = (0:360 * pi)/180
    x = sin(z) * r + cx
    y = cos(z) * r + cy
    lines(x,  y,  ...)
    invisible()
  }

  plot(datax, datay, xlim = range(c(datax, datay)), ylim = range(c(datax, datay)), ...)

  if (is.matrix(x)) {
    for (i in 1:nrow(x))
      circle(x[i, 1], x[i, 2], x[i, 3], col = ccol, lwd = clwd, lty = clty)
  } else {
    circle(x[1], x[2], x[3], col = ccol, lwd = clwd, lty = clty)
  }
}

print.circMclust = function(x, ...) {
  print(x[], ...)
}

# choose 'best' clusters
bestMclust = function(clust, nc = 1, crit = "value") {
  c = class(clust)

  if (!is.matrix(clust)) {
    result = clust
  } else {
    if (!(crit %in% colnames(clust))) {
      cat ("No row \"", crit, "\" found!\n")
      result = clust
    } else {
      result = rbind(clust[order(clust[, crit], decreasing = TRUE)[1:min(nc, nrow(clust))], ])
      class(result) = c
    }
  }
  result
}

# count projected points
projMclust = function(clust, x, y) {
  if (is.matrix(clust)) {
    cl = class(clust)

    if (all(is.null(dimnames(clust))) | all(dimnames(clust)[[2]]!="proj")) {
      proj = rep(0, nrow(clust))
      clust = cbind(clust, proj)
    }
    clust[, "proj"] = 0

    if (cl == "oregMclust") {
      a = clust[, 1]
      b = clust[, 2]
      for (i in 1:min(length(x), length(y))) {
        xp = x[i] - (cos(a) * x[i] + sin(a) * y[i] - b) * cos(a)
        yp = y[i] - (cos(a) * x[i] + sin(a) * y[i] - b) * sin(a)
        d = sqrt((x[i] - xp)^2 + (y[i] - yp)^2)
        c = order(d)[1]
        clust[c, "proj"] = clust[c, "proj"] + 1
      }
    } else {
      a1 = clust[, 1]
      a2 = clust[, 2]
      b = clust[, 3]
      for (i in 1:min(length(x), length(y))) {
        d = abs(sqrt((x[i] - a1)^2 + (y[i] - a2)^2) - b)
        c = order(d)[1]
        clust[c, "proj"] = clust[c, "proj"] + 1
      }

      projrel = clust[, "proj"]/clust[, "r"]
      if (all(is.null(dimnames(clust))) | all(dimnames(clust)[[2]] != "projrel")) {
        clust = cbind(clust, projrel)
      } else {
        clust[, "projrel"] = projrel
      }
    }
    class(clust) = cl
  }
  clust
}

# count points in environment
envMclust = function(clust, x, y, dist = 0) {
  if (is.matrix(clust)) {
    cl = class(clust)
    if (all(is.null(dimnames(clust))) | all(dimnames(clust)[[2]] != "env")) {
      env = rep(0, nrow(clust))
      clust=cbind(clust, env)
    }

    clust[, "env"] = 0

    if (cl == "oregMclust") {
      for (j in 1:nrow(clust)) {
        a = clust[j, 1]
        b = clust[j, 2]
        xp = x - (cos(a) * x + sin(a) * y - b) * cos(a)
        yp = y - (cos(a) * x + sin(a) * y - b) * sin(a)
        d = sqrt((x - xp)^2 + (y - yp)^2)
        clust[j, "env"] = clust[j, "env"] + sum(d <= dist)
      }
    } else {
      for (j in 1:nrow(clust)) {
        a1 = clust[j, 1]
        a2 = clust[j, 2]
        b = clust[j, 3]
        d = abs(sqrt((x - a1)^2 + (y - a2)^2) - b)
        clust[j, "env"] = clust[j, "env"] + sum(d <= dist)
      }
      envrel = clust[, "env"]/clust[, "r"]
      if (all(is.null(dimnames(clust))) | all(dimnames(clust)[[2]] != "envrel")) {
        clust = cbind(clust, envrel)
      } else {
        clust[, "envrel"] = envrel
      }
    }

    class(clust) = cl
  }
  clust
}
