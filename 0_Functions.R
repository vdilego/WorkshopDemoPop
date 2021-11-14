# --------------------------------------------------------------------------------
# Useful functions used in this project
# --------------------------------------------------------------------------------

### Life table

life.table <- function(mx){
  ax <- c(0.14, rep(0.5, length(mx)-1))
  qx <- mx/(1+(1-ax)*mx)
  qx[length(qx)] <- 1
  qx[qx > 1] <- 1
  px <- 1-qx
  lx <- c(100000, (cumprod(px)*100000)[1:(length(px)-1)])
  dx <- c(-diff(lx), lx[length(lx)])
  Lx1 <- lx[-1]+ax[-length(ax)]*dx[-length(dx)]
  open.Lx <-  ifelse( mx[length(mx)] == 0, 0, dx[length(dx)]/mx[length(mx)])
  Lx <- c(Lx1, open.Lx)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx/lx

  return(data.frame(qx=qx, px = px, ax = ax, lx = lx , dx = dx, Lx= Lx,
                    Tx = Tx, ex = ex))
}


# 5-years age group

age.cat <- function(x, lower = 0, upper, by = 10,
                    sep = "-", above.char = "+") {

  labs <- c(paste(seq(lower, upper - by, by = by),
                  seq(lower + by - 1, upper - 1, by = by),
                  sep = sep),
            paste(upper, above.char, sep = ""))

  cut(floor(x), breaks = c(seq(lower, upper, by = by), Inf),
      right = FALSE, labels = labs)
}





