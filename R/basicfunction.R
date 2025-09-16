# MCP and soft thresholding functions -----------------------------

mcp <- function(x, lam, a = 3) {
  b <- abs(x)
  z <- soft(x, lam) / (1 - 1 / a)
  z[b > (a * lam)] <- x[b > (a * lam)]
  z
}

soft <- function(a, b) {
  c <- abs(a) - b
  c[c < 0] <- 0
  c * sign(a)
}
