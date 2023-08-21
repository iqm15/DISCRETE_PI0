p <- function(u) {
  if (u < 1/2) {res = 0}
  else if (u >= 1/2 & u <1) {res = 1/2}
  else if (u >= 1) {res = 1}
  return(res)
  
}

p_new <- function(u) {
  if (u < 1/2) {res = 0}
  else if (u >= 1/2 & u <4/3) {res = 3/8}
  else if (u >= 4/3) {res = 1}
  return(res)
  
}