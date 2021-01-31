### relative entropy function
### https://link.springer.com/article/10.1007%2FBF01246098 - relative entropy

### assumes an object of type poLCA LCA result

relative_entropy <- function(lcares) {

  lcares2 <- lcares

  entropy <- function(p) sum(-p * log(p))

  error_prior <- entropy(lcares2[["P"]]) # Class proportions

  error_post <- mean(apply(lcares2[["posterior"]], 1, entropy))

  R2_entropy <- (error_prior - error_post) / error_prior

  R2_entropy
  
}