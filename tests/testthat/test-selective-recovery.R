test_that("dropout mask preserves observed values",{
 x <- matrix(c(5,0,4,3,0,2),nrow=2)
 membership <- c(1,1,1)
 s <- local_alra_score(x,membership,rank=1)
 m <- select_dropout_mask(s,threshold=0.95)
 y <- recover_dropout_expression(x,m,membership)
 expect_equal(y[matrix(FALSE,nrow(x),ncol(x))],numeric(0))
 expect_equal(dim(y),dim(x))
})

test_that("dropout killer returns complete object",{
 x <- matrix(rpois(50,2),nrow=5)
 membership <- rep(1,10)
 z <- dropout_killer(x,membership,rank=1)
 expect_true(all(c("expression","score","mask") %in% names(z)))
 expect_equal(dim(z$expression),dim(x))
})
