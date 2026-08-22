test_that("weighted recovery preserves dimensions",{
 x <- matrix(rpois(200,1),20,10)
 emb <- matrix(rnorm(30),10,3)
 mem <- rep(1:2,each=5)
 y <- weighted_neighbor_prediction(x,emb,mem)
 expect_equal(dim(y),dim(x))
})

test_that("mask only modifies selected entries",{
 x <- matrix(1:9,3,3)
 m <- matrix(FALSE,3,3);m[1,1]<-TRUE
 y <- selective_recovery(x,m,matrix(2,3,3),matrix(3,3,3))
 expect_equal(y[2:3,],x[2:3,])
})
