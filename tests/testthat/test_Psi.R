
n <- 5
X <- r_unif_sph(n = n, p = 2, M = 2)
Theta <- X_to_Theta(X)
dim(Theta) <- c(n, 1, 2)
set.seed(98720222)
A <- matrix(rnorm(4 * 7), nrow = 4, ncol = 7)
B1 <- tcrossprod(X[, , 1])
B2 <- tcrossprod(X[, , 2])
b <- acos(c(B1[upper.tri(B1)], B2[upper.tri(B2)]))

test_that("upper_tri_ind", {

  expect_equal(drop(upper_tri_ind(5)) + 1,
               {I <- matrix(1:5^2, 5, 5); I[upper.tri(I)]})
  expect_equal(drop(upper_tri_ind(7)) + 1,
               {I <- matrix(1:7^2, 7, 7); I[upper.tri(I)]})

})

test_that("sort_each_col and sort_index_each_col", {

  expect_equal(sphunif:::sort_each_col(A), apply(A, 2, sort))
  expect_equal(sphunif:::sort_index_each_col(A), apply(A, 2, order))

})

test_that("Psi with X/Theta and use_ind_tri = TRUE/FALSE", {

  expect_equal(Psi_mat(X, use_ind_tri = FALSE),
               Psi_mat(Theta, use_ind_tri = FALSE))
  expect_equal(Psi_mat(X, use_ind_tri = FALSE),
               Psi_mat(X, ind_tri = upper_tri_ind(n), use_ind_tri = TRUE))
  expect_equal(Psi_mat(Theta, use_ind_tri = FALSE),
               Psi_mat(Theta, ind_tri = upper_tri_ind(n), use_ind_tri = TRUE))
  expect_equal(Psi_mat(X, use_ind_tri = FALSE),
               Psi_mat(Theta, ind_tri = upper_tri_ind(n), use_ind_tri = TRUE))
  expect_equal(Psi_mat(Theta, use_ind_tri = FALSE),
               Psi_mat(X, ind_tri = upper_tri_ind(n), use_ind_tri = TRUE))

})

test_that("Psi with angle_diff = TRUE/FALSE", {

  expect_equal(Psi_mat(X, angles_diff = TRUE),
               Psi_mat(X, angles_diff = FALSE))
  expect_equal(Psi_mat(X, angles_diff = TRUE),
               Psi_mat(Theta, angles_diff = FALSE))
  expect_equal(acos(cos(Psi_mat(Theta, angles_diff = TRUE))),
               Psi_mat(Theta, angles_diff = FALSE))
  expect_equal(acos(cos(Psi_mat(Theta, angles_diff = TRUE))),
               Psi_mat(X))
  expect_equal(acos(cos(Psi_mat(Theta, angles_diff = TRUE,
                                ind_tri =  upper_tri_ind(n),
                                use_ind_tri = TRUE))),
               Psi_mat(X))

})

test_that("Psi with scalar_prod = TRUE/FALSE", {

  expect_equal(Psi_mat(Theta, scalar_prod = TRUE),
               Psi_mat(Theta, scalar_prod = FALSE))
  expect_equal(Psi_mat(X, scalar_prod = FALSE),
               Psi_mat(Theta, scalar_prod = TRUE))
  expect_equal(acos(Psi_mat(X, scalar_prod = TRUE)),
               Psi_mat(X, scalar_prod = FALSE))
  expect_equal(acos(Psi_mat(X, scalar_prod = TRUE)),
               Psi_mat(Theta))
  expect_equal(acos(Psi_mat(X, scalar_prod = TRUE,
                            ind_tri =  upper_tri_ind(n),
                            use_ind_tri = TRUE)),
               Psi_mat(X))

})

test_that("Psi reconstruction", {

  expect_equal(c(Psi_mat(Theta)), b)
  expect_equal(c(Psi_mat(X)), b)
  expect_equal(Psi_mat(X)[, 1],
               {D <- acos(1 - as.matrix(dist(X[, , 1]))^2 / 2);
               D[upper.tri(D)]})
  expect_equal(Psi_mat(X)[, 2],
               {D <- acos(1 - as.matrix(dist(X[, , 2]))^2 / 2);
               D[upper.tri(D)]})
  expect_equal(Psi_mat(Theta)[, 1],
               {D <- acos(1 - as.matrix(dist(X[, , 1]))^2 / 2);
               D[upper.tri(D)]})
  expect_equal(Psi_mat(Theta)[, 2],
               {D <- acos(1 - as.matrix(dist(X[, , 2]))^2 / 2);
               D[upper.tri(D)]})
})

