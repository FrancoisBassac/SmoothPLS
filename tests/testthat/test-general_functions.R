#### Basis orthonormalization ####

test_that("create_bspline_basis creates a valid basis object", {
  start <- 0
  end <- 100
  nbasis <- 10
  norder <- 4

  basis <- create_bspline_basis(start, end, nbasis, norder)

  expect_s3_class(basis, "basisfd")
  expect_equal(basis$nbasis, nbasis)
  expect_equal(as.numeric(basis$rangeval), c(start, end))
})

test_that("gram_schmidt_orthonormalize produces an orthonormal basis", {
  # 1. Create a non-orthogonal B-spline basis
  basis <- create_bspline_basis(0, 10, nbasis = 5, norder = 4)

  # Ensure it's not already orthonormal (B-splines are not)
  expect_false(is_orthonormal(basis))

  # 2. Orthonormalize
  ortho_basis_list <- gram_schmidt_orthonormalize(basis,
                                                  output_type = "fdlist")

  # 3. Check properties
  expect_type(ortho_basis_list, "list")
  expect_length(ortho_basis_list, 5)
  expect_s3_class(ortho_basis_list[[1]], "fd")

  # Verification using your helper function
  expect_true(is_orthonormal(ortho_basis_list))
})

test_that("is_orthonormal handles edge cases and invalid inputs", {
  # Should stop if input is not valid
  expect_error(is_orthonormal(matrix(1:4, 2)), "needs 'basisfd' or 'fd' object list")
})

test_that("evaluate_metric returns correct dimensions and positive definiteness", {
  basis <- create_bspline_basis(0, 1, nbasis = 4)
  metric <- evaluate_metric(basis)

  # Metric matrix should be nbasis x nbasis
  expect_equal(dim(metric), c(4, 4))

  # Symmetry check
  expect_equal(metric, t(metric))

  # Positive eigenvalues (positive definite matrix)
  expect_true(all(eigen(metric)$values > 0))
})

#### assemble_basis_metric  ####

test_that("assemble_basis_metric builds correct block-diagonal matrix", {
  # 1. Setup two different bases
  basis1 <- create_bspline_basis(0, 100, nbasis = 4, norder = 4)
  basis2 <- create_bspline_basis(0, 100, nbasis = 5, norder = 4)
  basis_list <- list(basis1, basis2)

  # 2. Assemble metric
  combined_metric <- assemble_basis_metric(basis_list)

  # Check dimensions: 4 (basis1) + 5 (basis2) = 9
  expect_equal(dim(combined_metric), c(9, 9))

  # Check block-diagonal structure (off-diagonal blocks should be 0)
  expect_equal(sum(combined_metric[1:4, 5:9]), 0)
  expect_equal(sum(combined_metric[5:9, 1:4]), 0)

  # 3. Test with selection (curves_to_keep)
  # Keeping only the first basis
  selected_metric <- assemble_basis_metric(basis_list, curves_to_keep = list(1))
  expect_equal(dim(selected_metric), c(4, 4))
})

test_that("obj_list_creation correctly replicates objects", {
  basis <- create_bspline_basis(0, 10, nbasis = 5)
  n_rep <- 3

  res_list <- obj_list_creation(n_rep, basis)

  expect_length(res_list, n_rep)
  expect_s3_class(res_list[[1]], "basisfd")
  expect_identical(res_list[[1]], res_list[[2]])
})

test_that("from_fd_to_func creates a working R function", {
  basis <- create_bspline_basis(0, 10, nbasis = 4)
  coef <- c(1, 0, 0, 0) # First basis function

  # From explicit coef and basis
  func <- from_fd_to_func(coef = coef, basisobj = basis)
  expect_type(func, "closure")

  # Check evaluation at t=0
  val_at_0 <- func(0)
  expect_type(val_at_0, "double")

  # From fd object
  fd_obj <- fda::fd(coef, basis)
  func_from_fd <- from_fd_to_func(fd_obj = fd_obj)
  expect_equal(func(5), func_from_fd(5))
})
