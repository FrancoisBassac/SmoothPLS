#### evaluate_lambda ####

test_that("evaluate_lambda throws error if curve_type is missing", {
  df_test <- data.frame(id = 1, time = c(0, 10), state = c(1, 1))
  basis <- create_bspline_basis(0, 10, nbasis = 4)

  # Check input validation
  expect_error(evaluate_lambda(df = df_test, basis = basis,
                               curve_type = NULL, int_mode = 2))
  expect_error(evaluate_lambda(df = df_test, basis = basis,
                               curve_type = NULL, int_mode = 1))
})


test_that("evaluate_lambda correctly calculates projection for CFD", {
  # 1. Individual staying in state 1 from t=0 to t=100
  df_test <- data.frame(id = 1, time = c(0, 100), state = c(1, 1))
  # Constant basis (nbasis=1, value=1)
  basis_const <- fda::create.constant.basis(rangeval = c(0, 100))

  res <- evaluate_lambda(df = df_test, basis = basis_const, curve_type = 'cat')

  # Integral of 1 * 1 over [0, 100] = 100
  expect_equal(as.numeric(res[1,1]), 100, tolerance = 1e-5)

  # 2. Individual changing state: 0 -> 1 -> 0
  # Active in state 1 only between t=20 and t=50 (duration 30)
  df_change <- data.frame(id = 1, time = c(0, 20, 50, 100),
                          state = c(0, 1, 0, 0))
  res_change <- evaluate_lambda(df = df_change,
                                basis = basis_const, curve_type = 'cat')

  # Integral should be 1 * 30 = 30
  expect_equal(as.numeric(res_change[1,1]), 30, tolerance = 1e-5)
})

test_that("evaluate_lambda correctly calculates projection for SFD (Scalar)", {
  # Create a linear signal X(t) = t
  time_vec <- seq(0, 10, length.out = 11)
  df_num <- data.frame(id = 1, time = time_vec, value = time_vec)

  # Constant basis (phi(t) = 1)
  basis_const <- fda::create.constant.basis(rangeval = c(0, 10))

  # mode = 2 uses pracma::trapz as recommended for SFD
  res <- evaluate_lambda(df = df_num, basis = basis_const,
                         curve_type = 'num', int_mode = 2)

  # Integral of t * 1 from 0 to 10 = [t^2/2] = 100/2 = 50
  expect_equal(as.numeric(res[1,1]), 50, tolerance = 1e-5)
})

test_that("evaluate_lambda handles multiple individuals correctly", {
  # Two individuals with different active durations
  df_multi <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(0, 100, 0, 100),
    state = c(1, 1, 1, 1) # Ind 1: active 100, Ind 2: will be modified below
  )
  # Modify Ind 2 to be active only half the time
  df_multi <- data.frame(
    id = rep(c(1, 2), each = 2),
    time = rep(c(0, 100), 2),
    state = c(1, 1, 0, 0) # Ind 1: full active, Ind 2: full inactive
  )

  basis_const <- fda::create.constant.basis(rangeval = c(0, 100))
  res <- evaluate_lambda(df = df_multi, basis = basis_const, curve_type = 'cat')

  expect_equal(nrow(res), 2)
  expect_equal(as.numeric(res[1,1]), 100) # Ind 1
  expect_equal(as.numeric(res[2,1]), 0)   # Ind 2
})

# tests/testthat/test-lambda.R

test_that("evaluate_lambda sensitivity to integration parameters", {
  # Testing if changing subdivisions in mode 1 (integrate)
  # doesn't break the function
  df_test <- data.frame(id = 1, time = c(0, 100), state = c(1, 1))
  basis <- fda::create.bspline.basis(c(0, 100), nbasis = 4)
  expect_silent(
    evaluate_lambda(
      df = df_test,
      basis = basis,
      curve_type = 'cat',
      int_mode = 1,
      subdivisions = 200
    )
  )
})

test_that("evaluate_lambda stops with appropriate messages for invalid inputs",{
  df_test <- data.frame(id = 1, time = c(0, 10), state = c(1, 1))
  basis <- fda::create.constant.basis(rangeval = c(0, 10))
  # Test for missing curve_type (stop is caught by expect_error)
  # The message must match the stop() message in the R code
  expect_error(
    evaluate_lambda(df = df_test, basis = basis, curve_type = NULL),
    "curve_type should be 'cat' or 'num'"
  )
  expect_error(
    evaluate_lambda(df = df_test, basis = basis, curve_type = "invalid"),
    "curve_type should be 'cat' or 'num'"
  )
})
