test_that("smoothPLS_predict handles time range mismatch", {
  # Training on [0, 100]
  df_train <- generate_X_df(nind = 5, start = 0, end = 100, curve_type = 'cat')
  basis <- create_bspline_basis(0, 100, nbasis = 5)
  Y <- rnorm(5)
  res <- smoothPLS(df_list = df_train, Y = Y, basis_obj = basis,
                   curve_type_obj = 'cat', print_steps = FALSE,
                   plot_rmsep = FALSE, print_nbComp = FALSE,
                   plot_reg_curves = FALSE)

  # Predict on an individual with data up to 110 (outside basis range)
  df_out <- data.frame(id = 1, time = c(0, 110), state = c(1, 1))

  # This should ideally warn the user or handle the truncation
  expect_error(smoothPLS_predict(df_predict_list = df_out,
                                 delta_list = res$reg_obj,
                                 curve_type_obj = 'cat'))
})
