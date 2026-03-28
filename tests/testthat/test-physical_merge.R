test_that("empty input returns empty data frame", {
  df  <- data.frame(position = numeric(0), value = numeric(0))
  res <- physical_merge(df, sig_th = 0.05, window = 500)
  expect_equal(nrow(res), 0L)
  expect_named(res, c("serial", "start", "end", "rps_BP", "rps_value"))
})

test_that("single significant SNP produces one block", {
  df  <- data.frame(position = c(1000, 2000, 3000),
                    value    = c(0.5,  0.01, 0.5))
  res <- physical_merge(df, sig_th = 0.05, window = 500)
  expect_equal(nrow(res), 1L)
  expect_equal(res$rps_BP,    2000)
  expect_equal(res$rps_value, 0.01)
})

test_that("two distant peaks produce two blocks", {
  df  <- data.frame(position = c(1000, 100000),
                    value    = c(0.01, 0.01))
  res <- physical_merge(df, sig_th = 0.05, window = 500)
  expect_equal(nrow(res), 2L)
})

test_that("adjacent blocks do not overlap", {
  set.seed(1)
  df <- data.frame(
    position = sort(sample(1:1e6, 1000)),
    value    = runif(1000)
  )
  df$value[c(100, 110, 600, 610)] <- c(1e-9, 1e-8, 1e-9, 1e-8)
  res <- physical_merge(df, sig_th = 5e-8, window = 50000)
  if (nrow(res) > 1L) {
    overlaps <- res$end[-nrow(res)] > res$start[-1L]
    expect_false(any(overlaps), info = "Blocks must not overlap")
  }
})

test_that("reward = 'max' picks larger value as representative", {
  df  <- data.frame(position = c(100, 200, 300),
                    value    = c(5,   10,  7))
  res <- physical_merge(df, sig_th = 4, window = 500, reward = "max")
  expect_equal(res$rps_BP,    200)
  expect_equal(res$rps_value, 10)
})

test_that("invalid inputs throw errors", {
  expect_error(physical_merge(list(), sig_th = 0.05, window = 500))
  df <- data.frame(position = 1:3, value = 1:3)
  expect_error(physical_merge(df, sig_th = 0.05, window = 500,
                              reward = "wrong"))
  expect_error(physical_merge(df, sig_th = 0.05, window = -1))
})
