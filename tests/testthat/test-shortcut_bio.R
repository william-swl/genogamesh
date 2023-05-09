test_that("nt2aa", {
  expect_identical(
    nt2aa(c("ATGAAA", "TTGCCC", "CTGTTT")), c("MK", "LP", "LF")
  )
})
