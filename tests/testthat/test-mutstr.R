test_that("mutstr", {
  raw_mut_string <- c(
    variant1 = "T10I,D20N,Q30E,A40T,P50L,G60R",
    variant2 = "T10I,D20-,Q30E,A40T,P50L,G60R,S80R",
    variant3 = "T10A,D20G,Q30E,A40T,P50L,G60R"
  )
  m <- mutstr(raw_mut_string, sep = ",")

  expect_identical(names(m), names(raw_mut_string))
  expect_identical(mstr(m), raw_mut_string)
  expect_identical(
    mut(m),
    list(
      variant1 = c("T10I", "D20N", "Q30E", "A40T", "P50L", "G60R"),
      variant2 = c("T10I", "D20-", "Q30E", "A40T", "P50L", "G60R", "S80R"),
      variant3 = c("T10A", "D20G", "Q30E", "A40T", "P50L", "G60R")
    )
  )
  expect_true(is(m[1:2], "mutstr"))
  expect_identical(m[[1]], c("T10I", "D20N", "Q30E", "A40T", "P50L", "G60R"))
  expect_snapshot(m)
})
