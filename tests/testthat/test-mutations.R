test_that("parse_aa_mut", {
  mut <- "G1:T10I,G2:D20N,G3:Q30E,G1:A40T,G3:P50L,G2:G60R"
  expect_snapshot(parse_aa_mut(mut))
})
