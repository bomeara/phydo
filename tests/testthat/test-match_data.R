test_that("check_continuous works", {
  funky <- data.frame(
    runif(3),
    round(runif(3)),
    as.character(runif(3)),
    as.character(round(runif(3))),
    c("a", "b", "3"),
    c(NA, 1.4, 2.5),
    c(NA, "a", "b"),
    c(NA,4,5),
    stringsAsFactors=FALSE)
  funky[,1] <- as.numeric(funky[,1])
  funky[,2] <- as.numeric(funky[,2])
  funky[,6] <- as.numeric(funky[,6])
  colnames(funky) <- NULL
  expect_true(all(check_continuous(funky)[c(1,3,6)]))

})

test_that("phydo_fitGeiger works", {
  data(geospiza, package="geiger")
  geospiza$dat[4,3] <- NA
  phydo <- match_data(geospiza$phy, geospiza$dat)
  fc_result <- phydo_fitGeiger(phydo, keep="continuous")
  expect_true(inherits(fc_result[[1]][[1]], "gfit"))
})

test_that("gbif_taxon_query works", {
  result <- gbif_taxon_query("Puma", gbif_limit=600)
  expect_true(inherits(result$data$decimalLatitude), "numeric")
})
