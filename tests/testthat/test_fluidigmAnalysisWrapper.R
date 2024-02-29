test_that("fluidigmAnalysisWrapper processes files correctly", {
  # Define the paths to the input and expected output files
    file_path_csv <- system.file("extdata", "example_data.csv", package = "Fluidigm")
    file_path_map <- system.file("extdata", "example_data_withY.map", package = "Fluidigm")

    outdir <- tempdir()

    fluidigmAnalysisWrapper(file=file_path_csv, map=file_path_map, outdir=outdir,
                            keep.rep=2,
                            sexing=TRUE,
                            y.marker = c("Y_scaffoldY158711_762",
                                         "Y_scaffoldY42647_3017",
                                         "Y_scaffoldY42656_3986"),
                            x.marker = c("X_scaffold11905_7659",
                                         "X_scaffold17088_4621",
                                         "X_scaffold1915_14108",
                                         "X_scaffold4825_648",
                                         "X_scaffold5374_1437",
                                         "X_scaffold10171:3154"))

# The expected outputs was calculated using this function during the preparation of the test
  outfiles <- c("example_data.csv.map",
                "example_data.csv.ped",
                "example_data.csv.GOOD.map",
                "example_data.csv.GOOD.mibs",
                "example_data.csv.GOOD.pairs",
                "example_data.csv.GOOD.ped",
                "example_data.csv.GOOD.similar_0.85.genotypes.rout")

  expected_output_md5 <- c("0b54bc9392e625e87af8294bd129b574",
                           "b15387bc87408614b57659b7452bdd5e",
                           "7b027a2fdbc55201af39449ad90d971c",
                           "686629d9ae2b7b00ed12adeda793da3f",
                           "7ca9702338c74fd5bc4aa28aba4c0ff8",
                           "bbd30f6954941c7cb627be23a4953930",
                           "80a90c0f5e65393ff0a2b2ec60fb18b2")

  # Compute the MD5 checksums of the output files
  output_md5 <- unlist(lapply(file.path(outdir, outfiles), function(file) {
                          digest::digest(file, algo = "md5", file = TRUE)
                          }))

  # Check that the MD5 checksums match
  expect_equal(output_md5, expected_output_md5)

  # Clean up the temporary files
  unlink(file.path(outdir, outfiles))
})

