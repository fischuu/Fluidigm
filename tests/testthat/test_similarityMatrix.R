
test_that("similarityMatrix processes files correctly", {
  # Define the paths to the input and expected output files
  file_path_csv <- system.file("extdata", "example_data.csv", package = "Fluidigm")
  file_path_map <- system.file("extdata", "example_data_withY.map", package = "Fluidigm")

  outdir <- tempdir()
  fluidigm2PLINK(file=file_path_csv, map=file_path_map, outdir=outdir)

  # Define the parameter sets for the estimateErrors function
  params.ee <- list(file=file.path(outdir, "example_data.csv.ped"), keep.rep=2, sexing=TRUE,
                 y.marker=c("Y_scaffoldY158711_762", "Y_scaffoldY42647_3017", "Y_scaffoldY42656_3986"),
                 x.marker=c("X_scaffold11905_7659", "X_scaffold17088_4621", "X_scaffold1915_14108",
                            "X_scaffold4825_648", "X_scaffold5374_1437", "X_scaffold10171:3154"))
  do.call(estimateErrors, params.ee)

  # Now apply the calculatePairwiseSimilarities function
  params <- list(file=file.path(outdir, "example_data.csv.GOOD"), sexing=TRUE)
  do.call(calculatePairwiseSimilarities, params)

# Now apply the getPairwiseSimilarityLoci function
  params <- list(file=file.path(outdir, "example_data.csv.GOOD"), verbose=FALSE)
  do.call(getPairwiseSimilarityLoci, params)

  params <- list(
    list(file=file.path(outdir, "example_data.csv.GOOD"),
         expected_md5=c("80a90c0f5e65393ff0a2b2ec60fb18b2"))
  )

  # Apply the estimateErrors function with each parameter set
  for (param_set in params) {

    # Now apply the similarityMatrix function
    do.call(similarityMatrix, param_set[-length(param_set)])

    # The expected outputs was calculated using this function during the preparation of the test
    outfiles <- c("example_data.csv.GOOD.similar_0.85.genotypes.rout")

    # Compute the MD5 checksums of the output files
    output_md5 <- unlist(lapply(file.path(outdir, outfiles), function(file) {
      digest::digest(file, algo = "md5", file = TRUE)
    }))

    # Check that the MD5 checksums match
    expect_equal(output_md5, param_set$expected_md5)

    # Clean up the temporary files
    unlink(file.path(outdir, outfiles))
  }
})
