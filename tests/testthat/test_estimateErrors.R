
test_that("estimateErrors processes files correctly", {
  # Define the paths to the input and expected output files
  file_path_csv <- system.file("extdata", "example_data.csv", package = "Fluidigm")
  file_path_map <- system.file("extdata", "example_data_withY.map", package = "Fluidigm")

  outdir <- tempdir()
  fluidigm2PLINK(file=file_path_csv, map=file_path_map, outdir=outdir)

  # Define the parameter sets for the estimateErrors function
  params <- list(
    list(file=file.path(outdir, "example_data.csv.ped"), keep.rep=2,
         expected_md5=c("0b54bc9392e625e87af8294bd129b574",
                        "030be3d9bea8de6433cf4e60b9e9678e",
                        "5f938e1707b13f1ac790b7ec5362714b",
                        "f55c71551e3e8826e055cbb074834a1c")),
    list(file=file.path(outdir, "example_data.csv.ped"), keep.rep=2, sexing=TRUE,
         y.marker=c("Y_scaffoldY158711_762", "Y_scaffoldY42647_3017", "Y_scaffoldY42656_3986"),
         x.marker=c("X_scaffold11905_7659", "X_scaffold17088_4621", "X_scaffold1915_14108",
                    "X_scaffold4825_648", "X_scaffold5374_1437", "X_scaffold10171:3154"),
         expected_md5=c("7b027a2fdbc55201af39449ad90d971c",
                        "bbd30f6954941c7cb627be23a4953930",
                        "5f938e1707b13f1ac790b7ec5362714b",
                        "021db5f2a7ce05f832a1dd71e6899fec"))
  )

  # Apply the estimateErrors function with each parameter set
  for (param_set in params) {

    # Now apply the estimateErrors function
    do.call(estimateErrors, param_set[-length(param_set)])

    # The expected outputs was calculated using this function during the preparation of the test
    outfiles <- c("example_data.csv.GOOD.map", "example_data.csv.GOOD.ped",
                  "example_data.csv.ped_samples_to_RERUN.txt", "example_data.csv.ped_summary_individuals.csv")

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
