
test_that("fluidigm2PLINK processes files correctly", {
  # Define the paths to the input and expected output files
    file_path_csv <- system.file("extdata", "example_data.csv", package = "Fluidigm")
    file_path_map <- system.file("extdata", "PlateD_withY.map", package = "Fluidigm")

    outdir <- tempdir()

    fluidigm2PLINK(file=file_path_csv, map=file_path_map, outdir=outdir)

# The expected outputs was calculated using this function during the preparation of the test
  outfiles <- c("example_data.csv.map",
                "example_data.csv.ped")

  # Compute the MD5 checksum of each file
  # Compute the MD5 checksum of each file's contents
  # expected_output_md5 <- unlist(lapply(file.path(outdir, outfiles), function(file) {
  #   digest::digest(file, algo = "md5", file = TRUE)
  # }))


  expected_output_md5 <- c("7843662622f1bb7ba1855ca4ea6aca84",
                           "58f50b73062f245b52c8669c3018a538")

  # Compute the MD5 checksums of the output files
  output_md5 <- unlist(lapply(file.path(outdir, outfiles), function(file) {
                          digest::digest(file, algo = "md5", file = TRUE)
                          }))

  # Check that the MD5 checksums match
  expect_equal(output_md5, expected_output_md5)

  # Clean up the temporary files
  unlink(file.path(outdir, outfiles))
})

