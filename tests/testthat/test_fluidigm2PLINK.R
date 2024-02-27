
test_that("fluidigm2PLINK processes files correctly", {
  # Define the paths to the input and expected output files
    file_path_csv <- system.file("extdata", "example_data.csv", package = "Fluidigm")
    file_path_map <- system.file("extdata", "PlateD_withY.map", package = "Fluidigm")

    outdir <- tempdir()

    fluidigm2PLINK(file=file_path, map=file_path_map, outdir=outdir)

# The expected outputs was calculated using this function during the preparation of the test
  outfiles <- c("example_data.csv.additional_summary_stats.png",
                "example_data.csv.call_rate_markers.png",
                "example_data.csv.genotyping_success_samples.png",
                "example_data.csv.map",
                "example_data.csv.ped")

  # Compute the MD5 checksum of each file
  # expected_output_md5 <- unlist(lapply(file.path(outdir, outfiles), function(file) digest::digest(file, algo = "md5", file = TRUE)))
  expected_output_md5 <- c("2d876f1b4bf83f44d06f740586227231",
                           "7c222738fa14f4f5d014f9b60e6cc689",
                           "78fdf6a054f89dbcd2b8d52feb785902",
                           "7843662622f1bb7ba1855ca4ea6aca84",
                           "58f50b73062f245b52c8669c3018a538")

  # Compute the MD5 checksums of the output files
  output_md5 <- lapply(file.path(outdir, outfiles), function(file) digest::digest(file, algo = "md5", file = TRUE))

  # Check that the MD5 checksums match
  expect_equal(output_md5, expected_output_md5)

  # Clean up the temporary files
  unlink(file.path(outdir, outfiles))
})

