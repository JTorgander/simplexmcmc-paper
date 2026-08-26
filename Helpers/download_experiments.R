# Download and extract fitted RDHMC models from Zenodo

record_id <- "22108727"              # Zenodo record ID
zip_filename <- "Experiments.zip"  # Name of the uploaded ZIP file
output_dir <- "."

dir.create(output_dir, showWarnings = FALSE)

zip_path <- file.path(output_dir, zip_filename)
url <- sprintf(
  "https://zenodo.org/records/%s/files/%s?download=1",
  record_id,
  zip_filename
)

if (!file.exists(zip_path)) {
  message("Downloading fitted models")
  download.file(url, zip_path, mode = "wb")
}

unzip(zip_path, exdir = output_dir)

message("Files downloaded and extracted")
