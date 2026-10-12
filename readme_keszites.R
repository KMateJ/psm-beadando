# A README.md-t a psm_szimulacio.Rmd-bol generaljuk (ugyanaz a szoveg, kod nelkul).
# A Monte Carlo eredmenyt az eredmenyek/ mappabol olvassa, igy nem fut ujra.
rmarkdown::render("psm_szimulacio.Rmd", output_format = "github_document",
                  output_file = "README.md", encoding = "UTF-8")
