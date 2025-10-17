
<!-- README.md is generated from README.Rmd. Please edit that file -->

# readxl <a href="https://readxl.tidyverse.org"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/readxl)](https://cran.r-project.org/package=readxl)
[![R-CMD-check](https://github.com/tidyverse/readxl/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tidyverse/readxl/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/tidyverse/readxl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/tidyverse/readxl?branch=main)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

## Overview

The readxl package makes it easy to get data out of Excel and into R.
Compared to many of the existing packages (e.g. gdata, xlsx,
xlsReadWrite) readxl has no external dependencies, so it’s easy to
install and use on all operating systems. It is designed to work with
*tabular* data.

readxl supports both the legacy `.xls` format and the modern xml-based
`.xlsx` format. The [libxls](https://github.com/libxls/libxls) C library
is used to support `.xls`, which abstracts away many of the complexities
of the underlying binary format. To parse `.xlsx`, we use the
[RapidXML](https://rapidxml.sourceforge.net/) C++ library.

## Installation

The easiest way to install the latest released version from CRAN is to
install the whole tidyverse.

``` r
install.packages("tidyverse")
```

NOTE: you will still need to load readxl explicitly, because it is not a
core tidyverse package loaded via `library(tidyverse)`.

Alternatively, install just readxl from CRAN:

``` r
install.packages("readxl")
```

Or install the development version from GitHub:

``` r
#install.packages("pak")
pak::pak("tidyverse/readxl")
```

## Cheatsheet

You can see how to read data with readxl in the **data import
cheatsheet**, which also covers similar functionality in the related
packages readr and googlesheets4.

<a href="https://github.com/rstudio/cheatsheets/blob/main/data-import.pdf"><img src="https://raw.githubusercontent.com/rstudio/cheatsheets/main/pngs/thumbnails/data-import-cheatsheet-thumbs.png" width="630" height="252"/></a>

## Usage

``` r
library(readxl)
```

readxl includes several example files, which we use throughout the
documentation. Use the helper `readxl_example()` with no arguments to
list them or call it with an example filename to get the path.

``` r
readxl_example()
#>  [1] "clippy.xls"    "clippy.xlsx"   "datasets.xls"  "datasets.xlsx"
#>  [5] "deaths.xls"    "deaths.xlsx"   "geometry.xls"  "geometry.xlsx"
#>  [9] "type-me.xls"   "type-me.xlsx"
readxl_example("clippy.xls")
#> [1] "/private/tmp/Rtmpzymw54/temp_libpath117f779f0191c/readxl/extdata/clippy.xls"
```

`read_excel()` reads both xls and xlsx files and detects the format from
the extension.

``` r
xlsx_example <- readxl_example("datasets.xlsx")
read_excel(xlsx_example)
#> # A tibble: 32 × 11
#>     mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  21       6   160   110  3.9   2.62  16.5     0     1     4     4
#> 2  21       6   160   110  3.9   2.88  17.0     0     1     4     4
#> 3  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1
#> # ℹ 29 more rows

xls_example <- readxl_example("datasets.xls")
read_excel(xls_example)
#> # A tibble: 32 × 11
#>     mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  21       6   160   110  3.9   2.62  16.5     0     1     4     4
#> 2  21       6   160   110  3.9   2.88  17.0     0     1     4     4
#> 3  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1
#> # ℹ 29 more rows
```

List the sheet names with `excel_sheets()`.

``` r
excel_sheets(xlsx_example)
#> [1] "mtcars"   "chickwts" "quakes"
```

Specify a worksheet by name or number.

``` r
read_excel(xlsx_example, sheet = "chickwts")
#> # A tibble: 71 × 2
#>   weight feed     
#>    <dbl> <chr>    
#> 1    179 horsebean
#> 2    160 horsebean
#> 3    136 horsebean
#> # ℹ 68 more rows
read_excel(xls_example, sheet = 3)
#> # A tibble: 1,000 × 5
#>     lat  long depth   mag stations
#>   <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1 -20.4  182.   562   4.8       41
#> 2 -20.6  181.   650   4.2       15
#> 3 -26    184.    42   5.4       43
#> # ℹ 997 more rows
```

There are various ways to control which cells are read. You can even
specify the sheet here, if providing an Excel-style cell range.

``` r
read_excel(xlsx_example, n_max = 3)
#> # A tibble: 3 × 11
#>     mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  21       6   160   110  3.9   2.62  16.5     0     1     4     4
#> 2  21       6   160   110  3.9   2.88  17.0     0     1     4     4
#> 3  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1
read_excel(xlsx_example, range = "C1:E4")
#> # A tibble: 3 × 3
#>    disp    hp  drat
#>   <dbl> <dbl> <dbl>
#> 1   160   110  3.9 
#> 2   160   110  3.9 
#> 3   108    93  3.85
read_excel(xlsx_example, range = cell_rows(1:4))
#> # A tibble: 3 × 11
#>     mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  21       6   160   110  3.9   2.62  16.5     0     1     4     4
#> 2  21       6   160   110  3.9   2.88  17.0     0     1     4     4
#> 3  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1
read_excel(xlsx_example, range = cell_cols("B:D"))
#> # A tibble: 32 × 3
#>     cyl  disp    hp
#>   <dbl> <dbl> <dbl>
#> 1     6   160   110
#> 2     6   160   110
#> 3     4   108    93
#> # ℹ 29 more rows
read_excel(xlsx_example, range = "mtcars!B1:D5")
#> # A tibble: 4 × 3
#>     cyl  disp    hp
#>   <dbl> <dbl> <dbl>
#> 1     6   160   110
#> 2     6   160   110
#> 3     4   108    93
#> # ℹ 1 more row
```

If `NA`s are represented by something other than blank cells, set the
`na` argument.

``` r
read_excel(xlsx_example, na = "0")
#> # A tibble: 32 × 11
#>     mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  21       6   160   110  3.9   2.62  16.5    NA     1     4     4
#> 2  21       6   160   110  3.9   2.88  17.0    NA     1     4     4
#> 3  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1
#> # ℹ 29 more rows
```

If you are new to the tidyverse conventions for data import, you may
want to consult the [data import
chapter](https://r4ds.had.co.nz/data-import.html) in R for Data Science.
readxl will become increasingly consistent with other packages, such as
[readr](https://readr.tidyverse.org/).

## Articles

Broad topics are explained in [these
articles](https://readxl.tidyverse.org/articles/index.html):

- [Cell and Column
  Types](https://readxl.tidyverse.org/articles/cell-and-column-types.html)
- [Sheet
  Geometry](https://readxl.tidyverse.org/articles/sheet-geometry.html):
  how to specify which cells to read
- [readxl
  Workflows](https://readxl.tidyverse.org/articles/articles/readxl-workflows.html):
  Iterating over multiple tabs or worksheets, stashing a csv snapshot

We also have some focused articles that address specific aggravations
presented by the world’s spreadsheets:

- [Column
  Names](https://readxl.tidyverse.org/articles/articles/column-names.html)
- [Multiple Header
  Rows](https://readxl.tidyverse.org/articles/articles/multiple-header-rows.html)

## Features

- No external dependency on, e.g., Java or Perl.

- Re-encodes non-ASCII characters to UTF-8.

- Loads datetimes into POSIXct columns. Both Windows (1900) and
  Mac (1904) date specifications are processed correctly.

- Discovers the minimal data rectangle and returns that, by default.
  User can exert more control with `range`, `skip`, and `n_max`.

- Column names and types are determined from the data in the sheet, by
  default. User can also supply via `col_names` and `col_types` and
  control name repair via `.name_repair`.

- Returns a
  [tibble](https://tibble.tidyverse.org/reference/tibble.html), i.e. a
  data frame with an additional `tbl_df` class. Among other things, this
  provide nicer printing.

## Other relevant packages

Here are some other packages with functionality that is complementary to
readxl and that also avoid a Java dependency.

**Writing Excel files**: The example files `datasets.xlsx` and
`datasets.xls` were created with the help of
[openxlsx](https://CRAN.R-project.org/package=openxlsx) (and Excel).
openxlsx provides “a high level interface to writing, styling and
editing worksheets”.

``` r
l <- list(mtcars = mtcars, chickwts = chickwts, quakes = quakes)
openxlsx::write.xlsx(l, file = "inst/extdata/datasets.xlsx")
```

[writexl](https://cran.r-project.org/package=writexl) is a new option in
this space, first released on CRAN in August 2017. It’s a portable and
lightweight way to export a data frame to xlsx, based on
[libxlsxwriter](https://github.com/jmcnamara/libxlsxwriter). It is much
more minimalistic than openxlsx, but on simple examples, appears to be
about twice as fast and to write smaller files.

**Non-tabular data and formatting**:
[tidyxl](https://cran.r-project.org/package=tidyxl) is focused on
importing awkward and non-tabular data from Excel. It also “exposes cell
content, position and formatting in a tidy structure for further
manipulation”.
