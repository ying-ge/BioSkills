###
### $Id: tictoc.R 22 2022-05-30 18:03:47Z proebuck $
###


##-----------------------------------------------------------------------------
test.toc <- function(delay, expected) {
    Sys.sleep(delay)
    output <- matlab::toc(FALSE)
    all.equal(output,
              expected,
              tolerance = 0.10)
}

matlab::tic()
try(test.toc(4, 4))
try(test.toc(2, 4+2))

