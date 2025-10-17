###
### $Id: strcmp.R 22 2022-05-30 18:03:47Z proebuck $
###


##-----------------------------------------------------------------------------
test.strcmp <- function(input, expected) {
    output <- do.call(getFromNamespace("strcmp", "matlab"), input)
    identical(output, expected)
}

test.strcmp(list(S = "foo", T = "foo"), TRUE)
test.strcmp(list(S = "foo", T = "bar"), FALSE)
test.strcmp(list(S = c("foo", "bar"),
                 T = c("foo", "bar")), TRUE)
# Case matters...
test.strcmp(list(S = c("foo", "bar"),
                 T = c("FOO", "BAR")), FALSE)
# Number of elements of each must match...
test.strcmp(list(S = c("foo", "bar"),
                 T = c("foo", "bar", "baz")), FALSE)
test.strcmp(list(S = c("foo", "bar", "baz"),
                 T = c("xxx", "bar", "xxx")), FALSE)

