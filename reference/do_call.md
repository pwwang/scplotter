# Call a function with a list of arguments

Different from `do.call`, this is a faster version, especially when
there are big objects in the list.

## Usage

``` r
do_call(fn, args, quote = FALSE, envir = parent.frame())
```

## Arguments

- fn:

  A function to call

- args:

  A list of arguments to pass to the function

- quote:

  Whether to quote the arguments

- envir:

  The environment to evaluate the function in

## Value

The result of the function call
