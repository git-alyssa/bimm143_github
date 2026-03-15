# class06
Alyssa Duran (PID: A18550696)

## Background

All functions in R have at least 3 things:

- A **name** that we use to call the function
- One or more input **arguments**
- The **body** the lines of R code that do the work

## Our First Function

Let’s write a function called add() to add some numbers (the input
argument).

``` r
add <- function(x, y) {
x + y
}
```

Now we can use this function:

``` r
add(100, 1)
```

    [1] 101

``` r
add(c(100, 1, 100), 1)
```

    [1] 101   2 101

``` r
add(10, 10)
```

    [1] 20

> Q. What if I give a multiple element vector to x and y?

``` r
add(x = c(100, 1), y = c(100, 1))
```

    [1] 200   2

> Q. What if I give three inputs to the function?

``` r
#add(x = c(100, 1), y = 1, z = 1)
```

> Q. What if I give only one input to the function?

``` r
addnew <- function(x, y = 1) {
  x + y
  }
addnew(x = 100)
```

    [1] 101

``` r
addnew(c(100, 0), 100)
```

    [1] 200 100

If we write our function with input arguments having no default value
then the user will be required to set them when they use the function.
We can give our input arguments “default” values by setting them equal
to some sensible value - e.g. y=1 in the `addnew()` function.

## A Second Function

Let’s try something more interesting: Make a sequence generating tool…
The `sample()` function can be a useful starting point here:

``` r
sample(1:10, size = 4)
```

    [1] 6 1 8 7

> Q. Generate 9 random numbers taken from the input vector x=1:10?

``` r
sample(1:10, size = 9)
```

    [1]  8  5  3  9 10  4  6  2  7

> Q. Generate 12 random numbers taken from the input vector x=1:10?

``` r
sample(1:10, size = 12, replace = TRUE, prob = NULL)
```

     [1] 4 5 3 2 8 6 7 6 9 6 9 4

> Q. Write code for the sample() function that generates nucleotide
> sequences of length 6.

``` r
sample(c("A","T","C","G"), size = 6, replace = TRUE, prob = NULL)
```

    [1] "C" "G" "T" "C" "G" "G"

> Q. Write a first function generate_dna() that returns a user specified
> length DNA sequence.

``` r
generate_dna <- function(len) {
  sample(c("A","T","C","G"), size = len, replace = TRUE)
  }
generate_dna(len = 6)
```

    [1] "T" "A" "C" "A" "T" "T"

> **Key-Points** Every function in R lookds fundamentally the same in
> terms of its structure. Basically 3 things: name, input, and body

name \<- function(input) { body }

> Functions can have multiple inputs. These can be required arguments or
> optional arguments. With optional arguments having a set default
> value.

> Q. Modify and improve our generate_dna() function to return it’s
> generated sequence in a more standard format like “AGTAGTA” rather
> thant he vector “A”, “G”, “T”, “A”

``` r
generate_dna <- function(len = 6, fasta = TRUE) {
ans <- sample(c("A","T","C","G"),
              size = len,
              replace = TRUE)
if(fasta) {
  cat("Single-element vector output")
ans <- paste(ans, collapse = "")
} else {
  cat("Multi-element vector output")
  }
return(ans)
}
generate_dna(fasta = TRUE)
```

    Single-element vector output

    [1] "GGGGCT"

The `paste()` function - it’s job is to join up or stick together
(a.k.a. paste) input strings together.

``` r
paste("alice", "loves R", sep = " ")
```

    [1] "alice loves R"

Flow control means where the R brain goes in your code.

``` r
good_mood <- TRUE
if(good_mood) {
  cat("Great!")
  } else {
    cat("Bummer!")
    }
```

    Great!

Great!

## A Protein Generating Function

> Q. Write a function, called generate_protein(), that generates a user
> specified length protein sequence.

``` r
generate_protein <- function(len) {
# The amino-acids to sample from
aa <- c("A","R","N","D","C","Q",
        "E","G","H","I","L","K",
        "M","F","P","S","T","W","Y","V")
# Draw n=len amino acids to make our sequence
ans <- sample(aa,
              size = len,
              replace = TRUE)
ans <- paste(ans, collapse = "")
return(ans)
}
generate_protein(6)
```

    [1] "RTLFDK"

> Q. Use that function to genreate random protein sequences between
> length 6 and 12.

``` r
# Labor intensive, done manually
generate_protein(6)
```

    [1] "WECWPH"

``` r
generate_protein(7)
```

    [1] "LCMMTYR"

``` r
generate_protein(8)
```

    [1] "VERPSYIT"

``` r
generate_protein(9)
```

    [1] "TPWIIELTY"

``` r
generate_protein(10)
```

    [1] "WPTCKIAFDQ"

``` r
generate_protein(11)
```

    [1] "PTVLVTHPIQR"

``` r
generate_protein(12)
```

    [1] "WHHEQNLCQTIV"

``` r
# Using loops
for (i in 6:12) {
  # FASTA ID line ">id"
  cat(">", i, sep = ""
      ,
      "\n")
  # Protein sequence line
  cat(generate_protein(i), "\n")
  }
```

    >6
    RGHDEF 
    >7
    HMRDTCG 
    >8
    RAKLCCYA 
    >9
    GITLNQPLG 
    >10
    DLAQKHMRTP 
    >11
    YKEAHIGYRRH 
    >12
    TYGPWFYHAVLC 

> Q. Are any of your sequences unique i.e. not found anywhere else in
> nature?

Sequences 6 (QSIDKK ) and 7 (FPNHLID) were found in nature with 100% Ide
and 100% coverage. Sequences 8 (EMIWTQGG), 9 (CQFGFHHQW), 10
(CCHQEFPEHY), 11 (NKLFVF- PWMRW), and 12 (IRSTRPIPCMEF) were unique and
not found anywhere else in nature.
