# Data Discretization for Risk Score Models

## About this Repository

`discretize.R` contains the main function, `discretize()`. 

`example_usage.Rmd` decribes the functionality and usage of `discretize()`.

`test_adherence.Rmd` and `test_discretize.Rmd` are test files used during the development of `discretize()` to compare performance between `risk_mod()` models using more involved preprocessing and discretization, respectively. 

## Future Improvements

- [ ] Perform more error checking for arguments to `discretize()`. For example, check that `threshold` is numeric, or check that all strings in `continuous_cols` are indeed column names in `X`.
- [ ] Implement a method to determine the exact range of values in a bucket. For example, if a column called `age` was split into 3 buckets, we would be able to see the precise ranges that each bucket represented (say, `<18`, `18-24`, `25+`).