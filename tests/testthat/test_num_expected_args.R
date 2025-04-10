context(".num_expected_args")

test_that(
  desc = 'Number of arguments that a `fun` expects a user to supply matches
  expectations. That is, the number of non-... arguments that a function expects
  that don\'t have a default value.',
  code = {

    # Function has two arguments and no defaults.
    weighted_mean_no_defaults <- function(df, weighted) {
      if(!weighted) {
        out <- mean(df$value)
      } else {
        out <- weighted.mean(x = df$value, w = df$value)
      }
      out
    }

    # Function has two arguments, but a non-NULL default is supplied
    # for one of them.
    weighted_mean_non_null_defaults <- function(df, weighted = TRUE) {
      if(!weighted) {
        out <- mean(df$value)
      } else {
        out <- weighted.mean(x = df$value, w = df$value)
      }
      out
    }

    # Function has two arguments, but a NULL default is supplied for
    # one of them.
    weighted_mean_null_defaults <- function(df, weighted = NULL) {
      if(is.null(weighted)) {
        out <- mean(df$value)
      } else {
        out <- weighted.mean(x = df$value, w = df$value)
      }
      out
    }

    # Function has two arguments and no defaults. Should count both
    # arguments without defaults.
    expect_equal(.num_expected_args(weighted_mean_no_defaults), 2)

    # Function has two arguments, but a non-NULL default is supplied
    # for one of them. Should only count the one argument without a
    # default.
    expect_equal(.num_expected_args(weighted_mean_non_null_defaults), 1)

    # Function has two arguments, but a NULL default is supplied for
    # one of them. Should only count the one argument without a default.
    expect_equal(.num_expected_args(weighted_mean_null_defaults), 1)
  }
)


