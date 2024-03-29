% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LASSO2_reg.R
\name{LASSO2_reg}
\alias{LASSO2_reg}
\title{LASSO2 Variable Selection and Regular Regression Modeling}
\usage{
LASSO2_reg(
  data = NULL,
  standardization = FALSE,
  columnWise = TRUE,
  biomks = NULL,
  outcomeType = c("binary", "continuous", "time-to-event"),
  Y = NULL,
  time = NULL,
  event = NULL,
  nfolds = 10,
  outfile = "nameWithPath"
)
}
\arguments{
\item{data}{A data matrix or a data frame where samples are in rows, and features/traits are in columns.}

\item{standardization}{A logical variable to indicate if standardization is needed before variable selection. The default is FALSE.}

\item{columnWise}{A logical variable to indicate if column-wise or row-wise normalization is needed. The default is TRUE, which performs column-wise normalization. 
This is only meaningful when "standardization" is TRUE.}

\item{biomks}{A vector of potential biomarkers for variable selection. They should be a subset of column names in "data".}

\item{outcomeType}{The outcome variable type. There are three choices: "binary" (default), "continuous", and "time-to-event".}

\item{Y}{The outcome variable name when the outcome type is either "binary" or "continuous".}

\item{time}{The time variable name when the outcome type is "time-to-event".}

\item{event}{The event variable name when the outcome type is "time-to-event".}

\item{nfolds}{The number of folds for cross-validation. The default value is 10.}

\item{outfile}{A string for the output file including the path, if necessary, but without the file type extension.}
}
\value{
A list is returned with the same output as from confirmVars.
\item{fit}{A model with selected variables for the given outcome variable.}
\item{allplot}{A list with all plots.}
}
\description{
This function performs variable selection with LASSO2 but ignores the shrunk LASSO coefficients and builds a regular regression model.
}
\details{
The first part of LASSO2_reg involves variable selection with LASSO2, typically based on the mean lambda.1se from 10 iterations of 
n-fold cross-validation-based LASSO regression. In each iteration, a lambda.1se refers to the largest value of lambda such that the error 
is within 1 standard error of the minimum. However, if only one or no variable is selected, the cross-validation results are ignored, 
and a minimum of two remaining variables is ensured through full-data lambda simulations.
The second part of LASSO2_reg involves ignoring the shrunk LASSO coefficients and building a regular regression model.
It is suitable for three types of outcomes: continuous, binary, and time-to-event.
}
\examples{
# Load in data sets:
data("datlist", package = "csmpv")
tdat = datlist$training

# The function saves files locally. You can define your own temporary directory. 
# If not, tempdir() can be used to get the system's temporary directory.
temp_dir = tempdir()
# As an example, let's define Xvars, which will be used later:
Xvars = c("highIPI", "B.Symptoms", "MYC.IHC", "BCL2.IHC", "CD10.IHC", "BCL6.IHC")

# The function can work with three different outcome types. 
# Here, we use time-to-event as an example:
# tlr = LASSO2_reg(data = tdat, biomks = Xvars,
#                  outcomeType = "time-to-event",
#                  time = "FFP..Years.",event = "Code.FFP",
#                  outfile = paste0(temp_dir, "/survivalLASSO2_reg"))
# You might save the files to the directory you want.
# To delete the "temp_dir", use the following:
unlink(temp_dir)
}
\references{
Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.

Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.

Hastie, T. J. and Pregibon, D. (1992) Generalized linear models. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. Springer-Verlag, 2000.

Kassambara A, Kosinski M, Biecek P (2021). survminer: Drawing Survival Curves using 'ggplot2'_. R package version 0.4.9,
}
\author{
Aixiang Jiang
}
