\name{NEWS}
\title{News for Package 'ACDm'}

\section{Changes in ACDm version 1.0.4.2 (2022-11-14)}{
  \itemize{
    \item Minor fix to make comments in C code valid UTF-8.
  }
}

\section{Changes in ACDm version 1.0.4 (2016-07-15)}{

    \itemize{
      \item One can now use exogenous regressors when fitting ACD models using \code{acdFit()}. See the function's help page.
      \item Fixed a bug in the internal C functions getLL_ACDcall, getLL_LACD1call, and getLL_LACD2call that unabled estimating ACD(p, q), LACD1(p, q), and LACD2(p, q) models when p differed from q.
      \item Added the \code{predict()} method.
      \item Added the \code{plotLL()} function.
      \item Changed how the \code{trace} argument inside the \code{control} list in the \code{acdFit()} function works. It will now plot the MLE search path, recorded as the parameter values each time the optimization function calls the log likelihood.
    } 

}

\section{Changes in ACDm version 1.0.3 (2015-12-01)}{

    \itemize{
      \item The examples were set to "dontrun" for a few of the functions.
      \item Changed the format of the time columns in the example data \code{durData} and \code{adjDurData} from POSIXlt to character strings. The time zone feature of the POSIXlt format made some of the examples to fail in other time zones - preventing the package to be installed from source within those time zones.}
    
	  \subsection{For the function \code{computeDurations()}:}{
    	  \itemize{
          \item The function  now works properly for sub-second precision.
          \item Fixed bug where the duration of the second transaction of a day would sometimes be calculated as the duration from the opening time
          \item Fixed bug where the last duration would sometimes not be added.
          \item Having only one of the price/volume column should not work.
          \item Now checks if any negative durations were computed - if so, checks that the data were given in chronological order.
    }
    }
    \subsection{For the function \code{diurnalAdj()}:}{
    	  \itemize{
          \item for \code{method = "cubicSpline"}: now checks if the time stamps of the durations are within the range of the \code{nodes} argument. 
    }
    }
  
}
