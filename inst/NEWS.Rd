\name{NEWS}
\title{News for Package 'ACDm'}

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
          \item for \code{method = "cubicSpline"}: now checks if the timestamps of the durations are within the range of the \code{nodes} argument. 
    }
    }
  
}
