\name{joinreads}
\alias{joinreads}
\title{
Function to run the application FLASH to overlap paired-end reads
}
\description{
Takes a pair of sequence files and overlaps them to produce one contiguous sequence. Overlapping can be in any orientation if the user specifies (i.e. user can determine if read pairs can be combined opposite the sequencing direction or in the same direction of sequencing.)
}
\usage{
joinreads(File1, File2, out.dir, threads, output.prefix, allow.outies=TRUE, mismatch.ratio = 0.25, minimum.outie.overlap = 35)
}
\arguments{
  \item{File1}{character string giving the name of the first of the two paired-end files.}
  \item{File2}{character string giving the name of the second of the two paired-end files.}
  \item{out.dir}{path giving the name of the destination directory.}
  \item{threads}{integer: sets the number of working threads.}
  \item{output.prefix}{character string indicating the desired prefix for ouput files.}
  \item{allow.outies}{logical: indicates if overlapping end orientation will extend with or opposite sequencing orientation. If TRUE matches will be attempted with reads extending in either orientation, with the higher-scoring set chosen as the overlap.}
  \item{mismatch.ratio}{proportion: the maximum allowed ratio between the number of mismatched base pairs and the overlap length. Any occurrence of an 'N' in either read is ignored.}
  \item{minimim.outie.overlap}{integer: the minimum required overlap length between two reads in an 'outie' orientation to provide a confident overlap. Recommended minimum is 35.}
}
\details{
This function 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
