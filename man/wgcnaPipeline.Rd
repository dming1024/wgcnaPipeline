\name{wgcnaPipeline}
\alias{wgcnaPipeline}
\title{Pipeline of WGCNA}
\usage{
wgcnaPipeline(GeneMatrix ,probeAnnotation ,traits,method=mean,sftEstimate=20)
}
\description{
Using Gene expression data in GEO database with clinical traits to
    figure out the co-expression network associated with the initiation and development of diseases.

}

\arguments{
\item{GeneMatrix}{
GEO dataset must be Provided.Columns of gene expression matrix are Samples, and rows are ProbeID. And so for the parament only support GEO dataset
}
\item{probeAnnotation}{
Including 2 clomuns.The first one is ProbeID, which can be matched with rownames of GeneMatrix;and the second one is definitely offical gene Symbol.
}
\item{traits}{
The features of Samples in GeneMatrix.And also includes 1-N columns.The first one is the names of Samples which must be consistent with the columns of GeneMatrix;and the rest of columns include the Features or Traits.
}
\item{method}{
The duplicated genes names frequently occures in the GeneMatrix. Therefroe, it's necessary to chose the method including "median" or "mean" of those genes to dealt with them. The mean, means to average those genes expression; while, the median refers to replace them with median of those genes expression.
}
\item{sftEstimate}{
The default soft threshold, when there no suitable threshold.
}
}

\examples{
wgcnaPipeline(GeneMatrix = "GSE_series_matrix.txt.gz",probeAnnotation = "probe_annotation.tsv",traits = "numeric_labels.csv")
}

\references{
Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 \n
G Yu, LG Wang, Y Han, QY He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology 2012, 16(5):284-287

}
