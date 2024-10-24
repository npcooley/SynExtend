\name{ExoLabel}
\alias{ExoLabel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ExoLabel: Out of Memory Fast Label Propagation
}
\description{
Runs Fast Label Propagation using disk space for constant memory complexity.
}
\usage{
ExoLabel(edgelistfiles,
              outfile=tempfile(),
              mode=c("undirected", "directed"),
              add_self_loops=FALSE,
              ignore_weights=FALSE,
              normalize_weights=FALSE,
              iterations=0L,
              inflation=1.05,
              return_table=FALSE,
              consensus_cluster=FALSE,
              verbose=interactive(),
              sep='\t',
              tempfiledir=tempdir())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{edgelistfiles}{
Character vector of files to be processed. Each entry should be a machine-interpretable path to an edgelist file. See Details for expected format.
}

\item{outfile}{
File to write final clusters to. Optional, defaults to a temporary file.
}

\item{mode}{
String specifying whether edges should be interpreted as undirected (default) or directed. Can be "undirected", "directed", or an unambiguous abbreviation.
}

\item{add_self_loops}{
Should self loops be added to the network? If \code{TRUE}, adds self loops of weight 1.0 to all vertices. If set to numeric value \code{w}, adds self loops of weight \code{w} to all nodes.
}

\item{ignore_weights}{
Should weights be ignored? If \code{TRUE}, all edges will be treated as an edge of weight 1. Must be set to \code{TRUE} if any of \code{edgelistfiles} are two-column tables (start->end only, lacking a weights column).
}

\item{normalize_weights}{
Should weights be normalized? If \code{TRUE}, each vertex's edge weights are normalized such that the sum of outgoing edge weights is 1. This normalization is done after adding self loops.
}

\item{iterations}{
Maximum number of times to process each node. If set to zero or \code{NULL}, automatically uses the square root of the max node degree. See "Algorithm Convergence" for more information.
}

\item{inflation}{
Inflation parameter for edges. See "Algorithm Convergence" below for a description of this parameter. Higher values speed up algorithm convergence but produce smaller clusters. Defaults to 1.05; set to 1.0 to disable inflation.
}

\item{return_table}{
Should result of clustering be returned as a file, or a \code{data.frame} object? If \code{FALSE}, returns a character vector corresponding to the path of \code{outfile}. If \code{TRUE}, parses \code{outfile} using \code{\link{read.table}} and returns the result. Not recommended for very large graphs.
}

\item{consensus_cluster}{
Should consensus clustering be used? If \code{TRUE}, runs the clustering algorithm multiple times and forms a consensus clustering based on the agreement of each run. Can be set to a vector of doubles to control the number of iterations. See "Consensus Clustering" below for more information.
}

\item{verbose}{
Should status messages (output, progress, etc.) be displayed while running?
}

\item{sep}{
Character that separates entries on a line in each file in \code{edgelistfiles}. Defaults to tab, as would be expected in a \code{.tsv} formatted file. Set to \code{','} for a \code{.csv} file.
}
\item{tempfiledir}{
Character vector corresponding to the location where temporary files used during execution should be stored. Defaults to R's \code{\link{tempdir}}.
}
}
\details{
Very large graphs require too much RAM for processing on some machines. In a graph containing billions of nodes and edges, loading the entire structure into RAM is rarely feasible. This implementation uses disk space for storing representations of each graph. While this is slower than computing on RAM, it allows this algorithm to scale to graphs of enormous size while only using a comparatively small amount of memory. See "Memory Consumption" for details on the total disk/memory consumption of ExoLabel.

This function expects a set of edgelist files, provided as a vector of filepaths. Each entry in the file is expected to be in the following:

\code{VERTEX1<sep>VERTEX2<sep>WEIGHT<linesep>}

This line defines a single edge between vertices \code{VERTEX1} and \code{VERTEX2} with weight \code{WEIGHT}. \code{VERTEX1} and \code{VERTEX2} are strings corresponding to vertex names, \code{WEIGHT} is a numeric value that can be interpreted as a \code{double}. The separators \code{<sep>} and \code{<linesep>} correspond to the arguments \code{sep} and \code{linesep}, respectively. The default arguments work for standard \code{.tsv} formatting, i.e., a file of three columns of tab-separated values.

If \code{ignore_weight=TRUE}, the file can be formatted as:

\code{VERTEX1<sep>VERTEX2<linesep>}

Note that the \code{v1 v2 w} format is still accepted for \code{ignore_weight=FALSE}, but the specified weights will be ignored.
}

\section{Algorithm Convergence}{
One of the main issues of Label Propagation algorithms is that they can fail to converge. Consider an unweighted directed graph with four nodes connected in a loop. That is, \code{A->B, B->C, C->D, D->A}. If \code{A,C} are in cluster 1 and \code{B,D} are in cluster 2, this algorithm could keep processing all the nodes in a loop and never converge. To solve this issue, we introduce two measures for convergence: \code{inflation} and \code{iterations}.

\code{iterations} is the simpler parameter to understand. If \code{iterations=x}, then we only allow the algorithm to process each node \code{x} times. Once a given node has been seen \code{x} times, it is no longer updated.

\code{inflation} gradually increases the influence of stronger weighted edges as we see a node more. In other words, the more often we see a node, the more bias we add towards its strongest weighted edges. For each node, we use the following weighting:

\eqn{w' = w^{1+log_2(n-1)}}

Here \eqn{w} is the original edge weight, \eqn{w'} is the new edge weight, and \eqn{n} is the number of times we've already processed the node. After this transformation, the edge weights are renormalized, meaning that large weights tend to get larger, and small weights tend to get smaller. Logarithms prevent the exponents from growing too large, and base 2 is chosen for computational efficiency.
}

\section{Consensus Clustering}{
Consensus clustering can be enabled by setting \code{consensus_cluster=TRUE}. Consensus clustering runs ExoLabel on the input graph multiple times, transforming weight values according to a sigmoid function. By default, this runs nine times for sigmoids with scale 0.5 and shape \code{c(0,0.2,0.4,0.6,0.8,1.0,1.33,1.67,2.0)}, collapsing weights below 0.1 to zero. The resulting clusters form a network such that the edge weight between any two nodes connected in the initial graph is the proportion of clusters they shared over clustering runs. This network is used for a final label propagation run, which identifies the consensus clusters. Users can specify a numeric vector as input to \code{consensus_cluster}, which will override the default shape parameters and number of iterations.
}

\value{
If \code{return_table=TRUE}, returns a \code{\link{data.frame}} object with two columns. The first column contains the name of each vertex, and the second column contains the cluster it was assigned to.

If \code{return_table=FALSE}, returns a character vector of length 1. This vector contains the path to the file where clusters were written to. The file is formatted as a \code{.tsv}, with each line containing two tab separated columns (vertex name, assigned cluster)
}

\references{
Traag, V.A., Subelj, L. Large network community detection by fast label propagation. \emph{Sci Rep} \bold{13}, 2701 (2023). https://doi.org/10.1038/s41598-023-29610-z
}
\author{
Aidan Lakshman <AHL27@pitt.edu>
}
\section{Warning}{
While this algorithm can scale to very large graphs, it does have some internal limitations. First, nodes must be comprised of no more than 254 characters. If this limitation is restrictive, please feel free to contact me. Alternatively, you can increase the size yourself by changing the definition of \code{MAX_NODE_NAME_SIZE} in \code{src/outmem_graph.c}. This limitation is provided to decrease memory overhead and improve runtime, but arbitrary values are possible.

Second, nodes are indexed using 64-bit unsigned integers, with 0 reserved for other values. This means that the maximum possible number of nodes available is 2^64-2, which is about 18.5 quintillion.

Third, this algorithm uses disk space to store large objects. As such, please ensure you have sufficient disk space for the graph you intend to process. I've tried to put safeguards in the code itself, but funky stuff can happen when the OS runs out of space. See "Memory Consumption" for details on the total disk/memory consumption of ExoLabel.
}

\section{Memory Consumption}{
  Let \eqn{v} be the number of unique nodes, \eqn{d} the average outdegree of nodes, and \eqn{l} the average length of node labels.

  Specific calculations for memory/disk consumption are detailed below. In summary, the absolute worst case memory consumption is roughly \eqn{(24l+34)v} bytes, and the disk consumption during computation is \eqn{(8+12d)v} bytes. The final table returned consumes \eqn{(2+l+\log_{10}{v})v} bytes.

  ExoLabel builds a trie to keep track of vertex names. Each internal node of the trie consumes 24 bytes, and each leaf node consumes 16 bytes. The lowest possible RAM consumption of the trie (if every label is length \eqn{l} and shares the same prefix of length \eqn{l-1}) is roughly \eqn{40v} bytes, and the maximum RAM consumption (if no two node labels have any prefix in common) is \eqn{(24l + 16)v} bytes. We can generalize this to estimate the total memory consumption as roughly \eqn{(24(l-p)+16)v}, where \eqn{p} is the average length of common prefix between any two node labels.

  ExoLabel also uses a number of internal caches to speed up read/writes from files. These caches take less than 100MB of RAM in total. It also uses an internal queue for processing nodes, which consumes roughly \eqn{10v} bytes. It also uses an internal index of size \eqn{8v} bytes.

  As for disk space, ExoLabel transforms the graph into a CSR-compressed network, which is split across three files: a header, a neighbors list, and a weights list. The header file contains \eqn{v+1} entries of 8 bytes, and the other two files consume a total of 12 bytes  per outgoing edge. The number of edges to record is \eqn{vd}. Thus, the total disk consumption in bytes is \eqn{8(v+1) + 12vd \approx (8+12d)v}.

  The final table returned is a tab-separated table containing vertex names and cluster numbers in human-readable format. Each line consumes at most \eqn{l + 2 + \log_{10}{v}} bytes. In the worst case, the number of clusters is equal to the number of vertices, which have \eqn{\log_{10}{v}} digits. The average number of digits is close to the number of digits of the largest number due to how number of digits scale with numbers. The extra two bytes are for the tab and newline characters. Thus, the total size of the file is at most \eqn{(2+l+\log_{10}{v})v} bytes.
}

\examples{
num_verts <- 20L
num_edges <- 20L
all_verts <- sample(letters, num_verts)
all_edges <- vapply(seq_len(num_edges),
      \(i) paste(c(sample(all_verts, 2L),
                   as.character(round(runif(1),3))),
                 collapse='\t'),
                    character(1L))
edgefile <- tempfile()
if(file.exists(edgefile)) file.remove(edgefile)
writeLines(all_edges, edgefile)
res <- ExoLabel(edgefile, return_table=TRUE)
print(res)
}

\seealso{
\code{\link{EstimateExoLabel}}
}
