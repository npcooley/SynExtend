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
ExoLabel(edgelistfiles, outfile=tempfile(),
              mode=c("undirected", "directed"),
              add_self_loops=FALSE,
              ignore_weights=FALSE,
              normalize_weights=FALSE,
              iterations=0L, inflation=1.05,
              return_table=FALSE,
              consensus_cluster=FALSE,
              verbose=interactive(),
              sep='\t',
              tempfiledir=tempdir(),
              cleanup_files=TRUE)
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
Should self loops be added to the network? If \code{TRUE}, adds self loops of weight 1.0 to all vertices. If set to numeric value \code{w}, adds self loops of weight \code{w} to all nodes. Note that this takes place prior to edge weight normalization (if requested).
}

\item{ignore_weights}{
Should weights be ignored? If \code{TRUE}, all edges will be treated as an edge of weight 1. Must be set to \code{TRUE} if any of \code{edgelistfiles} are two-column tables (start->end only, lacking a weights column).
}

\item{normalize_weights}{
Should weights be normalized? If \code{TRUE}, each vertex's edge weights are normalized such that the sum of outgoing edge weights is 1. This normalization is done after adding self loops.
}

\item{iterations}{
Number of iterations to run fast label propagation algorithm for. Set to a value of 0 or less for infinite iterations.
}

\item{inflation}{
Inflation parameter for edges. Every other iteration, all edges are raised to the power of \code{inflation} and then renormalized. Higher values speed up algorithm convergence but produce smaller clusters. Defaults to 1.05; set to 1.0 to disable inflation.
}

\item{return_table}{
Should result of clustering be returned as a file, or a \code{data.frame} object? If \code{FALSE}, returns a character vector corresponding to the path of \code{outfile}. If \code{TRUE}, parses \code{outfile} using \code{\link{read.table}} and returns the result. Not recommended for very large graphs.
}

\item{consensus_cluster}{
Should consensus clustering be used? If \code{TRUE}, runs the clustering algorithm nine times and forms a consensus clustering based on the agreement of each run. Can be set to a double to control the number of iterations, see Details for more information.
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
\item{cleanup_files}{
Should intermediary files be deleted when the process completes? Note that \code{outfile} will only be deleted if \code{return_table=TRUE} AND \code{cleanup_files=TRUE}.
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

Consensus clustering can be enabled by setting \code{consensus_cluster=TRUE}. Consensus clustering runs this algorithm on the graph multiple times, transforming weight values according to a sigmoid function. By default, this runs nine times for sigmoids with scale 0.5 and shape \code{c(0,0.2,0.4,0.6,0.8,1.0,1.33,1.67,2.0)}, collapsing weights below 0.1 to zero. The resulting clusters form a network such that the edge weight between any two nodes connected in the initial graph is the proportion of clusters they shared over clustering runs. This network is used for a final label propagation run, which identifies the consensus clusters. Users can specify a numeric vector as input to \code{consensus_cluster}, which will override the default shape parameters and number of iterations.
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

  Specific calculations for memory/disk consumption are detailed below. In summary, the absolute worst case memory consumption is roughly \eqn{(24l+16)v} bytes, and the disk consumption during computation is \eqn{(41+12d)v} bytes. The final table returned consumes \eqn{(2+l+\log_{10}{v})v} bytes.

  ExoLabel builds a trie to keep track of vertex names. Each internal node of the trie consumes 24 bytes, and each leaf node consumes 16 bytes. The lowest possible RAM consumption of the trie (if every label is length \eqn{l} and shares the same prefix of length \eqn{l-1}) is roughly \eqn{40v} bytes, and the maximum RAM consumption (if no two node labels have any prefix in common) is \eqn{(24l + 16)v} bytes. We can generalize this to estimate the total memory consumption as roughly \eqn{(24(l-p)+16)v}, where \eqn{p} is the average length of common prefix between any two node labels.

  ExoLabel also uses a number of internal caches to speed up read/writes from files. These caches take less than 100MB of RAM in total.

  As for disk space, ExoLabel produces five total files during its runtime: a cluster file, a CSR-compressed network, two queues, a bitmap, and two temp files for sorting. The two queues together contain no more than a total of \eqn{v} entries of 8 bytes each, and the bitmap contains \eqn{v} entries of 1 byte each. The cluster file also contains \eqn{v} entries of 8 bytes each, and the two sorting files are at most each the same size as the cluster file. The CSR file contains a header of \eqn{v+1} entries of 8 bytes, and then one 12 byte entry per outgoing edge. The number of edges to record is \eqn{vd}. Thus, the total disk consumption in bytes is \eqn{8v + v + 3*8v + 8(v+1) + 12vd \approx (41+12d)v}.

  The final table returned is a tab-separated table containing vertex names and cluster numbers in human-readable format. Each line consumes at most \eqn{l + 2 + \log_{10}{v}} bytes. In the worst case, the number of clusters is equal to the number of vertices, which have \eqn{\log_{10}{v}} digits. The average number of digits is close to the number of digits of the largest number due to how number of digits scale with numbers. The extra two bytes are for the tab and newline characters. Thus, the total size of the file is at most \eqn{(2+l+\log_{10}{v})v} bytes.

  Average node outdegree can be difficult to predict. Many biological networks are approximately scale-free, meaning their degree distribution roughly follows a power law with exponent in the range 2-3. The average of a general power law distribution is not always well-defined, but for the purposes of the below example I'll conservatively assume a minimum value of 2 and shape parameter of 2, leading to an average node outdegree of 4 if the degree is distributed according to a Pareto distribution.

  What does this look like in practice? For a set of one million unique nodes and average label length 10, the RAM consumption is \code{140-356MB} (note an extra \code{100MB} is added due to cache sizes). Sequence similarity networks often have  The intermediate files consume a total of \code{89MB}, and the final file consumes roughly \code{19MB}. At one billion unique nodes and the same other parameters, the RAM consumption is \code{40-256GB}, the intermediate files consume \code{89GB}, and the final file is roughly \code{19GB} (\code{108GB} total). This sounds like a lot, but consider that an input file with two decimal places for each edge weight must necessarily be at least \eqn{(2l+3+4)vd} bytes (each of the \eqn{vd} edges contains two vertex labels, two tabs, newline, and 4 characters for weight). Thus the input file is \code{108GB}, which is equivalent to the total disk space consumed by ExoLabel.

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
