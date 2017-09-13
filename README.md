# CompactSeparableGraph

Implemented a compact representation of separable graphs in C++, which supported linear space to store the graph and constant time on graph queries (adjacency, degree, neighborhood queries), using METIS graph partitioning API and succinct representation of bit vectors.

The implementation reduces space usage by almost an order of magnitude, while supporting bread-first-search in acceptable running time when compared to the array-based representation.
