# Approximate graph embedding

Simple and hopefully fast 2d embeddings for networks

## Compilation

Compile with mpicc

Suggested compilation command:

mpicc embedding.c -O3 -o embedding

## Running

The code can be run with

`mpirun -np 8 ./embedding -i 600 -t 10 -n 1000 graph.edgelist > out.txt`

to run on 8 threads for 600 iterations on a graph with 1000 nodes, with edges given in the file `graph.edgelis`. Intermediate output will be produced to standard out every $t = 10$ iterations.

Calling options:

| flag | description |
|------|-------------|
| FILE | Input edge list file. Required. |
| i    | Number of iterations. Defaults to 1000 |
| k    | Number of closest and furthest neighbors to consider. Defaults to approx. log(n). |
| n    | Highest node id in input. Required. |
| p    | Penalty for following reverse edges. Defaults to 3. |
| t    | Output universe state every T iterations. |
| x    | Number of rows in universe. Defaults to approx. √n. |
| y    | Number of columns in universe. Defaults to approx. √n. |
| ?    | Outputs help on the command line. |

The **edge list file** should contain on edge per line with the source id and target id separated by a tab character. An example `edgelist` file can be found in the source.

The input graph needs to have exactly **n** nodes numbered from 1 to n (inclusive). For reliable results, the graph should be connected.

## Demo with a small graph over 500 iterations

http://saskeli.kapsi.fi/1000_embed.mp4
