# Parallel-De-Bruijn-Graph
This project is designed to provide practical experience with parallel programming, using C++ Threads and MPI. It is for the project of the course CMPT 431 (Parallel and Distributed Systems).


## Table of Contents
1. [Introduction](#intro)

2. [Installation](#installation)

3. [Reproducing this project](#repro)

4. [Acknowledgements](#ack)

<a name="intro"></a>
## 1. Introduction

One of the most effective and widely used approaches for the assembly of genomic sequences from short reads is the construction and traversal of De Bruijn graphs.

This project explores the application of De Bruijn graphs for genomic assembly, focusing on enhancing the efficiency of the algorithm through the implementation of threading and distributing.
Given the computational complexity and the massive volume of data involved in sequence assembly, a key objective of this project was to improve the performance and scalability of the De Bruijn graph approach through parallelization.

To explore the challenges, this project employs parallelism in various part of De Bruijn Graph construction. 

  - The dividing of the initial set of kmers (substrings of length k) derived from fasta or fastq files, are processed by separate thread.
  -  Construction of graph are divided into segments to perform operations like nodes insertion, and edge creation.
  -  Graph transversal is parallelized by letting each process to start DFS at a random node denoted as head node, if a process hits a head node of another process, the two path joins and the initial process starts DFS on a new random node.

### What to find where

Explain briefly what files are found where

```bash
repository
├── core                         ## Folder with required library
├── Makefile                     ## Compile the program
├── dbg_serial.cpp               ## Serial Version of De Bruijn Graph
├── dbg_parallel.cpp             ## Parallel Version of De Bruijn Graph using C++ Threads
├── dbg_distributed.cpp          ## Distributed Version of De Bruijn Graph using MPI
├── README.md                    ## You are here
```

<a name="installation"></a>
## 2. Installation

Ensure C++ compiler is installed

```bash
git clone $THISREPO
make
    
```

<a name="repro"></a>
## 3. Reproduction
```bash
./dbg_serial -f PATH_TO_FILE -k KMER_SIZE
./dbg_parallel -f PATH_TO_FILE -k KMER_SIZE -n NUMBER_OF_THREADS
mpirun -np NUMBER_OF_SYSTEM ./dbg_distributed -f PATH_TO_FILE -k KMER_SIZE
```

Output will be displayed in the terminal with information regarding graph size, if eulerian or not, sequence, and time performance

<a name="ack"></a>
## 4. Acknowledgement
The serial version of De Bruijn Graph implementation is inspired from https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_deBruijn.ipynb
