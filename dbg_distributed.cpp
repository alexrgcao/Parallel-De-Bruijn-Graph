#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <bits/stdc++.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cassert>
#include <mpi.h>

#define DEFAULT_FILE_NAME "NONE"
#define DEFAULT_NUM_KMER 10

// Serial Version of Debruijn Graph and Eulerian Path are inspired from https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master
class Node {
public:
    std::string km1mer;
    int nin;
    int nout;

    Node(std::string km1mer) : km1mer(km1mer), nin(0), nout(0) {}

    bool isSemiBalanced() const {
        return abs(nin - nout) == 1;
    }

    bool isBalanced() const {
        return nin == nout;
    }
};

std::vector<char> stringToChar(const std::string& str) {
    return std::vector<char>(str.begin(), str.end());
}

std::string charToString(const std::vector<char>& charVec) {
    return std::string(charVec.begin(), charVec.end());
}

class DeBruijnGraph {
private:
    std::unordered_map<std::string, Node*> nodes;
    std::unordered_map<Node*, std::vector<Node*>> G;
    int nsemi, nbal, nneither;
    Node* head;
    Node* tail;

public:
    /**
     * Constructs a DeBruijnGraph from a list of k-mers.
     * 
     * @param combinedKmersChunks A string containing all k-mers separated by spaces.
     * @param k The length of the k-mers.
     * @param size The number of processes to use for parallelization.
     * @param rank The rank of the current process.
     * @param circularize A boolean indicating whether to circularize the k-mers. Default is false.
     * 
     * Each of the Rank is responsible for producing its own copies of the Graph.
     * For each k-mer, it extracts the left and right (k-1)-mers, and creates or retrieves the corresponding nodes.
     * It updates the in-degree and out-degree of the nodes, and adds an edge from the left node to the right node.
     * After the entire graph has been constructed, it updates the counts of balanced, semi-balanced, and neither nodes, 
     * and sets the head and tail of the graph if it has an Eulerian walk.
     */
    DeBruijnGraph(char* combinedKmersChunks, int k, int size, int rank, bool circularize = false) 
        : nsemi(0), nbal(0), nneither(0), head(nullptr), tail(nullptr) {
        std::istringstream iss(combinedKmersChunks);
        std::string outerKmer;
        int n = 0;
        while (std::getline(iss, outerKmer, ' ')) {
            n++;
            int kmers_length = outerKmer.size();
            if (circularize) {
                outerKmer += outerKmer.substr(0, k-1);
                kmers_length += k-1;
            }
            for (int i = 0; i < kmers_length - (k-1); ++i) {
                std::string kmer = outerKmer.substr(i, k);
                std::string km1L = outerKmer.substr(i, k-1);
                std::string km1R = outerKmer.substr(i + 1, k-1);

                Node* nodeL = nullptr;
                Node* nodeR = nullptr;
                {
                    if (nodes.count(km1L)) {
                        nodeL = nodes[km1L];
                    } else {
                        nodeL = new Node(km1L);
                        nodes[km1L] = nodeL;
                    }
                    if (nodes.count(km1R)) {
                        nodeR = nodes[km1R];
                    } else {
                        nodeR = new Node(km1R);
                        nodes[km1R] = nodeR;
                    }
                }
                nodeL->nout += 1;
                nodeR->nin += 1;
                {
                    G[nodeL].push_back(nodeR);
                }
            }
        }
        
        for (auto& kv : nodes) {
            Node* node = kv.second;
            if (node->isBalanced()) {
                nbal += 1;
            } else if (node->isSemiBalanced()) {
                if (node->nin == node->nout + 1) {
                    tail = node;
                }
                if (node->nin == node->nout - 1) {
                    head = node;
                }
                nsemi += 1;
            } else {
                nneither += 1;
            }
        }
    }

    int nnodes() const {
        return nodes.size();
    }

    int nedges() const {
        return G.size();
    }

    bool hasEulerianWalk() const {
        return nneither == 0 && nsemi == 2;
    }

    bool hasEulerianCycle() const {
        return nneither == 0 && nsemi == 0;
    }

    bool isEulerian() const {
        return hasEulerianWalk() || hasEulerianCycle();
    }

    void exportToDot(const std::string& filename) const {
        std::ofstream file(filename);
        file << "digraph DeBruijnGraph {\n";
        for (const auto& kv : G) {
            Node* src = kv.first;
            for (Node* dst : kv.second) {
                file << "    \"" << src->km1mer << "\" -> \"" << dst->km1mer << "\";\n";
            }
        }
        file << "}\n";
    }

    /**
     * Finds an Eulerian walk or cycle in the graph using Hierholzer's algorithm and parallelizes the process using multiple MPI processes.
     * 
     * This method uses depth-first search (DFS) to find an Eulerian walk or cycle in the graph. 
     * It divides the graph into chunks, and each process processes one chunk of the graph.
     * Each process starts at a random head node and starts performing DFS. 
     * If a process hits the head node of another process, it joins and that process starts DFS from a new head node.
     * 
     * @param rank The rank of the current process.
     * @param size The number of processes to use for parallelization.
     * @return A vector of strings representing the Eulerian walk or cycle.
     * 
     */
    std::vector<std::string> eulerianWalkOrCycle(int rank, int size) {
        std::unordered_map<std::string, std::vector<std::string>> g;
        for (const auto& kv : G) {
            Node* src = kv.first;
            if (src == nullptr) continue;
            for (Node* dst : kv.second) {
                g[src->km1mer].push_back(dst->km1mer);
            }
        }

        if (g.empty()) return {};

        std::string src = g.begin()->first;

        if (isEulerian() && hasEulerianWalk()) {
            g[tail->km1mer].push_back(head->km1mer);
        }

        std::unordered_set<std::string> visited;
        std::vector<std::vector<std::string>> paths(size);
        std::vector<std::stack<std::string>> stacks(size);

        std::vector<int> stop(size, 0);
        std::vector<std::string> headNodes(size);

        MPI_Barrier(MPI_COMM_WORLD);

        std::string node;
        while (true) {
            bool found = false;

            for (const auto& kv : g) {
                if (visited.count(kv.first) == 0) {
                    node = kv.first;
                    visited.insert(node);
                    found = true;
                    break;
                }
            }

            if (!found) {
                break;
            }

            MPI_Barrier(MPI_COMM_WORLD);

            std::vector<char> node_char = stringToChar(node);
            MPI_Bcast(node_char.data(), node_char.size(), MPI_CHAR, rank, MPI_COMM_WORLD);
            std::string received_node = charToString(node_char);
            headNodes[rank] = received_node;
            dfs_parallel(g, node, visited, paths[rank], stacks[rank], rank, stop, headNodes, paths);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < size; ++i) {
            if (i != rank) {
                MPI_Status status;
                int flag;
                MPI_Iprobe(i, 0, MPI_COMM_WORLD, &flag, &status);
                if (flag) {
                    std::vector<std::string> path_i;
                    int count;
                    MPI_Get_count(&status, MPI_CHAR, &count);
                    path_i.resize(count);
                    MPI_Recv(&path_i[0], count, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    for (const auto& node : path_i) {
                        paths[rank].push_back(node);
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        return paths[rank];
    }

    /**
     * Performs a depth-first search (DFS) on the graph in parallel using MPI.
     * 
     * This method is called by each MPI process to process its chunk of the graph. It uses DFS to find a path in the graph.
     * If the current node is the head node of another process's path, it finds a new node that has not been visited yet, 
     * sets it as the new head node for the current process, and merges the paths of the two processes.
     * 
     * @param g The graph.
     * @param node The starting node of the DFS.
     * @param visited A reference to a set of visited nodes.
     * @param tour A reference to a vector where the path found by the DFS will be stored.
     * @param stack A reference to a stack used by the DFS.
     * @param rank The rank of the current process.
     * @param stop A reference to a vector of booleans indicating whether each process should stop.
     * @param headNodes A reference to a vector of the head nodes of each process's path.
     * @param paths A reference to a vector of the paths found by each process.
     * 
     */
    void dfs_parallel(std::unordered_map<std::string, std::vector<std::string>>& g, std::string node, std::unordered_set<std::string>& visited, std::vector<std::string>& tour, std::stack<std::string>& stack, int rank, std::vector<int>& stop, std::vector<std::string>& headNodes, std::vector<std::vector<std::string>>& paths) {
        stack.push(node);
        while (!stack.empty() && !stop[rank]) {
            node = stack.top();
            stack.pop();

            if (g.count(node) > 0) {
                for (int k = 0; k < g[node].size(); ++k) {
                    tour.push_back(node);
                }

                for (auto it = g[node].rbegin(); it != g[node].rend(); ++it) {
                    if (visited.count(*it) == 0) {
                        visited.insert(*it);
                        stack.push(*it);

                        MPI_Barrier(MPI_COMM_WORLD);
                        std::vector<char> it_char = stringToChar(*it);
                        MPI_Bcast(it_char.data(), it_char.size(), MPI_CHAR, rank, MPI_COMM_WORLD);
                    }
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);

            for (int j = 0; j < headNodes.size(); ++j) {
                bool signal = 0;
                int received_values[3] = {0};
                if (j != rank && node == headNodes[j]) {
                    signal = 1;
                    int values[3] = {signal, j, rank};

                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Bcast(values, 3, MPI_INT, rank, MPI_COMM_WORLD);

                    std::string newNode;
                    for (const auto& kv : g) {
                        if (visited.count(kv.first) == 0) {
                            newNode = kv.first;
                            break;
                        }
                    }
                    if (newNode.empty()) {
                        stop[rank] = 1;
                        break;
                    }
                    headNodes[rank] = newNode;
                    visited.insert(newNode);

                    std::vector<char> headNodes_char = stringToChar(headNodes[rank]);
                    MPI_Bcast(headNodes_char.data(), headNodes_char.size(), MPI_CHAR, rank, MPI_COMM_WORLD);
                    std::vector<char> visited_char = stringToChar(newNode);
                    MPI_Bcast(visited_char.data(), visited_char.size(), MPI_CHAR, rank, MPI_COMM_WORLD);

                    paths[j].insert(paths[j].begin(), paths[rank].begin(), paths[rank].end());
                    paths[rank].clear();
                    paths[rank].push_back(newNode);
                    stop[rank] = 0;
                    break;
                }
                if (received_values[0] == 1) {
                    int received_j = received_values[1];
                    int received_rank = received_values[2];
                    stop[received_values[2]] = 1;
                    paths[received_j].insert(paths[received_j].begin(), paths[received_rank].begin(), paths[received_rank].end());
                    paths[received_rank].clear();
                    break;
                }
            }
            if (stop[rank]) {
                break;
            }
        }
    }
};



char* getKmersMPI(std::vector<char*>& seqChunks, int k, int size, int rank) {
    
    std::string kmersChunks;
    for (const auto& seq : seqChunks) {
        int seq_length = strlen(seq);
        for (int l = 0; l <= seq_length - k; ++l) {
            char* kmer = new char[k + 1];
            std::strncpy(kmer, seq + l, k);
            kmer[k] = '\0';
            kmersChunks += kmer;
            kmersChunks += ' ';
            delete[] kmer;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    char* serializedKmers = new char[kmersChunks.size() + 1];
    std::strcpy(serializedKmers, kmersChunks.c_str());
    return serializedKmers;
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
    cxxopts::Options options("De Bruijn Graph Distributed",
                             "Distributed implementation of De Bruijn Graph algorithm");
    options.add_options(
        "custom",
        {
            {"f,file", "File name for the dot file",
            cxxopts::value<std::string>()->default_value(DEFAULT_FILE_NAME)},
            {"k,kmer", "Number of kmer for the graph",
            cxxopts::value<int>()->default_value(std::to_string(DEFAULT_NUM_KMER))},
        }
    );
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    auto cl_options = options.parse(argc, argv);
    std::string filename = cl_options["file"].as<std::string>();
    int kmer = cl_options["kmer"].as<int>();

    if (filename == DEFAULT_FILE_NAME) {
        std::cerr << "Please provide a file name" << std::endl;
        return 1;
    }
    if (kmer <= 1) {
        std::cerr << "Please provide a kmer bigger than 1" << std::endl;
        return 1;
    }

    std::vector<char*> sequences;

    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    std::streamsize total_size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::streamsize local_start = world_rank * total_size / world_size;
    std::streamsize local_end = (world_rank + 1) * total_size / world_size;

    if (world_rank != 0) {
        file.seekg(local_start);
        std::string temp;
        std::getline(file, temp);
        local_start = file.tellg();
    }
    if (world_rank != world_size - 1) {
        file.seekg(local_end);
        std::string temp;
        std::getline(file, temp);
        local_end = file.tellg();
    } else {
        local_end = total_size;
    }

    std::streamsize local_size = local_end - local_start;
    file.seekg(local_start);
    std::vector<char> local_lines(local_size);
    file.read(local_lines.data(), local_size);

    bool is_sequence_line = false;
    std::istringstream iss(std::string(local_lines.begin(), local_lines.end()));
    std::string line, sequence;

    bool is_fastq = false;
    if (filename.substr(filename.find_last_of(".") + 1) == "fastq") {
        is_fastq = true;
    }
            
    while (std::getline(iss, line)) {
        if (line.empty()) {
            continue;
        }
        if (!is_fastq) {
            if (line[0] == '>' || line[0] == '@') { // New sequence in FASTA
                if (!sequence.empty()) {
                    char* sequence_cstr = new char[sequence.length() + 1];
                    std::strcpy(sequence_cstr, sequence.c_str());
                    sequences.push_back(sequence_cstr);
                    sequence.clear();
                }
                is_sequence_line = true;
                continue;
            }
        } else {
            if (line[0] == '@') { // New sequence in FASTQ
                if (!sequence.empty()) {
                    char* sequence_cstr = new char[sequence.length() + 1];
                    std::strcpy(sequence_cstr, sequence.c_str());
                    sequences.push_back(sequence_cstr);
                    sequence.clear();
                }
                is_sequence_line = true;
                continue;
            } else if (line[0] == '+') { // Quality score line in FASTQ
                is_sequence_line = false;
                continue;
            }
        }
        if (is_sequence_line) { // Sequence line in FASTQ or FASTA
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        char* sequence_cstr = new char[sequence.length() + 1];
        std::strcpy(sequence_cstr, sequence.c_str());
        sequences.push_back(sequence_cstr);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    timer tKmers, tGraph, tEulerian, tTotal;
    double timeKmers = 0.0;
    double timeGraph = 0.0;
    double timeEulerian = 0.0;
    double timeTotal = 0.0;
    tTotal.start();
    tKmers.start();
    char* kmersChunks = getKmersMPI(sequences, kmer, world_size, world_rank);
    
    int length = strlen(kmersChunks) + 1;

    int* lengths = new int[world_size];
    MPI_Allgather(&length, 1, MPI_INT, lengths, 1, MPI_INT, MPI_COMM_WORLD);

    int* displacements = new int[world_size];
    int total_length = 0;
    for (int i = 0; i < world_size; ++i) {
        displacements[i] = total_length;
        total_length += lengths[i];
    }

    char* allKmersChunks = new char[total_length];
    allKmersChunks[total_length - 1] = '\0'; // null-terminate the string

    MPI_Allgatherv(kmersChunks, length, MPI_CHAR, allKmersChunks, lengths, displacements, MPI_CHAR, MPI_COMM_WORLD);

    for (int i = 0; i < total_length - 1; ++i) {
        if (allKmersChunks[i] == '\0') {
            allKmersChunks[i] = ' ';
        }
    }
    timeKmers = tKmers.stop();
    MPI_Barrier(MPI_COMM_WORLD);
    tGraph.start();

    DeBruijnGraph graph(allKmersChunks, kmer, world_size, world_rank);
    delete[] kmersChunks;

    if (world_rank == 0) {
        printf("Number of nodes: %d\n", graph.nnodes());
        printf("Number of edges: %d\n", graph.nedges());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    timeGraph = tGraph.stop();
    MPI_Barrier(MPI_COMM_WORLD);
    tEulerian.start();
    if (graph.isEulerian()) {
        std::vector<std::string> longestPath = graph.eulerianWalkOrCycle(world_rank, world_size);
        if (world_rank == 0) {
            printf("The graph is Eulerian\n");
            printf("Size of longestPath: %lu\n", longestPath.size());
            std::string losequence = longestPath[0];
            for (int i = 1; i < longestPath.size(); ++i) {
                losequence += longestPath[i].back();
            }
            printf("Sequence: %s\n", losequence.c_str());
        }
    } else {
        std::vector<std::string> longestPath = graph.eulerianWalkOrCycle(world_rank, world_size);
        if (world_rank == 0) {
            printf("The graph is not Eulerian\n");
            std::string losequence = longestPath[0];
            for (int i = 1; i < longestPath.size(); ++i) {
                losequence += longestPath[i].back();
            }
            printf("Sequence: %s\n", losequence.c_str());
        }
    }
    timeEulerian = tEulerian.stop();
    timeTotal = tTotal.stop();

    for (char* sequence : sequences) {
        delete[] sequence;
    }
    
    if (world_rank == 0) {
        printf("Rank, Kmer Time, Graph Time, Eulerian Time, Total Time\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("%d, %f, %f, %f, %f\n", world_rank, timeKmers, timeGraph, timeEulerian, timeTotal);

    MPI_Finalize();
    return 0;
}

