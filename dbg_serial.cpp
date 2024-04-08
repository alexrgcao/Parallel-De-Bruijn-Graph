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

#define DEFAULT_FILE_NAME "NONE"
#define DEFAULT_NUM_KMER 10

// All implementation of Debruijn Graph and Eulerian Path are inspired from https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_deBruijn.ipynb
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

class DeBruijnGraph {
private:
    std::unordered_map<std::string, Node*> nodes;
    std::unordered_map<Node*, std::vector<Node*>> G;
    int nsemi, nbal, nneither;
    Node* head;
    Node* tail;

public:
    /**
     * Construct a De Bruijn graph from a collection of strings.
     * 
     * @param strIter A vector of strings from which the De Bruijn graph is constructed.
     * @param k The length of the k-mers used to construct the De Bruijn graph.
     * @param circularize A boolean flag indicating whether to circularize the input strings. 
     *        If true, the function appends the first k-1 characters of each string to its end. 
     *        This is useful for constructing De Bruijn graphs of circular sequences. Default is false.
     * 
     * The constructor initializes the counters for the number of balanced nodes (nbal), 
     * semi-balanced nodes (nsemi), neither balanced nor semi-balanced nodes (nneither), 
     * head and tail pointers.
     * 
     * For each string in strIter, it generates all k-mers and (k-1)-mers. 
     * For each k-mer, it creates two nodes corresponding to the left and right (k-1)-mers 
     * and adds an edge from the left node to the right node in the graph.
     * 
     * After constructing the graph, it iterates over all nodes and checks their balance. 
     * If a node is balanced, it increments nbal. 
     * If a node is semi-balanced, it increments nsemi and sets the head and tail pointers if necessary. 
     * If a node is neither balanced nor semi-balanced, it increments nneither.
     */
    DeBruijnGraph(std::vector<std::string> strIter, int k, bool circularize = false) 
        : nsemi(0), nbal(0), nneither(0), head(nullptr), tail(nullptr) {
        for (auto& st : strIter) {
            if (circularize) {
                st += st.substr(0, k-1);
            }
            for (int i = 0; i < st.size() - (k-1); ++i) {
                std::string kmer = st.substr(i, k);
                std::string km1L = st.substr(i, k-1);
                std::string km1R = st.substr(i+1, k-1);
                Node* nodeL = nullptr;
                Node* nodeR = nullptr;
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
                nodeL->nout += 1;
                nodeR->nin += 1;
                G[nodeL].push_back(nodeR);
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
     * Recursive function to find an Eulerian tour in the graph.
     * 
     * @param kmer A reference to the adjacency list of the graph.
     * @param node The current node in the tour.
     * @param tour A reference to the vector storing the tour.
     * 
     * This function uses Hierholzer's algorithm, which is a depth-first search (DFS) algorithm, 
     * to find an Eulerian tour in the graph. 
     * 
     * The algorithm works as follows:
     * 1. Starting from the given node, follow edges until you come back to the same node, removing each edge as you follow it. 
     *    This is the DFS part of the algorithm.
     * 2. Add the node to the tour.
     * 3. If there are still edges left from the node, continue the process from the new node. 
     *    This ensures that all edges are followed exactly once.
     */
    std::vector<std::string> eulerianWalkOrCycle() {
        std::unordered_map<std::string, std::vector<std::string>> g;
        printf("Size of G: %ld\n", G.size());
        for (const auto& kv : G) {
            Node* src = kv.first;
            for (Node* dst : kv.second) {
                g[src->km1mer].push_back(dst->km1mer);
            }
        }

        std::string src = g.begin()->first;
        std::vector<std::string> tour;

        if (isEulerian()) {
            if (hasEulerianWalk()) {
                g[tail->km1mer].push_back(head->km1mer);
            }
        } else {
            src = g.begin()->first;
        }

        euler_r(g, src, tour);

        std::reverse(tour.begin(), tour.end());
        tour.pop_back();

        if (isEulerian() && hasEulerianWalk()) {
            auto it = std::find(tour.begin(), tour.end(), head->km1mer);
            std::rotate(tour.begin(), it, tour.end());
        }

        return tour;
    }
    
    void euler_r(std::unordered_map<std::string, std::vector<std::string>>& kmer, std::string node, std::vector<std::string>& tour) {
        while (!kmer[node].empty()) {
            std::string next_node = kmer[node].back();
            kmer[node].pop_back();
            euler_r(kmer, next_node, tour);
        }
        tour.push_back(node);
    }
};

int main(int argc, char *argv[]) {
  // Initialize command line arguments
    cxxopts::Options options("DFS_serial",
                             "serial implementation of depth first search algorithm");
    options.add_options(
        "custom",
        {
            {"f,file", "File name for the dot file",
            cxxopts::value<std::string>()->default_value(DEFAULT_FILE_NAME)},
            {"k,kmer", "Number of kmer for the graph",
            cxxopts::value<int>()->default_value(std::to_string(DEFAULT_NUM_KMER))},
            
        }
    );

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

    std::vector<std::string> sequences;
    std::ifstream file(filename);
    std::string line, sequence;
    bool is_sequence_line = false;

    bool is_fastq = false;
    if (filename.substr(filename.find_last_of(".") + 1) == "fastq" || filename.substr(filename.find_last_of(".") + 1) == "fq") {
        is_fastq = true;
    }

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        if (!is_fastq) {
            if (line[0] == '>' || line[0] == '@') {
                if (!sequence.empty()) {
                    sequences.push_back(sequence);
                    sequence.clear();
                }
                is_sequence_line = true;
                continue;
            }
        } else {
            if (line[0] == '@') {
                if (!sequence.empty()) {
                    sequences.push_back(sequence);
                    sequence.clear();
                }
                is_sequence_line = true;
                continue;
            } else if (line[0] == '+') {
                is_sequence_line = false;
                continue;
            }
        }
        if (is_sequence_line) {
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }

    timer totalTimer, graphKmerTimer, eulerianTimer;
    double total_time = 0.0;
    double graph_kmer_time = 0.0;
    double eulerian_time = 0.0;
    totalTimer.start();
    graphKmerTimer.start();
    DeBruijnGraph graph(sequences, kmer);
    graphKmerTimer.stop();

    std::cout << "Number of nodes: " << graph.nnodes() << std::endl;
    std::cout << "Number of edges: " << graph.nedges() << std::endl;
    
    if (graph.isEulerian()) {
        eulerianTimer.start();
        std::cout << "The graph is Eulerian" << std::endl;
        std::vector<std::string> eulerianPath = graph.eulerianWalkOrCycle();
        eulerian_time = eulerianTimer.stop();
        std::string sequence = eulerianPath[0];
        for (size_t i = 1; i < eulerianPath.size(); ++i) {
            sequence += eulerianPath[i].back();
        }
        std::cout << sequence << std::endl;
    } else {
        eulerianTimer.start();
        std::cout << "The graph is not Eulerian" << std::endl;
        std::vector<std::string> eulerianPath = graph.eulerianWalkOrCycle();
        eulerian_time = eulerianTimer.stop();
        std::string sequence = eulerianPath[0];
        for (size_t i = 1; i < eulerianPath.size(); ++i) {
            sequence += eulerianPath[i].back();
        }
        std::cout << sequence << std::endl;
    }
    total_time = totalTimer.stop();
    std::cout << "Kmer graph construction time: " << graph_kmer_time << " seconds" << std::endl;
    std::cout << "Eulerian path construction time: " << eulerian_time << " seconds" << std::endl;
    std::cout << "Total time: " << total_time << " seconds" << std::endl;
    
    return 0;
}
