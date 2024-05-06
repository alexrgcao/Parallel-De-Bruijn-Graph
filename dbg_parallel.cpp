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

#define DEFAULT_FILE_NAME "NONE"
#define DEFAULT_NUM_KMER 10
#define DEFAULT_NUMBER_OF_THREADS 4

// Serial version of the DeBruijnGraph class and the eulerianWalkOrCycle method inspired from https://colab.research.google.com/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_deBruijn.ipynb
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

/**
 * ConcurrentQueue is a thread-safe queue data structure.
 * 
 * This class is a template class, meaning it can hold elements of any type.
 * It uses a std::queue to hold the elements, and a std::mutex to ensure thread safety.
 * 
 * @tparam T The type of elements in the queue.
 */
template<typename T>
class ConcurrentQueue {
private:
    std::queue<T> queue;
    mutable std::mutex mutex;

public:
    /**
     * Adds an item to the end of the queue.
     * 
     * This function is thread-safe because it locks the mutex before modifying the queue.
     * 
     * @param item The item to add to the queue.
     */
    void push(const T& item) {
        std::lock_guard<std::mutex> lock(mutex);
        queue.push(item);
    }

    /**
     * Attempts to remove and return the item at the front of the queue.
     * 
     * This function is thread-safe because it locks the mutex before modifying the queue.
     * 
     * @param item A reference where the front item will be stored if the queue is not empty.
     * @return true if an item was successfully popped, false if the queue was empty.
     */
    bool try_pop(T& item) {
        std::lock_guard<std::mutex> lock(mutex);
        if (queue.empty()) {
            return false;
        }
        item = queue.front();
        queue.pop();
        return true;
    }

    /**
     * Checks if the queue is empty.
     * 
     * This function is thread-safe because it locks the mutex before accessing the queue.
     * 
     * @return true if the queue is empty, false otherwise.
     */
    bool empty() {
        std::lock_guard<std::mutex> lock(mutex);
        return queue.empty();
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
     * Constructs a DeBruijnGraph from a list of k-mers,
     * uses multithreading to speed up the construction of the graph.
     * 
     * @param kmersChunks A vector of vectors of strings, where each inner vector is a chunk of k-mers.
     * @param k The length of the k-mers.
     * @param nThreads The number of threads to use for parallelization.
     * @param times A reference to a vector of doubles where the time taken by each thread will be stored.
     * @param circularize A boolean indicating whether to circularize the k-mers. Default is false.
     * 
     * It divides the list of k-mers into chunks, and each thread processes one chunk of k-mers.
     * For each k-mer, extract the left and right (k-1)-mers, and create or retrieve the corresponding nodes.
     * Update the in-degree and out-degree of the nodes, and add an edge from the left node to the right node.
     * The access to shared resources (the nodes and the adjacency list) is protected by the mutex.
     * The constructor uses a mutex to synchronize access to shared resources, ensuring thread safety.
     * Update the counts of balanced, semi-balanced, and neither nodes, and set the head and tail of the graph if it has an Eulerian walk.
     */
    DeBruijnGraph(std::vector<std::vector<std::string>> seqChunks, int k, int nThreads, std::vector<double>& times, bool circularize = false) 
        : nsemi(0), nbal(0), nneither(0), head(nullptr), tail(nullptr) {
        std::vector<std::thread> threads(nThreads);
        std::mutex mtx; // Mutex for synchronizing access to shared resources
        std::vector<timer> timers(nThreads);

        for (int t = 0; t < nThreads; ++t) {
            threads[t] = std::thread([&, t]() {
                timers[t].start();
                for (auto& seq : seqChunks[t]) {
                    if (circularize) {
                        seq += seq.substr(0, k-1);
                    }
                    for (int i = 0; i < seq.size() - (k-1); ++i) {
                        std::string kmer = seq.substr(i, k);
                        std::string km1L = seq.substr(i, k-1);
                        std::string km1R = seq.substr(i+1, k-1);
                        Node* nodeL = nullptr;
                        Node* nodeR = nullptr;
                        {
                            std::lock_guard<std::mutex> lock(mtx); // Lock the mutex while accessing shared resources
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
                            std::lock_guard<std::mutex> lock(mtx); // Lock the mutex while accessing shared resources
                            G[nodeL].push_back(nodeR);
                        }
                    }
                }
                times[t] = timers[t].stop();
            });
        }
        for (auto& thread : threads) {
            thread.join();
        }

        // Update nsemi, nbal, and nneither after the entire graph has been constructed
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
     * Finds an Eulerian walk or cycle in the graph using Hierholzer's algorithm and parallelizes the process using multiple threads.
     * 
     * This method uses depth-first search (DFS) to find an Eulerian walk or cycle in the graph. 
     * It divides the graph into chunks, and each thread processes one chunk of the graph.
     * Each thread is set to starts at a random head node and start performing DFS. 
     * If a thread hits the head node of another thread, it joins and that thread starts at a new head node.
     * 
     * @param kmersChunks A vector of vectors of strings, where each inner vector is a chunk of k-mers.
     * @param nThreads The number of threads to use for parallelization.
     * @param times A reference to a vector of doubles where the time taken by each thread will be stored.
     * @return A vector of strings representing the Eulerian walk or cycle.
     */

    std::vector<std::string> eulerianWalkOrCycle(int nThreads, std::vector<double>& times) {
        std::unordered_map<std::string, std::vector<std::string>> g;
        std::vector<timer> timers(nThreads);

        for (const auto& kv : G) {
            Node* src = kv.first;
            if (src == nullptr) continue;
            for (Node* dst : kv.second) {
                g[src->km1mer].push_back(dst->km1mer);
            }
        }

        std::string src = g.begin()->first;

        if (isEulerian()) {
            if (hasEulerianWalk()) {
                g[tail->km1mer].push_back(head->km1mer);
            }
        } else {
            src = g.begin()->first;
        }

        std::unordered_set<std::string> visited;
        std::mutex visitedMutex;
        std::vector<std::thread> workers;
        std::vector<std::vector<std::string>> paths(nThreads);
        std::vector<std::mutex> pathMutexes(nThreads);

        std::vector<std::stack<std::string>> stacks(nThreads);

        std::vector<bool> stop(nThreads, false);
        std::vector<std::string> headNodes(nThreads);

        
        for (int i = 0; i < nThreads; ++i) {
            workers.emplace_back([&, i]() {
                timers[i].start();
                std::string node;
                while (true) {
                    bool found = false;
                    {
                        std::unique_lock<std::mutex> lock(visitedMutex);
                        for (const auto& kv : g) {
                            if (visited.count(kv.first) == 0) {
                                node = kv.first;
                                visited.insert(node);
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found) {
                        break;
                    }
                    headNodes[i] = node;
                    euler_parallel(g, node, visited, paths[i], visitedMutex, pathMutexes[i], stacks[i], i, stop, headNodes, paths);
                }
                times[i] = timers[i].stop();
            });
        }
        for(auto& t : workers) {
            t.join();
        }

        std::vector<std::string> tour;
        std::unordered_set<std::string> included;

        for (int i = 0; i < nThreads; ++i) {
            for (const auto& node : paths[i]) {
                tour.push_back(node);
                included.insert(node);
                
            }
        }

        return tour;
    }

    /**
     * Performs a depth-first search (DFS) on the graph in parallel.
     * 
     * This method is called by each thread to process its chunk of the graph. It uses DFS to find a path in the graph.
     * If the current node is the head node of another thread's path, it finds a new node that has not been visited yet, 
     * sets it as the new head node for the current thread, and merges the paths of the two threads.
     * 
     * @param g The graph.
     * @param node The starting node of the DFS.
     * @param visited A reference to a set of visited nodes.
     * @param tour A reference to a vector where the path found by the DFS will be stored.
     * @param visitedMutex A reference to a mutex for synchronizing access to the visited nodes.
     * @param tourMutex A reference to a mutex for synchronizing access to the path.
     * @param stack A reference to a stack used by the DFS.
     * @param i The index of the current thread.
     * @param stop A reference to a vector of booleans indicating whether each thread should stop.
     * @param headNodes A reference to a vector of the head nodes of each thread's path.
     * @param paths A reference to a vector of the paths found by each thread.
     */
    void euler_parallel(std::unordered_map<std::string, std::vector<std::string>>& g, const std::string& node, std::unordered_set<std::string>& visited, std::vector<std::string>& tour, std::mutex& visitedMutex, std::mutex& tourMutex, std::stack<std::string>& stack, int i, std::vector<bool>& stop, std::vector<std::string>& headNodes, std::vector<std::vector<std::string>>& paths) {
        stack.push(node);
        while (!stack.empty() && !stop[i]) {
            std::string node;
            {
                std::unique_lock<std::mutex> lock(visitedMutex);
                if (stack.empty()) continue;
                node = stack.top();
                stack.pop();
            }
            {
                std::unique_lock<std::mutex> lock(tourMutex);
                for (int k = 0; k < g[node].size(); ++k) {
                    tour.push_back(node);
                }
            }
            for (auto it = g[node].rbegin(); it != g[node].rend(); ++it) {
                std::unique_lock<std::mutex> lock(visitedMutex);
                if (visited.count(*it) == 0) {
                    visited.insert(*it);
                    stack.push(*it);
                }
                for (int j = 0; j < headNodes.size(); ++j) {
                    if (j != i && node == headNodes[j]) {
                        std::string newNode;
                        for (const auto& kv : g) {
                            if (visited.count(kv.first) == 0) {
                                newNode = kv.first;
                                break;
                            }
                        }
                        if (newNode.empty()) {
                            stop[i] = true;
                            break;
                        }
                        headNodes[i] = newNode;
                        paths[j].insert(paths[j].begin(), paths[i].begin(), paths[i].end());
                        paths[i].clear();
                        paths[i].push_back(newNode);
                        stack.push(newNode);
                        stop[i] = false;
                        break;
                    }
                }
                if (stop[i]) {
                    break;
                }
            }
        }
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
            {"nThreads", "Number of Threads",
            cxxopts::value<int>()->default_value(std::to_string(DEFAULT_NUMBER_OF_THREADS))},
        }
    );

    auto cl_options = options.parse(argc, argv);
    std::string filename = cl_options["file"].as<std::string>();
    int kmer = cl_options["kmer"].as<int>();
    int nThreads = cl_options["nThreads"].as<int>();
    std::vector<std::string> sequences;

    if (filename == DEFAULT_FILE_NAME) {
        std::cerr << "Please provide a file name" << std::endl;
        return 1;
    }
    if (kmer <= 1) {
        std::cerr << "Please provide a kmer bigger than 1" << std::endl;
        return 1;
    }
    if (nThreads <= 0) {
        std::cerr << "Please provide a number of threads bigger than 0" << std::endl;
        return 1;
    }

    std::ifstream file(filename);
    std::string line, sequence;
    bool is_sequence_line = false;

    bool is_fastq = false;
    if (filename.substr(filename.find_last_of(".") + 1) == "fastq") {
        is_fastq = true;
    }

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        if (!is_fastq) {
            if (line[0] == '>' || line[0] == '@') { // New sequence in FASTA
                if (!sequence.empty()) {
                    sequences.push_back(sequence);
                    sequence.clear();
                }
                is_sequence_line = true;
                continue;
            }
        } else {
            if (line[0] == '@') { // New sequence in FASTQ
                if (!sequence.empty()) {
                    sequences.push_back(sequence);
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
        sequences.push_back(sequence);
    }
    
    timer tGraph, tEulerian, tTotal;
    double totalGraph = 0.0;
    double totalEulerian = 0.0;
    double timeTotal = 0.0;
    tTotal.start();
    
    std::vector<double> graphTimes = {0.0, 0.0, 0.0, 0.0};
    int chunkSize = sequences.size() / nThreads;
    int remainder = sequences.size() % nThreads;
    std::vector<std::vector<std::string>> seqChunks(nThreads);
    for (int i = 0; i < nThreads; ++i) {
        int start = i * chunkSize + std::min(i, remainder);
        int end = start + chunkSize + (i < remainder);
        for (int j = start; j < end; ++j) {
            seqChunks[i].push_back(sequences[j]);
        }
    }
    tGraph.start();
    DeBruijnGraph graph(seqChunks, kmer, nThreads, graphTimes);
    totalGraph = tGraph.stop();
    std::cout << "Number of nodes: " << graph.nnodes() << std::endl;
    std::cout << "Number of edges: " << graph.nedges() << std::endl;

    std::vector<double> eulerTimes = {0.0, 0.0, 0.0, 0.0};
    std::vector<int> nodesAddedByProcess(nThreads, 0);
    tEulerian.start();
    if (graph.isEulerian()) {
        std::vector<std::string> longestPath = graph.eulerianWalkOrCycle(nThreads, eulerTimes);
        printf("Size of longest path: %ld\n", longestPath.size());
        std::cout << "The graph is Eulerian" << std::endl;
        std::string completeSequence = longestPath[0];
        for (size_t i = 1; i < longestPath.size(); ++i) {
            completeSequence += longestPath[i].back();
        }
        std::cout << completeSequence << std::endl;

    } else {
        std::cout << "The graph is not Eulerian" << std::endl;
    }
    totalEulerian = tEulerian.stop();
    timeTotal = tTotal.stop();

    std::cout << "Threads, Graph, Eulerian" << std::endl;
    for (int i = 0; i < nThreads; ++i) {
        std::cout << i << ", " << graphTimes[i] << ", " << eulerTimes[i] << std::endl;
    }
    std::cout << "Total, Total Graph, Total Eulerian" << std::endl;
    std::cout << timeTotal << ", " << totalGraph << ", " << totalEulerian << std::endl;

    return 0;
}
