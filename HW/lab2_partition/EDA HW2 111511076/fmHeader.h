#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <climits>
#include <algorithm>
#include <functional>
#include <numeric>
#include <random>

using namespace std;

// Forward declarations
struct Node;
struct Net;
struct Bucket;

class Partitioner {
private:
    // Data structures
    struct Node {
        int id;
        int group; // 0 or 1
        int gain = 0; // New field
        vector<int> connected_nets;
    };

    struct Net {
        int id;
        vector<int> connected_nodes;
        int count_in_group[2] = {0,0};
    };

    // Bucket structure for FM algorithm
    // Optimized Bucket: array-based singly-linked lists per gain index
    class Bucket {
    private:
        int max_gain;
        int num_buckets;
        vector<int> head;           // head[node_gain + max_gain] -> node_id or -1
        vector<int> next_node;      // next node in bucket for each node_id
        vector<int> node_bucket;    // current bucket index for each node_id
        vector<int> node_gain;      // current gain for each node_id

    public:
        // Constructor: need max_gain and total number of nodes
        Bucket(int max_gain_val, int num_nodes)
            : max_gain(max_gain_val),
                num_buckets(2 * max_gain_val + 1),
                head(num_buckets, -1),
                next_node(num_nodes, -1),
                node_bucket(num_nodes, -1),
                node_gain(num_nodes, INT_MIN)
        {}

        void insert(int node_id, int gain) {
            gain = max(-max_gain, min(max_gain, gain)); //caps the gain to prevent accessing out of range
            int idx = gain + max_gain;
            node_gain[node_id] = gain;
            node_bucket[node_id] = idx;
            next_node[node_id] = head[idx];
            head[idx] = node_id;
        }

        void remove(int node_id) {
            int idx = node_bucket[node_id];
            if (idx < 0) return;
            int curr = head[idx], prev = -1;
            while (curr != -1 && curr != node_id) {
                prev = curr;
                curr = next_node[curr];
            }
            if (curr == -1) return; // not found
            if (prev == -1) head[idx] = next_node[curr];
            else next_node[prev] = next_node[curr];
            node_bucket[node_id] = -1;
            next_node[node_id] = -1;
        }

        void updateGain(int node_id, int new_gain) {
            remove(node_id);
            insert(node_id, new_gain);
        }

        int getMaxGainNode() const {
            for (int idx = num_buckets - 1; idx >= 0; --idx) {
                if (head[idx] != -1)
                    return head[idx];
            }
            return -1;
        }

        int getGain(int node_id) const {
            return node_gain[node_id];
        }
    };

    // Member variables
    vector<Node> nodes;
    vector<Net> nets;
    vector<Node> bestNodes;
    vector<Net> bestNets;
    vector<bool> isCutFlag;
    int cut_count, numOfNets, numOfNodes, maxGain = 0;
    int upperLimit, lowerLimit;
    int group_sizes[2];

    // Helper methods
    bool is_cut(const Net& net) const {
        return (net.count_in_group[0] > 0 && net.count_in_group[1] > 0);
    }

    int compute_incremental_gain_change(int node_id, int net_id) { //Used during the iteration
        const Net& net = nets[net_id];
        const Node& node = nodes[node_id];
        int group = node.group;
        int other_group = 1 - group;
    
        // Before move
        bool was_cut = (net.count_in_group[0] > 0 && net.count_in_group[1] > 0);
        int before_from = net.count_in_group[group];
        int before_to = net.count_in_group[other_group];
    
        // Simulate one move of 'node' to other group
        int after_from = before_from - 1;
        int after_to = before_to + 1;
        bool is_now_cut = (after_from > 0 && after_to > 0);
    
        if (!was_cut && is_now_cut) return -1;      //Gain -1 because moving led to increase in cut
        else if (was_cut && !is_now_cut) return 1;  //Gain +1 because moving led to decrease in cut
        else if (was_cut && is_now_cut) {
            if (before_from == 1) return 1;
            //Gain +1 because old FROM only had one node(this node)
            //Moving it decreases cut by one
            else if (before_to == 0 && after_to == 1) return -1;          
            //Gain -1 because old TO had zero node(no cut)
            //Moving it leads to increase in cut
        }
        return 0;
    }

public:
    // Constructor
    Partitioner(const string& filename) : cut_count(0)
    {
        group_sizes[0] = group_sizes[1] = 0;
        parse_input(filename);
        move_half_to_other_group_initialization();
    }

    // Input parsing
    void parse_input(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            return;
        }
        string str;
        if (getline(file, str)) {
            stringstream ss(str);
            ss >> numOfNets >> numOfNodes;
            
            nets.resize(numOfNets);
            nodes.resize(numOfNodes);

            // Initialize nets and nodes
            for (int i = 0; i < numOfNets; ++i) {
                //nets[i].connected_nodes.reserve(50); // Somehow not reserving runs faster
                nets[i].id = i;
            }
            // Initialize node IDs first
            for (int i = 0; i < numOfNodes; ++i) {
                nodes[i].id = i;
                nodes[i].group = 0; 
                // Assign all nodes to Group 0 initially
            }

            group_sizes[0] = numOfNodes;
            upperLimit = ceil(0.55 * double(numOfNodes));
            lowerLimit = floor(0.45 * double(numOfNodes));

            // Read net connections
            int currentNetIdx = 0;
            int nodeId;
            while (getline(file, str)) {
                stringstream ss(str);
                while (ss >> nodeId) {
                    nodeId -= 1;
                    nets[currentNetIdx].connected_nodes.push_back(nodeId);
                    nets[currentNetIdx].count_in_group[0]++;
                    // Becasue we assign all nodes to Group 0, 
                    // there's no need for the following if condition

                    // if (nodes[nodeId].group == 0)
                    //     nets[currentNetIdx].count_in_group[0]++;
                    // else 
                    //     nets[currentNetIdx].count_in_group[1]++;
                    
                    nodes[nodeId].connected_nets.push_back(currentNetIdx);
                    if (nodes[nodeId].connected_nets.size() > maxGain){
                        maxGain = nodes[nodeId].connected_nets.size();
                    }
                }
                cut_count += is_cut(nets[currentNetIdx]);
                ++currentNetIdx;
            }
            isCutFlag.resize(numOfNets);
            for (int i = 0; i < numOfNets; ++i) {
                isCutFlag[i] = (nets[i].count_in_group[0] > 0 && nets[i].count_in_group[1] > 0);
            }
            cut_count = accumulate(isCutFlag.begin(), isCutFlag.end(), 0);
        } else {
            cerr << "Error: File is empty" << endl;
        }
    }

    void move_half_to_other_group_initialization(){
        Bucket bucket(maxGain, numOfNodes);    // All nodes are in Group 0 initially 
        vector<bool> locked(numOfNodes, false); 
        for (auto& node : nodes) {
            node.gain = calculate_gain(node.id);
            bucket.insert(node.id, node.gain);
        }

        for (size_t i = 0; i < numOfNodes / 2; ++i) {
            int selected_node = -1;
            int group = -1;
            
            auto node = bucket.getMaxGainNode();
            if (!locked[node]) selected_node = node;

            // Bucket returns -1 if empty

            if (selected_node == -1) break;
            
            // Move node, lock the moved node, and remove it from the bucket
            move_node(selected_node);
            locked[selected_node] = 1;
            bucket.remove(selected_node);

            // Update neighbors of moved node
            const auto& moved_node = nodes[selected_node];
            for (int net_id : moved_node.connected_nets) {
                const Net& net = nets[net_id];
                for (int neighbor_id : net.connected_nodes) {
                    if (locked[neighbor_id]) continue;
                    Node& neighbor = nodes[neighbor_id];
                    int delta = compute_incremental_gain_change(neighbor_id, net_id);
                    if (delta != 0) {
                        neighbor.gain += delta;
                        bucket.updateGain(neighbor.id, neighbor.gain);
                    }
                }
            }
        }

        int cut = 0;

        for (auto& net : nets) {
            net.count_in_group[0] = net.count_in_group[1] = 0;
        } //Reset Group 0 & 1 Count

        for (int net_id = 0; net_id < numOfNets; ++net_id) {
            for (int node_id : nets[net_id].connected_nodes) {
                nets[net_id].count_in_group[nodes[node_id].group]++;
            } //Recount Group 0 & 1 Elements
            if (is_cut(nets[net_id])) cut++; //Recalc cut_count
        }
        cut_count = cut;
    }

    //Partitioning operations
    void move_node(int node_id) {
        Node& node = nodes[node_id];
        int old_group = node.group;
        int new_group = 1 - old_group;
        node.group = new_group;
        group_sizes[old_group]--;
        group_sizes[new_group]++;

        for (int net_id : node.connected_nets) {
            Net& net = nets[net_id];
            bool was_cut = isCutFlag[net_id];
            
            net.count_in_group[old_group]--;
            net.count_in_group[new_group]++;
            bool is_now_cut = (net.count_in_group[0] > 0 && net.count_in_group[1] > 0);
            if (was_cut != is_now_cut) {
                isCutFlag[net_id] = is_now_cut;
                cut_count += is_now_cut ? +1 : -1;
            }
        }
    }

    void revert_node(int node_id) { //Essentially just move_node without the gain update for neighbors
        Node& node = nodes[node_id];
        int old_group = node.group;
        int new_group = 1 - old_group;
        node.group = new_group;
        group_sizes[old_group]--;
        group_sizes[new_group]++;
        for (int net_id : node.connected_nets) {
            Net& net = nets[net_id];
            net.count_in_group[old_group]--;
            net.count_in_group[new_group]++;
        }
    }

    int calculate_gain(int node_id) const { //Only used at the very beginning 
        int gain = 0;
        const Node& node = nodes[node_id];
        int from = node.group;
        int to = 1 - from;

        for (int net_id : node.connected_nets) {
            const Net& net = nets[net_id];
            bool was_cut = is_cut(net);

            int from_count = net.count_in_group[from];
            int to_count = net.count_in_group[to];

            from_count--;
            to_count++;

            bool is_now_cut = (from_count > 0 && to_count > 0);

            if (!was_cut && is_now_cut) gain -= 1;
            else if (was_cut && !is_now_cut) gain += 1;
            else if (was_cut && is_now_cut) {
                if (net.count_in_group[from] == 1) gain += 1;
                else if (net.count_in_group[to] == 1) gain -= 1;
            }
        }
        return gain;
    }

    void run_fm_pass() {
        Bucket bucket0(maxGain, numOfNodes), bucket1(maxGain, numOfNodes);
        vector<bool> locked(numOfNodes, false);
        vector<int> move_sequence;    
        move_sequence.reserve(nodes.size());

        int cut = 0;
        for (auto& net : nets) {
            net.count_in_group[0] = net.count_in_group[1] = 0;
        } //Reset Group 0 & 1 Count

        for (int net_id = 0; net_id < numOfNets; ++net_id) {
            for (int node_id : nets[net_id].connected_nodes) {
                nets[net_id].count_in_group[nodes[node_id].group]++;
            } //Recount Group 0 & 1 Elements
            if (is_cut(nets[net_id])) cut++; //Recalc cut_count
        }
        cut_count = cut;


        int best_cut = cut_count;
        int best_index = -1;
        int tolerance = 0;
        if (numOfNodes > 10000)
            tolerance = 10;
        else
            tolerance = 100;
        int badMoves = 0;

        for (auto& node : nodes) {
            node.gain = calculate_gain(node.id);
            if (node.group == 0) bucket0.insert(node.id, node.gain);
            else                 bucket1.insert(node.id, node.gain);
        }

        for (size_t i = 0; i < nodes.size(); ++i) {
            int selected_node = -1;
            int group = -1;

            int node0 = bucket0.getMaxGainNode();
            int node1 = bucket1.getMaxGainNode();
            int gain0 = (node0 != -1 && !locked[node0]) ? bucket0.getGain(node0) : INT_MIN;
            int gain1 = (node1 != -1 && !locked[node1]) ? bucket1.getGain(node1) : INT_MIN;

            if (gain0 == INT_MIN && gain1 == INT_MIN) break;

            int selected_gain;

            // Check if balance is outside allowed range to force selection from larger group
            if (group_sizes[0] > upperLimit) {
                // Must move from group0 to reduce its size
                selected_node = node0;
                selected_gain = gain0;
                group = 0;
            } else if (group_sizes[0] < lowerLimit) {
                // Must move from group1 to reduce its size
                selected_node = node1;
                selected_gain = gain1;
                group = 1;
            } else {
                // Balance is within allowed range: select max gain node
                if (gain0 > gain1) {
                    selected_node = node0;
                    selected_gain = gain0;
                    group = 0;
                } else if (gain1 > gain0) {
                    selected_node = node1;
                    selected_gain = gain1;
                    group = 1;
                } else {
                    // Gains are equal; choose the group with larger size to balance
                    if (group_sizes[0] > group_sizes[1]) {
                        selected_node = node0;
                        selected_gain = gain0;
                        group = 0;
                    } else {
                        selected_node = node1;
                        selected_gain = gain1;
                        group = 1;
                    }
                }
            }

            if (selected_node == -1) break;

            if (selected_gain < 0){
                ++badMoves;     //Increment the counter if the gain is negative
            } else {
                badMoves = 0;   //Reset badMoves counter if we get positive gain
            }

            if (badMoves == tolerance) break;

            move_node(selected_node);
            locked[selected_node] = 1;
            if (group == 0) bucket0.remove(selected_node);
            else            bucket1.remove(selected_node);

            move_sequence.push_back(selected_node);

            bool is_balanced = true;
            if (numOfNodes <= 20) {
                // Allow difference of at most 1 node
                is_balanced = abs(group_sizes[0] - group_sizes[1]) <= 1;
            } else {
                is_balanced = (group_sizes[0] > lowerLimit && group_sizes[0] < upperLimit);
            }
            
            if (cut_count < best_cut && is_balanced) {
                best_cut = cut_count;
                best_index = move_sequence.size() - 1;
            }

            const auto& moved_node = nodes[selected_node];
            for (int net_id : moved_node.connected_nets) {
                const Net& net = nets[net_id];
                for (int neighbor_id : net.connected_nodes) {
                    if (locked[neighbor_id]) continue;
                    Node& neighbor = nodes[neighbor_id];
                    int delta = compute_incremental_gain_change(neighbor_id, net_id);
                    if (delta != 0) {
                        neighbor.gain += delta;
                        if (neighbor.group == 0)
                            bucket0.updateGain(neighbor.id, neighbor.gain);
                        else
                            bucket1.updateGain(neighbor.id, neighbor.gain);
                    }
                }
            }
        }

        if (best_index == -1) return;

        for (size_t i = move_sequence.size() - 1; i > best_index; --i) {
            revert_node(move_sequence[i]); //Revert the moves
        }
        return;
    }

    void print_status(string outputFile) {
        ofstream outfile(outputFile);
        for (const Node& node : nodes) {
            outfile << node.group << "\n";
        }
    }
    
    // Accessors
    int get_cut_count() const { return cut_count; }
    const vector<Node>& get_nodes() const { return nodes; }
    const vector<Net>& get_nets() const { return nets; }
    int get_group1count() const { return group_sizes[1]; }
    int get_group0count() const { return group_sizes[0]; }
};
