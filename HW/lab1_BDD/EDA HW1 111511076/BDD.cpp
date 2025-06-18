#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <memory>
#include <vector>
#include <limits>
#include "ReadTXT.h"
#include "SortVariablesByOrder.h"
#include "BoolExprTranslate.h"
#include "BoolExprEval.h"

using namespace std;

// ------------------------ BDD Node ----------------------------
struct BDDNode {
    char var; // Variable name (ex: 'a', 'b')
    shared_ptr<BDDNode> low;  // Pointer when var = 0
    shared_ptr<BDDNode> high; // Pointer when var = 1
    int value = -1; // Terminal node: 0 or 1 (non-terminal = -1)

    BDDNode(char v, shared_ptr<BDDNode> l, shared_ptr<BDDNode> h)
        : var(v), low(l), high(h) {}

    BDDNode(int val) : value(val) {} // Terminal node

    bool is_terminal(){
        return (value == -1) ? false : true;
    }
};

// ------------------------ Unique Table for Reduction ----------------------------
struct NodeHasher {
    size_t operator()(const tuple<char, shared_ptr<BDDNode>, shared_ptr<BDDNode>>& t) const {
        return hash<char>()(get<0>(t)) ^
               hash<shared_ptr<BDDNode>>()(get<1>(t)) ^
               hash<shared_ptr<BDDNode>>()(get<2>(t));
    }
};

// ------------------------ Build BDD Recursively ----------------------------
shared_ptr<BDDNode> buildBDD(BooleanExpression &expr, const string vars, int index, vector<bool>& inputs,
    unordered_map<tuple<char, shared_ptr<BDDNode>, shared_ptr<BDDNode>>, shared_ptr<BDDNode>, NodeHasher>& uniqueTable,
    const shared_ptr<BDDNode>& terminal0, const shared_ptr<BDDNode>& terminal1
) {
    if (index == vars.size()) {
        return expr.evaluate(inputs) ? terminal1 : terminal0;  
    }

    char var = vars[index];
    // Low branch: var = 0
    inputs[index] = false;
    auto low = buildBDD(expr, vars, index + 1, inputs, uniqueTable, terminal0, terminal1);

    // High branch: var = 1
    inputs[index] = true;
    auto high = buildBDD(expr, vars, index + 1, inputs, uniqueTable, terminal0, terminal1);

    // Reduce: if low == high, skip this node
    if (low == high) return low;

    auto key = make_tuple(var, low, high);
    if (uniqueTable.count(key)) {
        return uniqueTable[key];
    }

    auto node = make_shared<BDDNode>(var, low, high);
    uniqueTable[key] = node;
    return node;
}

// ------------------------ Count BDD Nodes ----------------------------
void count_nodes_helper(const shared_ptr<BDDNode>& node, unordered_set<BDDNode*>& visited) {
    if (!node || node->is_terminal()) return;
    if (visited.insert(node.get()).second) {
        count_nodes_helper(node->low, visited);
        count_nodes_helper(node->high, visited);
    }
}

int count_nodes(const shared_ptr<BDDNode>& node) {
    unordered_set<BDDNode*> visited;
    count_nodes_helper(node, visited);
    return visited.size() + 2;
    //+2 for 0 and 1 terminals
}

// ------------------------ Print BDD Nodes for Debugging ----------------

void printBDD(shared_ptr<BDDNode> node, int depth = 0) {
    if (!node) return;

    for (int i = 0; i < depth; ++i) cout << "  ";

    if (node->value != -1) {
        cout << "Leaf: " << node->value << endl;
        return;
    }

    cout << "Var: " << node->var << endl;

    for (int i = 0; i < depth; ++i) cout << "  ";
    cout << " 0 → ";
    printBDD(node->low, depth + 1);

    for (int i = 0; i < depth; ++i) cout << "  ";
    cout << " 1 → ";
    printBDD(node->high, depth + 1);
}


// ------------------------ Main ----------------------------
int main(int argc, char* argv[]) {    
    string boolExpr;
    vector<string> orders;
    readFileAndParse(argv[1], boolExpr, orders);
    int varCount = orders[0].size();
    int bestAns = 2147483647;
    for (string order : orders){
        vector<bool> inputs(varCount, false);
    
        // Fresh terminals for this BDD
        auto terminal0 = make_shared<BDDNode>(0);
        auto terminal1 = make_shared<BDDNode>(1);
    
        // Fresh unique table for this order
        unordered_map<tuple<char, shared_ptr<BDDNode>, shared_ptr<BDDNode>>, shared_ptr<BDDNode>, NodeHasher> uniqueTable;
    
        // Sorted & translated expression
        string sortedExprCpp = convertToCppLogic(sortExprByOrder(boolExpr, order));
        BooleanExpression expr(sortedExprCpp, order);
    
        // Build and count this BDD
        auto root = buildBDD(expr, order, 0, inputs, uniqueTable, terminal0, terminal1);
        int localAns = count_nodes(root);
        if (localAns < bestAns) bestAns = localAns;
    }
    ofstream outfile(argv[2]);
    outfile << bestAns;
    outfile.close();
    return 0 ;
}
