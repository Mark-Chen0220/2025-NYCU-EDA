#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <queue>
#include <utility>
#include <limits>
#include <algorithm>
#include <chrono> // For timing
#include <numeric> // For accumulate
#include <random>

using namespace std;

struct Gcell{
    int topEdgeCapacity, bottomEdgeCapacity, leftEdgeCapacity, rightEdgeCapacity;
    double M1cost, M2cost;
    double cost;
    int reachedUsingMetal;
    int topEdgeNet, bottomEdgeNet, leftEdgeNet, rightEdgeNet;
    int coordX, coordY;
};

struct RoutingArea {
    int offsetX, offsetY, width, height;
    vector<vector<Gcell>> GcellGrid;
};

struct GridDimensions {
    int Xlength, Yheight;
};

struct Chip {
    int dx, dy;       // Chip offset relative to routing area
    int width, height; // Chip dimensions (unused in routing)
    vector<pair<int, pair<int, int>>> bumps; // Bumps: (id, (bx, by))
};

vector<int> parseIntegers(const string& line) {
    vector<int> result;
    stringstream ss(line);
    int num;
    while (ss >> num) result.push_back(num);
    return result;
}

double alpha, myBeta, myGamma, delta, viaCost, overflowCost = 0;
double overflowWeight = 1; 
int gridCols, gridRows;
RoutingArea ra;
GridDimensions grid;
vector<Chip> chips;
//Mark these as global variables

void parseGMP(string filename){
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file." << endl;
        return;
    }

    vector<string> lines;
    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        lines.push_back(line);
    }

    for (int i = 0; i < lines.size();) {
        string current = lines[i];
        if (current == ".ra") {
            vector<int> values = parseIntegers(lines[++i]);
            ra = {values[0], values[1], values[2], values[3]};
            i++;
        } else if (current == ".g") {
            vector<int> values = parseIntegers(lines[++i]);
            grid = {values[0], values[1]};
            i++;
        } else if (current == ".c") {
            Chip chip;
            vector<int> cValues = parseIntegers(lines[++i]);
            chip.dx = cValues[0];
            chip.dy = cValues[1];
            chip.width = cValues[2];
            chip.height = cValues[3];
            i++;
            // Parse bumps
            while (i < lines.size() && lines[i] != ".c" && lines[i] != ".ra" && lines[i] != ".g") {
                if (lines[i] == ".b") {
                    i++;
                    while (i < lines.size() && lines[i][0] != '.') {
                        vector<int> bValues = parseIntegers(lines[i]);
                        chip.bumps.push_back({bValues[0], {bValues[1], bValues[2]}});
                        i++;
                    }
                } else i++;
            }
            chips.push_back(chip);
        } else i++;
    }
    gridCols = ra.width / grid.Xlength;
    gridRows = ra.height / grid.Yheight;
    ra.GcellGrid = vector<vector<Gcell>>(gridCols, vector<Gcell>(gridRows));
    file.close();
}

void parseGCL(string filename){
    ifstream file(filename);
    string line;
    getline(file, line); //get .ec
    for (int y = 0; y < gridRows; y++){
        for (int x = 0; x < gridCols; x++){
            getline(file, line); 
            vector<int>values = parseIntegers(line);
            ra.GcellGrid[x][y].topEdgeCapacity = values[1];
            ra.GcellGrid[x][y].bottomEdgeCapacity = values[1];
            ra.GcellGrid[x][y].leftEdgeCapacity = values[0];
            ra.GcellGrid[x][y].rightEdgeCapacity = values[0];
        }
    }
    file.close();
}

void parseCST(string filename){
    ifstream file(filename);
    string line;
    stringstream ss;
    getline(file, line);
    ss << line;
    ss >> line; 
    ss >> alpha;
    ss.clear();
    getline(file, line);
    ss << line;
    ss >> line; 
    ss >> myBeta;
    ss.clear();
    getline(file, line);
    ss << line;
    ss >> line; 
    ss >> myGamma;
    ss.clear();
    getline(file, line);
    ss << line;
    ss >> line; 
    ss >> delta;
    ss.clear();
    getline(file, line);
    getline(file, line);
    ss << line;
    ss >> viaCost;
    ss.clear();
    getline(file, line);
    for (int y = 0; y < gridRows; y++){
        getline(file, line);
        ss << line;
        for (int x = 0; x < gridCols; x++){
            double cost;
            ss >> cost;
            ra.GcellGrid[x][y].M1cost = cost;
            if (cost > overflowCost) overflowCost = cost;
        }
        ss.clear();
    }
    getline(file, line);
    for (int y = 0; y < gridRows; y++){
        getline(file, line);
        ss << line;
        for (int x = 0; x < gridCols; x++){
            double cost;
            ss >> cost;
            ra.GcellGrid[x][y].coordX = x;
            ra.GcellGrid[x][y].coordY = y;
            ra.GcellGrid[x][y].M2cost = cost;
            if (cost > overflowCost) overflowCost = cost;
        }
        ss.clear();
    }
    file.close();
}

double gCostFunction(int dir, bool layerSwap, Gcell &parent, Gcell &cell){
    
    double cost = 0;

    if (layerSwap){
        cost += myGamma * ((parent.M1cost + parent.M2cost) / 2) + delta * viaCost;
        if (parent.reachedUsingMetal == 1) cost -= myGamma * parent.M1cost;
        else cost -= myGamma * parent.M2cost;
    }

    if (cell.reachedUsingMetal == 1) cost += myGamma * cell.M1cost;
    else cost += myGamma * cell.M2cost;

    switch(dir){
        case 0: //up
            if (parent.topEdgeNet  + 1 - parent.topEdgeCapacity  > 0) cost += myBeta * (1 + parent.topEdgeNet - parent.topEdgeCapacity)   * overflowCost * overflowWeight;
            if (cell.bottomEdgeNet + 1 - cell.bottomEdgeCapacity > 0) cost += myBeta * (1 + cell.bottomEdgeNet - cell.bottomEdgeCapacity) * overflowCost * overflowWeight;
            break;
        case 1: //down
            if (cell.topEdgeNet  + 1 - cell.topEdgeCapacity  > 0) cost += myBeta * (1 + cell.topEdgeNet - cell.topEdgeCapacity)                   * overflowCost * overflowWeight;
            if (parent.bottomEdgeNet + 1 - parent.bottomEdgeCapacity > 0) cost += myBeta * (1 + parent.bottomEdgeNet - parent.bottomEdgeCapacity) * overflowCost * overflowWeight;
            break;
        case 2: //left
            if (parent.leftEdgeNet   + 1 - parent.leftEdgeCapacity   > 0) cost += myBeta * (1 + parent.leftEdgeNet - parent.leftEdgeCapacity)     * overflowCost * overflowWeight;
            if (cell.rightEdgeNet  + 1 - cell.rightEdgeCapacity  > 0) cost += myBeta * (1 + cell.rightEdgeNet - cell.rightEdgeCapacity)   * overflowCost * overflowWeight;
            break;
        case 3: //right
            if (cell.leftEdgeNet   + 1 - cell.leftEdgeCapacity   > 0) cost += myBeta * (1 + cell.leftEdgeNet - cell.leftEdgeCapacity)     * overflowCost * overflowWeight;
            if (parent.rightEdgeNet  + 1 - parent.rightEdgeCapacity  > 0) cost += myBeta * (1 + parent.rightEdgeNet - parent.rightEdgeCapacity)   * overflowCost * overflowWeight;
            break;
    }
    
    return cost;
}

vector<vector<pair<int,int>>> performRouting(const vector<int>& bumpOrder, double& totalCostOut) {
    vector<vector<pair<int,int>>> solutions(bumpOrder.size());
    totalCostOut = 0;

    // Directions: up, down, left, right
    const int dirs[4][2] = {{0,1}, {0,-1}, {-1,0}, {1,0}};

    for (int i : bumpOrder) {
        // Compute start and target cell coordinates
        int startX = (chips[0].dx + chips[0].bumps[i].second.first) / grid.Xlength;
        int startY = (chips[0].dy + chips[0].bumps[i].second.second) / grid.Yheight;
        int targetX = (chips[1].dx + chips[1].bumps[i].second.first) / grid.Xlength;
        int targetY = (chips[1].dy + chips[1].bumps[i].second.second) / grid.Yheight;

        // Precompute Manhattan-based heuristic cost on grid cells
        // Manhattan distance without obstacles can be computed without BFS

        // A* data structures
        vector<vector<pair<int,int>>> cameFrom(gridCols, vector<pair<int,int>>(gridRows, {-1,-1}));
        vector<vector<bool>> closedSet(gridCols, vector<bool>(gridRows, false));
        const double INF = numeric_limits<double>::infinity();
        vector<vector<double>> gScore(gridCols, vector<double>(gridRows, INF));

        // Min-heap entry: (fScore, pointer to Gcell)
        struct CellNode {
            double f;
            Gcell* cell;
            int metal;
            bool operator>(const CellNode& other) const { return f > other.f; }
        };

        priority_queue<CellNode, vector<CellNode>, greater<CellNode>> openHeap;

        // Initialize start
        gScore[startX][startY] = 0;
        CellNode startingNode = {0, &ra.GcellGrid[startX][startY], 1 };
        openHeap.push(startingNode);
        bool found = false;

        while (!openHeap.empty()) {
            auto current = openHeap.top();
            openHeap.pop();
            int cx = current.cell->coordX;
            int cy = current.cell->coordY;
            int cm = current.metal;

            if (closedSet[cx][cy]) continue;
            closedSet[cx][cy] = true;

            if (cx == targetX && cy == targetY) {
                found = true;
                break;
            }

            // Explore neighbors
            for (int i = 0; i < 4; i++) {
                auto &d = dirs[i];
                int nx = cx + d[0], ny = cy + d[1];
                if (nx < 0 || nx >= gridCols || ny < 0 || ny >= gridRows) continue;
                if (closedSet[nx][ny]) continue;

                // Determine metal change and calculate cost
                bool switchMetal = false;
                int nextMetal = cm;
                if (cm == 1 && d[0] != 0) switchMetal = true;
                else if (cm == 2 && d[1] != 0) switchMetal = true;
                if (switchMetal){
                    if (cm == 1) nextMetal = 2;
                    else nextMetal = 1;
                }

                ra.GcellGrid[nx][ny].reachedUsingMetal = nextMetal;
                double nextMoveGCost = gScore[cx][cy] + gCostFunction(i, switchMetal, ra.GcellGrid[cx][cy], ra.GcellGrid[nx][ny]);
                double fScore = nextMoveGCost +  alpha * (abs(targetX - nx) * grid.Xlength + abs(targetY - ny) * grid.Yheight);

                //Gcost accurate modeling

                if (nextMoveGCost < gScore[nx][ny]) {
                    cameFrom[nx][ny] = {cx, cy};
                    gScore[nx][ny] = nextMoveGCost;

                    ra.GcellGrid[nx][ny].cost = fScore; //allows accurate summing                
                    openHeap.push({fScore, &ra.GcellGrid[nx][ny], nextMetal});
                }
            }
        }

        // Reconstruct path if found
        vector<pair<int,int>> path;
        if (found) {
            pair<int,int> cur = {targetX, targetY};
            while (!(cur.first == startX && cur.second == startY)) {
                
                path.push_back(cur);
                cur = cameFrom[cur.first][cur.second];
            }
            path.push_back({startX, startY});
            reverse(path.begin(), path.end());
        }

        for (const auto& [x, y] : path) {
            totalCostOut += ra.GcellGrid[x][y].cost;
        }
        //Add path to solutions
        solutions[i] = path;

        //Update nets in the grid
        for (int i = 1; i < path.size(); ++i) { 
            auto [x1, y1] = path[i - 1];
            auto [x2, y2] = path[i];

            if (x1 == x2) {
                if (y2 > y1) {
                    // Moving up
                    ra.GcellGrid[x1][y1].topEdgeNet++;
                    ra.GcellGrid[x2][y2].bottomEdgeNet++;
                } else {
                    // Moving down
                    ra.GcellGrid[x1][y1].bottomEdgeNet++;
                    ra.GcellGrid[x2][y2].topEdgeNet++;
                }
            } else if (y1 == y2) {
                if (x2 > x1) {
                    // Moving right
                    ra.GcellGrid[x1][y1].rightEdgeNet++;
                    ra.GcellGrid[x2][y2].leftEdgeNet++;
                } else {
                    // Moving left
                    ra.GcellGrid[x1][y1].leftEdgeNet++;
                    ra.GcellGrid[x2][y2].rightEdgeNet++;
                }
            }
        }
    }

    return solutions;
}

string deduceMovementDir(int x1, int y1, int x2, int y2) {
    if (x2 > x1) return "right";
    if (x2 < x1) return "left";
    if (y2 > y1) return "up";
    return "down";
}

void printSolutionFile(vector<vector<pair<int,int>>> solutions, string filename){
    ofstream file(filename);
    int bumpIndex = 1;
    string move, nextMove = "start";
    for (auto solution : solutions){
        file << "n" << bumpIndex << "\n";
        auto [startX, startY] = solution[0];
        auto [nextX, nextY] = solution[1];
        nextMove = deduceMovementDir(startX, startY, nextX, nextY);
        if (nextMove == "left" || nextMove == "right"){
            file << "via" << "\n";
            file << "M2 "; 
        }else{
            file << "M1 ";
        }
        file << startX * grid.Xlength + ra.offsetX << " " << startY * grid.Yheight + ra.offsetY << " ";
        move = nextMove;

        for (int i = 2; i < solution.size(); i++){
            startX = nextX; startY = nextY;
            nextX = solution[i].first; nextY = solution[i].second;
            nextMove = deduceMovementDir(startX, startY, nextX, nextY);
            if (move == nextMove) continue;
            // This prints the endpoint of a metal
            file << startX * grid.Xlength + ra.offsetX << " " << startY * grid.Yheight + ra.offsetY << " \n";
            if (nextMove == "left" || nextMove == "right"){
                file << "via" << "\n";
                file << "M2 "; 
            }else{
                file << "via" << "\n";
                file << "M1 ";
            }
            // This prints the starting point of a metal
            file << startX * grid.Xlength + ra.offsetX << " " << startY * grid.Yheight + ra.offsetY << " ";
            move = nextMove;
        }

        file << solution.back().first * grid.Xlength + ra.offsetX << " " << solution.back().second * grid.Yheight + ra.offsetY << "\n";
        if (move == "left" || move == "right")
            file << "via" << "\n";
        file << ".end" << "\n";
        bumpIndex++;
    }
    return;
}

int main() {
    parseGMP("testcase0.gmp");
    parseGCL("testcase0.gcl");
    parseCST("testcase0.cst");

    int bumpCount = chips[0].bumps.size();
    vector<int> bumpOrder(bumpCount);
    iota(bumpOrder.begin(), bumpOrder.end(), 0); // [0, 1, ..., n-1]

    auto startTime = chrono::steady_clock::now();
    const auto timeLimit = chrono::minutes(20);
    const int maxTries = 1;

    double bestTotalCost = numeric_limits<double>::infinity();
    vector<vector<pair<int,int>>> bestSolution;

    random_device rd;
    mt19937 rng(rd());  // Good random number generator

    int tries = 0;

    while (tries < maxTries) {
        auto now = chrono::steady_clock::now();
        if (now - startTime >= timeLimit) break;

        //shuffle(bumpOrder.begin(), bumpOrder.end(), rng);

        double currentCost = 0.0;
        vector<vector<pair<int,int>>> solution = performRouting(bumpOrder, currentCost);

        for (int y = 0; y < gridRows; y++){
            for (int x = 0; x < gridCols; x++){
                ra.GcellGrid[x][y].topEdgeNet    = 0;
                ra.GcellGrid[x][y].bottomEdgeNet = 0;
                ra.GcellGrid[x][y].leftEdgeNet   = 0;
                ra.GcellGrid[x][y].rightEdgeNet  = 0;
            }
        }

        if (currentCost < bestTotalCost) {
            bestTotalCost = currentCost;
            bestSolution = solution;
        }

        tries++;
    }

    auto endTime = chrono::steady_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(endTime - startTime);

    cout << "Best total cost: " << bestTotalCost << endl;
    cout << "Total random tries: " << tries << endl;
    cout << "Execution time: " << duration.count() << " seconds" << endl;

    printSolutionFile(bestSolution, "testcase0.lg");
    return 0;
}