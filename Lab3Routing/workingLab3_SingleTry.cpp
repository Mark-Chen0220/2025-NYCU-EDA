#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <queue>
#include <utility>
#include <limits>
#include <algorithm>

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

vector<vector<int>> markManhattanDistance(int gridCols, int gridRows, int gridSizeX, int gridSizeY, int targetX, int targetY) {
    // Initialize distance grid with -1 (unvisited)
    vector<vector<int>> distance(gridCols, vector<int>(gridRows, -1));
    queue<pair<int, int>> q;

    // Start BFS from the target
    distance[targetX][targetY] = 0;
    q.push({targetX, targetY});

    // Directions: up, down, left, right
    int dirs[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    while (!q.empty()) {
        auto [x, y] = q.front();
        q.pop();

        for (auto [dx, dy] : dirs) {
            int nx = x + dx, ny = y + dy;
            // Check boundaries and if cell is unvisited
            if (nx >= 0 && nx < gridCols && ny >= 0 && ny < gridRows && distance[nx][ny] == -1) {
                if (dy == 0) distance[nx][ny] = distance[x][y] + gridSizeX;
                else distance[nx][ny] = distance[x][y] + gridSizeY; //dx = 0 
                q.push({nx, ny});
            }
        }
    }
    return distance;
}

double alpha, beta, myGamma, delta, viaCost, overflowCost = 0;
int overflowWeight = 100;
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
    ss >> beta;
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

void costFunction(bool layerSwap, Gcell &cell){
    double cost = 0;
    if (cell.reachedUsingMetal == 1) cost += myGamma * cell.M1cost;
    else cost += myGamma * cell.M2cost;
    if (layerSwap) cost += myGamma * ((cell.M1cost + cell.M2cost) / 2) + delta * viaCost;
    if (cell.topEdgeNet + 1 - cell.topEdgeCapacity > 0) cost += (cell.topEdgeNet - cell.topEdgeCapacity) * overflowCost * overflowWeight;
    if (cell.bottomEdgeNet + 1 - cell.bottomEdgeCapacity > 0) cost += (cell.bottomEdgeNet - cell.bottomEdgeCapacity) * overflowCost * overflowWeight;
    if (cell.leftEdgeNet + 1 - cell.leftEdgeCapacity > 0) cost += (cell.leftEdgeNet - cell.leftEdgeCapacity) * overflowCost * overflowWeight;
    if (cell.rightEdgeNet + 1 - cell.rightEdgeCapacity > 0) cost += (cell.rightEdgeNet - cell.rightEdgeCapacity) * overflowCost * overflowWeight;

    cell.cost += cost;
    return;
}

vector<vector<pair<int,int>>> performRouting() {
    vector<vector<pair<int,int>>> solutions;

    // Directions: up, down, left, right
    const int dirs[4][2] = {{0,1}, {0,-1}, {-1,0}, {1,0}};

    for (size_t i = 0; i < chips[0].bumps.size(); ++i) {
        // Compute start and target cell coordinates
        int startX = (chips[0].dx + chips[0].bumps[i].second.first) / grid.Xlength;
        int startY = (chips[0].dy + chips[0].bumps[i].second.second) / grid.Yheight;
        int targetX = (chips[1].dx + chips[1].bumps[i].second.first) / grid.Xlength;
        int targetY = (chips[1].dy + chips[1].bumps[i].second.second) / grid.Yheight;

        // Precompute Manhattan-based heuristic cost on grid cells
        auto manDist = markManhattanDistance(gridCols, gridRows, gridCols, gridRows, targetX, targetY);
        for (int x = 0; x < gridCols; ++x) {
            for (int y = 0; y < gridRows; ++y) {
                ra.GcellGrid[x][y].cost = alpha * manDist[x][y];
            }
        }

        // A* data structures
        vector<vector<pair<int,int>>> cameFrom(gridCols, vector<pair<int,int>>(gridRows, {-1,-1}));
        vector<vector<bool>> closedSet(gridCols, vector<bool>(gridRows, false));
        const double INF = numeric_limits<double>::infinity();
        vector<vector<double>> gScore(gridCols, vector<double>(gridRows, INF));

        // Min-heap entry: (fScore, pointer to Gcell)
        struct CellNode {
            double f;
            Gcell* cell;
            string metal;
            bool operator>(const CellNode& other) const { return f > other.f; }
        };

        priority_queue<CellNode, vector<CellNode>, greater<CellNode>> openHeap;

        // Initialize start
        gScore[startX][startY] = 0;
        CellNode startingNode = { ra.GcellGrid[startX][startY].cost, &ra.GcellGrid[startX][startY], "M1" };
        openHeap.push(startingNode);
        bool found = false;

        while (!openHeap.empty()) {
            auto current = openHeap.top();
            openHeap.pop();
            int cx = current.cell->coordX;
            int cy = current.cell->coordY;
            string cm = current.metal;

            if (closedSet[cx][cy]) continue;
            closedSet[cx][cy] = true;

            if (cx == targetX && cy == targetY) {
                found = true;
                break;
            }

            // Explore neighbors
            for (auto &d : dirs) {
                int nx = cx + d[0], ny = cy + d[1];
                if (nx < 0 || nx >= gridCols || ny < 0 || ny >= gridRows) continue;
                if (closedSet[nx][ny]) continue;

                // Determine metal change and calculate cost
                bool switchMetal = (d[0] != 0);
                int nextMetal = (cm == "M1" ^ switchMetal) ? 2 : 1;
                ra.GcellGrid[nx][ny].reachedUsingMetal = nextMetal;
                costFunction(switchMetal, ra.GcellGrid[nx][ny]);

                double tentativeG = gScore[cx][cy] + ra.GcellGrid[nx][ny].cost;
                if (tentativeG < gScore[nx][ny]) {
                    cameFrom[nx][ny] = {cx, cy};
                    gScore[nx][ny] = tentativeG;
                    double fScore = tentativeG + ra.GcellGrid[nx][ny].cost;
                    openHeap.push({fScore, &ra.GcellGrid[nx][ny], nextMetal == 1 ? "M1" : "M2"});
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

        solutions.push_back(path);
    }

    return solutions;
}


// M1: Vertical
// M2: Horizontal

string deduceMovementDir(int x1, int y1, int x2, int y2) {
    if (x2 > x1) return "right";
    if (x2 < x1) return "left";
    if (y2 > y1) return "up";
    return "down";
}

void printSolutionDebug(vector<vector<pair<int,int>>> solutions){
    int bumpIndex = 1;
    string move, nextMove = "start";
    for (auto solution : solutions){
        cout << "n" << bumpIndex << "\n";
        auto [startX, startY] = solution[0];
        auto [nextX, nextY] = solution[1];
        nextMove = deduceMovementDir(startX, startY, nextX, nextY);
        if (nextMove == "left" || nextMove == "right"){
            cout << "via" << "\n";
            cout << "M2 "; 
        }else{
            cout << "M1 ";
        }
        cout << startX << " " << startY << " ";
        move = nextMove;

        for (int i = 2; i < solution.size(); i++){
            startX = nextX; startY = nextY;
            nextX = solution[i].first; nextY = solution[i].second;
            nextMove = deduceMovementDir(startX, startY, nextX, nextY);
            if (move == nextMove) continue;
            // This prints the endpoint of a metal
            cout << startX << " " << startY << " \n";
            if (nextMove == "left" || nextMove == "right"){
                cout << "via" << "\n";
                cout << "M2 "; 
            }else{
                cout << "via" << "\n";
                cout << "M1 ";
            }
            // This prints the starting point of a metal
            cout << startX << " " << startY << " ";
            move = nextMove;
        }

        cout << solution.back().first << " " << solution.back().second << "\n";
        if (move == "left" || move == "right")
            cout << "via" << "\n";
        cout << ".end" << "\n";
        bumpIndex++;
    }
    return;
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
    string filename = "testcase2.gmp";
    parseGMP(filename);
    filename = "testcase2.gcl";
    parseGCL(filename);
    filename = "testcase2.cst";
    parseCST(filename);
    vector<vector<pair<int,int>>> solutions = performRouting();
    filename = "testcase2.lg";
    printSolutionFile(solutions, filename);
    return 0;
}
