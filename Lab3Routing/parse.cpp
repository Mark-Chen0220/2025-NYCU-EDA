#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <queue>
#include <utility>

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
    int x1, y1, x2, y2, width, height;
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
            ra = {values[0], values[1], values[0] + values[2], values[1] + values[3], values[2], values[3]};
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
    for (int i = 0; i < gridRows; i++){
        for (int j = 0; j < gridCols; j++){
            getline(file, line); 
            vector<int>values = parseIntegers(line);
            ra.GcellGrid[i][j] = {values[1], values[1], values[0], values[0]};
            //Very Bad Code don't do this
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
    if (cell.topEdgeNet + 1 - cell.topEdgeCapacity > 0) cost += (cell.topEdgeNet - cell.topEdgeCapacity) * overflowCost;
    if (cell.bottomEdgeNet + 1 - cell.bottomEdgeCapacity > 0) cost += (cell.bottomEdgeNet - cell.bottomEdgeCapacity) * overflowCost;
    if (cell.leftEdgeNet + 1 - cell.leftEdgeCapacity > 0) cost += (cell.leftEdgeNet - cell.leftEdgeCapacity) * overflowCost;
    if (cell.rightEdgeNet + 1 - cell.rightEdgeCapacity > 0) cost += (cell.rightEdgeNet - cell.rightEdgeCapacity) * overflowCost;

    cell.cost += cost;
    return;
}

string deduceMovementDir(int startX, int startY, int endX, int endY) {
    int dx = endX - startX;
    int dy = endY - startY;
    if (dx == 0 && dy == 1) {
        return "up";
    } else if (dx == 0 && dy == -1) {
        return "down";
    } else if (dx == -1 && dy == 0) {
        return "left";
    } else if (dx == 1 && dy == 0) {
        return "right";
    } else {
        return "invalid";  // optional: handle diagonal or unexpected movement
    }
}

vector<vector<pair<int,int>>> performRouting(){
    vector<vector<pair<int,int>>> solutions;
    for (int i = 0; i < chips[0].bumps.size(); i++){
        int chip0OffsetX = chips[0].dx, chip0OffsetY = chips[0].dy;
        int chip1OffsetX = chips[1].dx, chip1OffsetY = chips[1].dy;
        int startX = (chip0OffsetX + chips[0].bumps[i].second.first) / grid.Xlength, startY = (chip0OffsetY + chips[0].bumps[i].second.second) / grid.Yheight; 
        int targetX = (chip1OffsetX + chips[1].bumps[i].second.first) / grid.Xlength, targetY = (chip1OffsetY + chips[1].bumps[i].second.second) / grid.Yheight; 
        //Run Manhattan Dist BFS
        vector<vector<int>> ManhattanDist = markManhattanDistance(gridCols, gridRows, grid.Xlength, grid.Yheight, targetX, targetY);
        vector<vector<bool>> beenThere(gridCols, vector<bool>(gridRows, false));
        // Copy Manhattan Distance 
        for (int i = 0; i < gridCols; i++)
            for (int j = 0; j < gridRows; j++)
                ra.GcellGrid[i][j].cost = alpha * ManhattanDist[i][j];

        // Lambda comparator for min-heap
        auto cmp = [](Gcell* lhs, Gcell* rhs) {
            return lhs->cost > rhs->cost;  // min heap
        };

        // Define the priority queue with the lambda comparator
        priority_queue<Gcell*, vector<Gcell*>, decltype(cmp)> minHeap(cmp);
        
        bool targetFound = false; 
        int dirs[4][2] = {{0, 1}, {0, -1}, {-1, 0}, {1, 0}};  
        int currentX = startX, currentY = startY;
        beenThere[startX][startY] = true;
        string currentMetal = "M1";
        int cameFromX = currentX, cameFromY = currentY;
        vector<pair<int,int>> solution;
        solution.push_back(make_pair(startX, startY));

        while (!targetFound){
            for (auto [dx, dy] : dirs) {
                int nx = currentX + dx, ny = currentY + dy;
                // Check boundaries and if cell is unvisited
                if (nx >= 0 && nx < gridCols && ny >= 0 && ny < gridRows ) {
                    if (beenThere[nx][ny]) continue;
                    if (currentMetal == "M1"){
                        if (dx == 0){//Stay as M1
                            ra.GcellGrid[nx][ny].reachedUsingMetal = 1;
                            costFunction(false, ra.GcellGrid[nx][ny]);
                        } else { //Switch to M2
                            ra.GcellGrid[nx][ny].reachedUsingMetal = 2;
                            costFunction(true, ra.GcellGrid[nx][ny]);
                        }
                    }else{ //currentMetal == "M2"
                        if (dx == 0){//Switch to M1
                            ra.GcellGrid[nx][ny].reachedUsingMetal = 1;
                            costFunction(true, ra.GcellGrid[nx][ny]);
                        } else { //Stay as M2
                            ra.GcellGrid[nx][ny].reachedUsingMetal = 2;
                            costFunction(false, ra.GcellGrid[nx][ny]);
                        }
                    }
                    minHeap.push(&ra.GcellGrid[nx][ny]);
                }
            }
            cameFromX = solution.back().first; cameFromY = solution.back().second;
            currentX = minHeap.top()->coordX; currentY = minHeap.top()->coordY;
            int Xmove = cameFromX-currentX, Ymove = cameFromY-currentY;
            solution.push_back(make_pair(currentX, currentY));
            beenThere[currentX][currentY] = true;
            minHeap.pop();
            if (currentX == targetX && currentY == targetY)
                targetFound = true;
            //Update pass-through net count in Gcells
            //dirs[4][2] = {{0, 1}, {0, -1}, {-1, 0}, {1, 0}};  上下左右
            if (Xmove == dirs[0][0] && Ymove == dirs[0][1]){ //上
                ra.GcellGrid[cameFromX][cameFromY].topEdgeNet++;
                ra.GcellGrid[currentX][currentY].bottomEdgeNet++;
                if (currentMetal != "M1") currentMetal = "M1";
            }else if (Xmove == dirs[1][0] && Ymove == dirs[1][1]){ //下
                ra.GcellGrid[cameFromX][cameFromY].bottomEdgeNet++;
                ra.GcellGrid[currentX][currentY].topEdgeNet++;
                if (currentMetal != "M1") currentMetal = "M1";
            }else if (Xmove == dirs[2][0] && Ymove == dirs[2][1]){ //左
                ra.GcellGrid[cameFromX][cameFromY].leftEdgeNet++;
                ra.GcellGrid[currentX][currentY].rightEdgeNet++;
                if (currentMetal != "M2") currentMetal = "M2";
            }else { //右
                ra.GcellGrid[cameFromX][cameFromY].rightEdgeNet++;
                ra.GcellGrid[currentX][currentY].leftEdgeNet++;
                if (currentMetal != "M2") currentMetal = "M2";
            }
        }
        solutions.push_back(solution);
    }
    return solutions;
}

// M1: Vertical
// M2: Horizontal

void printSolution(vector<vector<pair<int,int>>> solutions){
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

        for (int i = 1; i < solution.size() - 1; i++){
            startX = nextX; startY = nextY;
            nextX = solution[i].first; nextY = solution[i].second;
            nextMove = deduceMovementDir(startX, startY, nextX, nextY);

            if (move == nextMove) continue;
            
            cout << nextX << " " << nextY << " \n";
            if (nextMove == "left" || nextMove == "right"){
                cout << "via" << "\n";
                cout << "M2 "; 
            }else{
                cout << "via" << "\n";
                cout << "M1 ";
            }
            cout << nextX << " " << nextY << " ";
            move = nextMove;
        }
        if (move == "left" || move == "right")
            cout << "via" << "\n";
        cout << ".end" << "\n";
        bumpIndex++;
    }
    return;
}

void debugPrint();

int main() {
    string filename = "testcase0.gmp";
    parseGMP(filename);
    filename = "testcase0.gcl";
    parseGCL(filename);
    filename = "testcase0.cst";
    parseCST(filename);
    vector<vector<pair<int,int>>> solutions = performRouting();
    printSolution(solutions);
    return 0;
}

void debugPrint(){
    cout << "Routing Area: (" << ra.x1 << "," << ra.y1 << ") to (" 
         << ra.x2 << "," << ra.y2 << ")\n";
    cout << "Grid: " << grid.Xlength << "x" << grid.Yheight << "\n";
    for (const Chip& chip : chips) {
        int realX = ra.x1 + chip.dx;
        int realY = ra.y1 + chip.dy;
        cout << "\nChip (dx=" << chip.dx << ", dy=" << chip.dy << "):\n";
        cout << "Real Position: (" << realX << "," << realY << ")\n";
        cout << "Bumps (Real Position):\n";
        for (const auto& bump : chip.bumps) {
            int bx = realX + bump.second.first;
            int by = realY + bump.second.second;
            cout << "  ID " << bump.first << ": (" << bx << "," << by << ")\n";
        }
    }
}