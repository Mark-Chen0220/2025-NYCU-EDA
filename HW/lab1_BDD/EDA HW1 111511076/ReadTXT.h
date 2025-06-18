#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void readFileAndParse(const std::string& filename, std::string& boolExpr, std::vector<std::string>& order) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    // Read the first line as boolExpr
    if (std::getline(file, boolExpr)) {
        boolExpr.pop_back();
        // Read remaining lines into order vector
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty()) {  // Skip empty lines if any
                line.pop_back();
                order.push_back(line);
            }
        }
    } else {
        std::cerr << "Error: File is empty" << std::endl;
    }
    
    file.close();
}