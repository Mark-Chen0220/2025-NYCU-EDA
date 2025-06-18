#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <unordered_map>

std::string sortExprByOrder(const std::string& expression, 
                              const std::string& variable_order) {
    // Create priority map for custom ordering (case-insensitive)
    std::unordered_map<char, int> priority;
    for (size_t i = 0; i < variable_order.size(); ++i) {
        priority[tolower(variable_order[i])] = i;
    }

    std::istringstream iss(expression);
    std::ostringstream oss;
    std::string product;
    bool first_term = true;

    while (std::getline(iss, product, '+')) {
        if (!first_term) {
            oss << "+";
        }
        first_term = false;

        // Sort literals in this product only
        std::sort(product.begin(), product.end(),
            [&priority](char a, char b) {
                char a_base = tolower(a);
                char b_base = tolower(b);
                
                // Sort by custom order first
                if (priority.at(a_base) != priority.at(b_base)) {
                    return priority.at(a_base) < priority.at(b_base);
                }
                // For same variable, lowercase before uppercase
                return islower(a) && isupper(b);
            });
        
        oss << product;
    }

    return oss.str();
}