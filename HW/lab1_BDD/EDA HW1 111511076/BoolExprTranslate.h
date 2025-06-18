#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cctype>

// Split string by a delimiter
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delimiter)) {
        tokens.push_back(item);
    }
    return tokens;
}

// Convert boolean expression string to C++ logic expression
std::string convertToCppLogic(const std::string& expr) {
    std::string cleaned;
    for (char c : expr) {
        if (!isspace(c)){
            if(isupper(c)){
                cleaned += '!';
                cleaned += tolower(c);
            }
            else cleaned += c;
        }
    }

    std::vector<std::string> terms = split(cleaned, '+');
    std::vector<std::string> convertedTerms;

    for (const std::string& term : terms) {
        std::string andExpr;
        for (size_t i = 0; i < term.size(); ++i) {
            if(term[i] == '!') {
                andExpr += '!';
                continue;
            }
            else andExpr += term[i];
            if (i != term.size() - 1) {
                andExpr += " & ";
            }
        }
        convertedTerms.push_back("(" + andExpr + ")");
    }

    std::string finalExpr;
    for (size_t i = 0; i < convertedTerms.size(); ++i) {
        finalExpr += convertedTerms[i];
        if (i != convertedTerms.size() - 1) {
            finalExpr += " | ";
        }
    }

    return finalExpr;
}
