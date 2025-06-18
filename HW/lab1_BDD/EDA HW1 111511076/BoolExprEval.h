#include <string>
#include <memory>
#include <stdexcept>
#include <cctype>
#include <vector>
#include <iostream>

class BooleanExpression {
private:

    std::string varOrder;
    int numofVars = 0; 

    enum class TokenType {
        VARIABLE, AND, OR, NOT,
        LPAREN, RPAREN, END
    };

    struct Token {
        TokenType type;
        char value; // Changed to char since variables are single letters
    };

    class ExprNode {
    public:
        virtual ~ExprNode() = default;
        virtual bool evaluate(const std::vector<bool>& values) const = 0;
    };

    class VariableNode : public ExprNode {
        char name;
        const std::string& varOrder; // Reference to the variable order
    public:
        VariableNode(char n, const std::string& order) : name(n), varOrder(order) {
            if (name < 'a' || name > 'z') {
                throw std::runtime_error("Variable must be lowercase a-z");
            }
        }
        bool evaluate(const std::vector<bool>& values) const override {
            size_t index = varOrder.find(name);
            if (index == std::string::npos) {
                throw std::runtime_error("Variable not found in order: " + std::string(1, name));
            }
            if (index >= values.size()) {
                throw std::runtime_error("Not enough values for variable " + std::string(1, name));
            }
            return values[index];
        }
    };

    class NotNode : public ExprNode {
        std::unique_ptr<ExprNode> child;
    public:
        NotNode(std::unique_ptr<ExprNode> c) : child(std::move(c)) {}
        bool evaluate(const std::vector<bool>& values) const override {
            return !child->evaluate(values);
        }
    };

    class BinaryNode : public ExprNode {
        std::unique_ptr<ExprNode> left, right;
        bool isAnd;
    public:
        BinaryNode(std::unique_ptr<ExprNode> l, std::unique_ptr<ExprNode> r, bool is_and)
            : left(std::move(l)), right(std::move(r)), isAnd(is_and) {}
        bool evaluate(const std::vector<bool>& values) const override {
            bool lval = left->evaluate(values);
            bool rval = right->evaluate(values);
            return isAnd ? (lval && rval) : (lval || rval);
        }
    };

    std::string expression;
    size_t pos;
    std::unique_ptr<ExprNode> root;

    Token nextToken() {
        while (pos < expression.length() && std::isspace(expression[pos])) pos++;
        if (pos >= expression.length()) return {TokenType::END, '\0'};

        char c = expression[pos++];
        switch (c) {
            case '&': return {TokenType::AND, '&'};
            case '|': return {TokenType::OR, '|'};
            case '!': return {TokenType::NOT, '!'};
            case '(': return {TokenType::LPAREN, '('};
            case ')': return {TokenType::RPAREN, ')'};
            default:
                if (c >= 'a' && c <= 'z') {
                    return {TokenType::VARIABLE, c};
                }
                throw std::runtime_error("Invalid character: " + std::string(1, c));
        }
    }

    std::unique_ptr<ExprNode> parseExpression() {
        auto left = parseTerm();
        Token token = nextToken();
        while (token.type == TokenType::OR) {
            auto right = parseTerm();
            left = std::make_unique<BinaryNode>(std::move(left), std::move(right), false);
            token = nextToken();
        }
        if (token.type != TokenType::END) {
            pos--; // Backtrack
        }
        return left;
    }

    std::unique_ptr<ExprNode> parseTerm() {
        auto left = parseFactor();
        Token token = nextToken();
        while (token.type == TokenType::AND) {
            auto right = parseFactor();
            left = std::make_unique<BinaryNode>(std::move(left), std::move(right), true);
            token = nextToken();
        }
        if (token.type != TokenType::END) {
            pos--; // Backtrack
        }
        return left;
    }

    std::unique_ptr<ExprNode> parseFactor() {
        Token token = nextToken();
        if (token.type == TokenType::NOT) {
            auto child = parseFactor();
            return std::make_unique<NotNode>(std::move(child));
        }
        else if (token.type == TokenType::LPAREN) {
            auto expr = parseExpression();
            token = nextToken();
            if (token.type != TokenType::RPAREN) {
                throw std::runtime_error("Missing closing parenthesis");
            }
            return expr;
        }
        else if (token.type == TokenType::VARIABLE) {
            return std::make_unique<VariableNode>(token.value, varOrder);
        }
        else {
            throw std::runtime_error("Unexpected token at position " + std::to_string(pos));
        }
    }

public:

    BooleanExpression(const std::string& expr) : expression(expr), pos(0) {
        root = parseExpression();
        Token last = nextToken();
        if (last.type != TokenType::END) {
            throw std::runtime_error("Extra tokens after expression");
        }
    }

    BooleanExpression(const std::string& expr, const std::string& order) 
        : expression(expr), pos(0), varOrder(order) {
        root = parseExpression();
        Token last = nextToken();
        if (last.type != TokenType::END) {
            throw std::runtime_error("Extra tokens after expression");
        }
    }

    bool evaluate(const std::vector<bool>& values) const {
        return root->evaluate(values);
    }
};
