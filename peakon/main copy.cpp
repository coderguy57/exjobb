#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <iomanip>

using namespace Eigen;

const int VARS = 4;

void printAsUint16Segments(uint64_t value) {
    // Extracting each 16-bit segment
    uint16_t segment1 = (value >> 48) & 0xFFFF; // Extracts bits 48-63
    uint16_t segment2 = (value >> 32) & 0xFFFF; // Extracts bits 32-47
    uint16_t segment3 = (value >> 16) & 0xFFFF; // Extracts bits 16-31
    uint16_t segment4 = value & 0xFFFF;         // Extracts bits 0-15

    // Printing each segment in columns
    const int columnWidth = 2;
    std::cout << std::setw(columnWidth) << std::setfill(' ') << segment1
              << std::setw(columnWidth) << segment2
              << std::setw(columnWidth) << segment3
              << std::setw(columnWidth) << segment4 << std::endl;
}


struct basis {
    uint64_t value = 0;

    basis() = default;
    basis(uint16_t v1, uint16_t v2, uint16_t v3, uint16_t v4) {
        value = (uint64_t(v4) << 48) + (uint64_t(v3) << 32) + (uint64_t(v2) << 16) + uint64_t(v1);
    }
    template<int size>
    basis static mask2basis(int i) {
        int mask = (1 << size) - 1;
        return basis(i & mask, (i >> size) & mask, (i >> (2 * size)) & mask, (i >> (3 * size)) & mask);
    }

    bool operator==(const basis& other) const {
        return value == other.value;
    }

    uint16_t& operator()(int index) {
        uint8_t* bytePointer = reinterpret_cast<uint8_t*>(&value);
        return *reinterpret_cast<uint16_t*>(bytePointer + 2 * index);
    }

    basis& operator+=(const basis& rhs) // compound assignment (does not need to be a member,
    {                           // but often is, to modify the private members)
        value += rhs.value;
        return *this; // return the result by reference
    }
 
    // friends defined inside class body are inline and are hidden from non-ADL lookup
    friend basis operator+(basis lhs,        // passing lhs by value helps optimize chained a+b+c
                       const basis& rhs) // otherwise, both parameters may be const references
    {
        lhs += rhs; // reuse compound assignment
        return lhs; // return the result by value (uses move constructor)
    }
};

struct basis_hash {
    std::size_t operator()(const basis& b) const {
        return std::hash<uint64_t>()(b.value);
    }
};

struct expr {
    std::unordered_map<basis, int, basis_hash> terms;

    expr() = default;
    expr(std::unordered_map<basis, int, basis_hash> t)
        : terms{std::move(t)} {}
};

expr mul(expr expr1, expr expr2) {
    expr res{};
    for (auto& key_value1 : expr1.terms) {
        for (auto& key_value2 : expr2.terms) {
            auto new_basis = key_value1.first + key_value2.first;
            auto factor = key_value1.second * key_value2.second;
            auto it = res.terms.find(new_basis);
            if (it == res.terms.end()) {
                res.terms.emplace(new_basis, factor);
            } else {
                it->second += factor;
            }
        }
    }
    return res;
}

MatrixXd createA(std::vector<expr> exprs) {
    std::unordered_map<basis, int, basis_hash> basis_enumeration;
    std::vector<basis> reverse_enumeration;
    for (auto& curr_expr : exprs) {
        for (auto& curr_term : curr_expr.terms) {
            if (basis_enumeration.count(curr_term.first) == 0) {
                basis_enumeration.insert({curr_term.first, basis_enumeration.size()});                
            }
        }
    }
    reverse_enumeration = std::vector<basis>(basis_enumeration.size());
    for (auto& key_value : basis_enumeration) {
        reverse_enumeration[key_value.second] = key_value.first;
    }
    MatrixXd A = MatrixXd::Zero(basis_enumeration.size(), exprs.size());
    for (int i = 0; i < exprs.size(); i++) {
        for (auto& curr_term : exprs[i].terms) {
            int j = basis_enumeration[curr_term.first];
            A(j, i) = curr_term.second;
        }
    }
    return A;
}

expr derive(std::vector<expr>& bases, basis used) {
    expr res{};
    for (int i = 0; i < bases.size(); i++) {
        if (used(i) == 0)
            continue;
        
        for (auto& key_value : bases[i].terms) {
            auto factor = key_value.second * used(i);
            auto new_basis = key_value.first + used;
            new_basis(i) -= 1;
            auto it = res.terms.find(new_basis);
            if (it == res.terms.end()) {
                res.terms.emplace(new_basis, factor);
            } else {
                it->second += factor;
            }
        }
    }
    return res;
}

std::vector<basis> build_uses;
std::vector<expr> build_exprs(std::vector<expr> bases) {
    std::vector<expr> res;
    for (int i = 1; i < (1 << (bases.size() * 2)); i++) {
        auto uses = basis::mask2basis<2>(i);
        build_uses.push_back(uses);
        expr new_expr = derive(bases, uses);
        res.emplace_back(std::move(new_expr));
    }
    return res;
}

void print_constants(MatrixXd null_space, std::vector<std::string> names) {
    std::cout << std::endl;
    for (int j = 0; j < null_space.cols(); j++) {
        bool found = false;
        for (int i = 0; i < null_space.rows(); i++) {
            auto factor = null_space(i, j);
            if (factor == 0)
                continue;
            found = true;
            
            std::cout << (factor > 0 ? "+" : "") << factor;
            auto combination = build_uses[i];
            for (int i = 0; i < VARS; i++) {
                if (combination(i) == 1) {
                    std::cout << names[i];
                } else if (combination(i) > 1) {
                    std::cout << names[i] << "^" << combination(i);
                }
            }
        }
        if (found)
            std::cout << std::endl;
    }
}

int main() {
    // MatrixXd A(2, 2);
    // A << 2, -2,
    //      -2, 2;

    std::vector<expr> bases;
    bases.emplace_back(expr({{basis(1, 2, 0, 1), 1}, {basis(1, 2, 1, 0), -1}}));
    bases.emplace_back(expr({{basis(2, 1, 0, 1), -1}, {basis(2, 1, 1, 0), 1}}));
    bases.emplace_back(expr({{basis(0, 2, 2, 0), 1}, {basis(0, 2, 0, 2), 1}, {basis(0, 2, 1, 1), -2}}));
    bases.emplace_back(expr({{basis(2, 0, 2, 0), 1}, {basis(2, 0, 0, 2), 1}, {basis(2, 0, 1, 1), -2}}));

    std::vector<expr> exprs = build_exprs(bases);

    auto A = createA(exprs);
    std::cout << "A:" << std::endl << A << std::endl;

    // Compute the null space using Eigen's FullPivLU decomposition
    FullPivLU<MatrixXd> lu_decomp(A);
    MatrixXd null_space = lu_decomp.kernel();

    // Display the null space matrix
    std::cout << "Null Space Matrix:" << std::endl << null_space << std::endl;
    print_constants(null_space, {"m_1", "m_2", "x_1", "x_2"} );

    return 0;
}
