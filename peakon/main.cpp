#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <cstdint>
#include <cassert>
#include <set>
#include <numeric>

using namespace Eigen;

template<size_t N>
struct basis {
    static_assert(N > 0, "Number of uint32_t values must be greater than 0.");
    static constexpr size_t uint32_per_uint64 = sizeof(uint64_t) / sizeof(uint32_t);
    static constexpr size_t size = N;

    std::vector<uint64_t> values;

    basis() : values((N + uint32_per_uint64 - 1) / uint32_per_uint64, 0) {}

    basis(std::vector<uint32_t> init) : basis() {
        assert(init.size() <= N);
        auto it = init.begin();
        for (size_t i = 0; i < N && it != init.end(); ++i, ++it) {
            (*this)(i) = *it;
        }
    }

    basis(std::initializer_list<uint32_t> init) : basis() {
        assert(init.size() <= N);
        auto it = init.begin();
        for (size_t i = 0; i < N && it != init.end(); ++i, ++it) {
            (*this)(i) = *it;
        }
    }

    template<int size>
    static basis mask2basis(int i) {
        static_assert(size > 0 && size <= 16, "Size must be between 1 and 16");
        std::vector<uint32_t> values;

        for (size_t j = 0; j < N; ++j) {
            int shiftAmount = j * size;
            uint32_t value = (i >> shiftAmount) & ((1 << size) - 1);
            values.push_back(value);
        }

        return basis(values);
    }

    uint32_t& operator()(size_t index) {
        assert(index < N);
        size_t vectorIndex = index / uint32_per_uint64;
        size_t offset = (index % uint32_per_uint64) * sizeof(uint32_t);
        uint8_t* bytePointer = reinterpret_cast<uint8_t*>(&values[vectorIndex]);
        return *reinterpret_cast<uint32_t*>(bytePointer + offset);
    }

    const uint32_t& operator()(size_t index) const {
        assert(index < N);
        size_t vectorIndex = index / uint32_per_uint64;
        size_t offset = (index % uint32_per_uint64) * sizeof(uint32_t);
        const uint8_t* bytePointer = reinterpret_cast<const uint8_t*>(&values[vectorIndex]);
        return *reinterpret_cast<const uint32_t*>(bytePointer + offset);
    }

    basis& operator+=(const basis& rhs) {
        assert(values.size() == rhs.values.size());
        for (size_t i = 0; i < values.size(); ++i) {
            values[i] += rhs.values[i];
        }
        return *this;
    }

    friend basis operator+(basis lhs, const basis& rhs) {
        lhs += rhs;
        return lhs;
    }

    bool operator==(const basis& other) const {
        return values == other.values;
    }
};

template<size_t N>
void print_basis(basis<N> b) {
    const int columnWidth = 3;
    for (size_t i = 0; i < N; ++i) {
        std::cout << std::setw(columnWidth) << std::setfill(' ') << b(i);
    }
    std::cout << std::endl;
}

template<size_t N>
struct basis_hash {
    std::size_t operator()(const basis<N>& b) const {
        std::size_t h = 0;
        for (const auto& val : b.values) {
            h ^= std::hash<uint64_t>()(val) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

const int VARS = 8;
using basisV = basis<VARS>;
using basis_hashV = basis_hash<VARS>;

struct expr {
    std::unordered_map<basisV, int, basis_hashV> terms;

    expr() = default;
    expr(std::unordered_map<basisV, int, basis_hashV> t)
        : terms{std::move(t)} {}
};

void print_expr(expr expression) {
    const int columnWidth = 3;
    std::cout << std::endl;
    for (auto& key_value : expression.terms) {
        std::cout << std::setw(columnWidth) << std::setfill(' ') << key_value.second;
        for (size_t i = 0; i < key_value.first.size; ++i) {
            std::cout << std::setw(columnWidth) << std::setfill(' ') << key_value.first(i);
        }
        std::cout << std::endl;
    }
}



void clear_unused(
    std::unordered_map<basisV, int, basis_hashV>& basis_count,
    std::unordered_map<basisV, int, basis_hashV>& basis_enumeration,
    std::vector<std::set<int>>& reverse_count,
    std::vector<expr>& exprs,
    int i)
{
    for (auto& curr_term : exprs[i].terms) {
        if (basis_count[curr_term.first] == 0) {
            std::cout << "Error: basis_count is 0" << std::endl;
        }
        basis_count[curr_term.first]--;
        if (reverse_count[basis_enumeration[curr_term.first]].count(i) == 0) {
            std::cout << "Error: reverse_count is 0" << std::endl;
        }
        reverse_count[basis_enumeration[curr_term.first]].erase(i);
    }
    for (auto& curr_term : exprs[i].terms) {
        if (basis_count[curr_term.first] == 1) {
            int id = *reverse_count[basis_enumeration[curr_term.first]].begin();
            clear_unused(basis_count, basis_enumeration, reverse_count, exprs, id);
        }
    }
}

std::unordered_map<int, int> indices_map;
MatrixXd createA(std::vector<expr> exprs) {
    std::unordered_map<basisV, int, basis_hashV> basis_enumeration;
    for (auto& curr_expr : exprs) {
        for (auto& curr_term : curr_expr.terms) {
            if (basis_enumeration.count(curr_term.first) == 0) {
                basis_enumeration.insert({curr_term.first, basis_enumeration.size()});                
            }
        }
    }
    std::unordered_map<basisV, int, basis_hashV> basis_count;
    std::vector<std::set<int>> reverse_count(basis_enumeration.size());
    for (int i = 0; i < exprs.size(); i++) {
        for (auto& curr_term : exprs[i].terms) {
            int j = basis_enumeration[curr_term.first];
            if(basis_count.contains(curr_term.first)) {
                basis_count[curr_term.first]++;
            } else {
                basis_count[curr_term.first] = 1;
            }
            reverse_count[j].insert(i);
        }
    }
    for (int i = 0; i < exprs.size(); i++) {
        bool is_lonely = false;
        for (auto& curr_term : exprs[i].terms) {
            if (basis_count[curr_term.first] == 1) {
                is_lonely = true;
                break;
            }
        }
        if (is_lonely)
            clear_unused(basis_count, basis_enumeration, reverse_count, exprs, i);
    }
    int lonely_basis = 0;
    for (auto& key_value : basis_count) {
        if (key_value.second == 0)
            lonely_basis++;
    }
    int lonely_expr = 0;
    for (int i = 0; i < exprs.size(); i++) {
        for (auto& curr_term : exprs[i].terms) {
            if (basis_count[curr_term.first] == 0) {
                lonely_expr++;
                break;
            }
        }
    }
    std::cout << "rows: " << basis_count.size() - lonely_basis << " cols: " << exprs.size() - lonely_expr << std::endl; 
    MatrixXd A = MatrixXd::Zero(basis_count.size() - lonely_basis, exprs.size() - lonely_expr);
    std::unordered_map<basisV, int, basis_hashV> basis_enumeration2;
    int expr_id = 0;
    for (int i = 0; i < exprs.size(); i++) {
        bool used = true;
        for (auto& curr_term : exprs[i].terms) {
            if(basis_count[curr_term.first] == 0) {
                used = false;
                break;
            }
        }
        if (!used)
            continue;
        for (auto& curr_term : exprs[i].terms) {
            if (basis_enumeration2.count(curr_term.first) == 0) {
                basis_enumeration2.insert({curr_term.first, basis_enumeration2.size()});                
            }
            int j = basis_enumeration2[curr_term.first];
            A(j, expr_id) = curr_term.second;
        }
        indices_map[expr_id] = i;
        expr_id++;
    }
    return A;
}

expr derive(std::vector<expr>& bases, basisV used) {
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
                if (it->second == 0)
                    res.terms.erase(it);
            }
        }
    }
    return res;
}

void generate_powers(int max_power, int from, int to, basisV& basis, std::function<void()> callback) {
    int size = max_power;
    std::vector<int> powers(size, 0);
    int depth = size;
    while (depth >= 0) {
        if (depth == size) {
            for (int i = from; i < to; i++)
                basis(i) = 0;
            for (int i = 0; i < size; i++)
                basis(powers[i] + from) += 1;

            callback();
            depth--;
            continue;
        }
        if (powers[depth] > (to - from) - 1) {
            powers[depth] = powers[depth - 1];
            depth++;
        } else if(powers[depth] < (to - from) - 1) {
            powers[depth++]++;
        } else {
            powers[depth]++;
            depth--;
        }
    }
}

std::vector<expr> build_exprs(std::vector<expr>& bases, std::vector<basisV>& build_uses, int m_power, int x_power) {
    std::vector<expr> res;

    auto basis = basisV{};
    auto K = VARS / 2;
    generate_powers(m_power, 0, K, basis, [&]() {
        generate_powers(x_power, K, 2*K, basis, [&]() {
            build_uses.push_back(basis);
            expr new_expr = derive(bases, basis);
            res.emplace_back(std::move(new_expr));
        });
    });
    std::cout << res.size() << std::endl;
    return res;
}

void print_constants(std::vector<expr> exprs, std::vector<std::string> names) {
    std::cout << std::endl;
    for (int j = 0; j < exprs.size(); j++) {
        auto expr = exprs[j];
        for (auto& key_value : expr.terms) {
            auto factor = key_value.second;
            if (std::abs(std::abs(factor) - 1) < 1e-10)
                std::cout << (factor > 0 ? "+" : "-");
            else
                std::cout << (factor > 0 ? "+" : "") << factor;

            auto combination = key_value.first;
            for (int k = 0; k < VARS; k++) {
                if (combination(k) == 1) {
                    std::cout << names[k];
                } else if (combination(k) > 1) {
                    std::cout << names[k] << "^" << combination(k);
                }
            }
        }
        std::cout << std::endl << std::endl;
    }
}

void print_python(std::vector<expr> exprs, std::vector<std::string> names) {
    std::cout << std::endl;
    for (int j = 0; j < exprs.size(); j++) {
        auto expr = exprs[j];
        for (auto& key_value : expr.terms) {
            auto factor = key_value.second;
            std::cout << (factor > 0 ? "+" : "") << factor;

            auto combination = key_value.first;
            for (int k = 0; k < VARS; k++) {
                if (combination(k) == 0)
                    continue;

                std::cout << "*";
                
                if (combination(k) == 1) {
                    std::cout << names[k];
                } else if (combination(k) > 1) {
                    std::cout << names[k] << "**" << combination(k);
                }
            }
        }
        std::cout << std::endl << std::endl;
    }
    std::cout << std::endl;
}

expr add(expr const& expr1, expr const& expr2) {
    expr res{};
    for (auto& key_value : expr1.terms) {
        auto it = expr2.terms.find(key_value.first);
        if (it == expr2.terms.end()) {
            res.terms.emplace(key_value.first, key_value.second);
        } else {
            int factor = key_value.second + it->second;
            if (factor != 0)
                res.terms.emplace(key_value.first, factor);
        }
    }
    for (auto& key_value : expr2.terms) {
        auto it = expr1.terms.find(key_value.first);
        if (it == expr1.terms.end()) {
            res.terms.emplace(key_value.first, key_value.second);
        }
    }
    return res;
}

expr mul(expr const& expr1, int factor) {
    expr res = expr1;
    for (auto& key_value1 : res.terms) {
        key_value1.second *= factor;
    }
    return res;
}

expr mul(expr const& expr1, expr const& expr2) {
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

bool equal(expr const& expr1, expr const& expr2) {
    if (expr1.terms.size() != expr2.terms.size())
        return false;
    if (expr1.terms.size() == 0)
        return true;
    
    auto it1 = expr1.terms.begin();
    auto it2 = expr2.terms.find(it1->first);
    if (it2 == expr2.terms.end())
        return false;
    
    auto factor = double(it1->second) / it2->second;
    for (auto& key_value : expr1.terms) {
        auto it = expr2.terms.find(key_value.first);
        if (it == expr2.terms.end())
            return false;
         
        if (std::abs(double(key_value.second) / it->second - factor) > 1e-10)
            return false;
    }
    return true;
}

std::vector<expr> remove_duplicates(MatrixXd& null_space, std::vector<basisV>& build_uses) {
    std::cout << std::endl;
    std::vector<expr> exprs;
    for (int j = 0; j < null_space.cols(); j++) {
        bool found = false;
        int lcm = 1;
        for (int i = 0; i < null_space.rows(); i++) {
            auto factor = null_space(i, j);
            if (factor == 0)
                continue;
            found = true;

            int int_factor = std::round(std::abs(1 / factor));
            lcm = std::lcm(lcm, int_factor);
        }
        expr res{};
        for (int i = 0; i < null_space.rows(); i++) {
            auto factor = null_space(i, j);
            if (factor == 0)
                continue;
            
            int int_factor = std::round(lcm * factor);

            auto combination = build_uses[indices_map[i]];
            res.terms.insert({combination, int_factor});
        }
        exprs.push_back(res);
    }
    return exprs;
    // int n = exprs.size();
    // for (int i = 0; i < n; i++) {
    //     for (int j = i; j < n; j++) {
    //         exprs.push_back(mul(exprs[i], exprs[j]));
    //     }
    // }
    // std::unordered_map<basisV, int, basis_hashV> basis_enumeration;
    // for (auto& curr_expr : exprs) {
    //     for (auto& curr_term : curr_expr.terms) {
    //         if (basis_enumeration.count(curr_term.first) == 0) {
    //             basis_enumeration.insert({curr_term.first, basis_enumeration.size()});                
    //         }
    //     }
    // }
    // MatrixXd A = MatrixXd::Zero(basis_enumeration.size(), exprs.size());
    // for (int i = 0; i < exprs.size(); i++) {
    //     for (auto& curr_term : exprs[i].terms) {
    //         int j = basis_enumeration[curr_term.first];
    //         A(j, i) = curr_term.second;
    //     }
    // }
    // FullPivLU<MatrixXd> lu_decomp(A);
    // auto& permMatrix = lu_decomp.permutationQ();
    // auto& kernel = lu_decomp.kernel();
    // std::cout << "Exprs size: " << exprs.size() << std::endl;
    // std::cout << "Rank: " << lu_decomp.rank() << std::endl;
    // std::cout << "Null space size: " << kernel.cols() << std::endl;

    // std::vector<int> pivotalIndices(exprs.size(), 0);
    // for (int i = 0; i < lu_decomp.rank(); ++i) {
    //     // Apply the permutation inverse to get original column indices
    //     int originalIndex = permMatrix.indices()(i);
    //     pivotalIndices[originalIndex] = 1;
    // }

    // std::vector<expr> valid_exprs;
    // for (int i = 0; i < exprs.size(); i++) {
    //     auto expr1 = exprs[i];
    //     bool found = false;
    //     for (int j = 0; j < exprs.size(); j++) {
    //         if (found)
    //             break;
    //         for (int k = j; k < exprs.size(); k++) {
    //             auto expr2 = mul(exprs[j], exprs[k]);
    //             if (equal(expr1, expr2)) {
    //                 found = true;
    //                 break;
    //             }
    //         }
    //     }
    //     if (!found && pivotalIndices[i] == 1) {
    //         valid_exprs.push_back(expr1);
    //     }
    // }
    // return valid_exprs;
}


std::vector<expr> build_bases3() {
    std::vector<expr> u_k;
    std::vector<expr> ux_k;
    constexpr int K = VARS / 2;
    for (int k = 0; k < K; k++) {
        expr u;
        expr ux;
        for (int i = 0; i < K; i++) {
            if (i == k)
                continue;

            int factor = k > i ? 1 : -1;

            auto coeffx = std::vector<uint32_t>(VARS, 0);
            coeffx[i] = 1;
            ux.terms.insert({basisV(coeffx), factor});

            int var_id = i + k - 1;
            auto coeff = coeffx;
            if (var_id == 1) {
                coeff[0 + K] = 1;
                u.terms.insert({basisV(coeff), 1});
                coeff[0 + K] = 0;
                coeff[1 + K] = 1;
                u.terms.insert({basisV(coeff), 1});
                coeff[1 + K] = 0;
                coeff[2 + K] = 1;
                u.terms.insert({basisV(coeff), 1});
            } else {
                coeff[var_id + K] = 1;
                u.terms.insert({basisV(coeff), 1});
            }

        }
        u_k.emplace_back(std::move(u));
        ux_k.emplace_back(std::move(ux));
    }
    std::vector<expr> bases;
    // m
    for (int k = 0; k < K; k++) {
        auto coeff = std::vector<uint32_t>(VARS, 0);
        coeff[k] = 1;
        auto product = mul(mul(ux_k[k], u_k[k]), expr({{basisV(coeff), -1}}));
        bases.push_back(std::move(product));
    }
    // p
    {
        auto x1 = mul(u_k[0], u_k[0]);
        auto x2 = mul(u_k[1], u_k[1]);
        auto x3 = mul(u_k[2], u_k[2]);

        bases.push_back(add(x2, mul(x1, -1)));
        // bases.push_back(add(x3, mul(x1, -1)));
        bases.push_back({});
        bases.push_back(add(x3, mul(x2, -1)));
    }
    return bases;
}

std::vector<expr> build_bases2() {
    std::vector<expr> u_k;
    std::vector<expr> ux_k;
    constexpr int K = VARS / 2;
    for (int k = 0; k < K; k++) {
        expr u;
        expr ux;
        for (int i = 0; i < K; i++) {
            if (i == k)
                continue;
            auto coeffx = std::vector<uint32_t>(VARS, 0);
            coeffx[i] = 1;

            int var_id = i + k - 1;
            auto coeff = coeffx;
            coeff[var_id + K] = 1;
            int factor = k > i ? 1 : -1;

            ux.terms.insert({basisV(coeffx), factor});
            u.terms.insert({basisV(coeff), 1});
        }
        u_k.emplace_back(std::move(u));
        ux_k.emplace_back(std::move(ux));
    }
    std::vector<expr> bases;
    // m
    for (int k = 0; k < K; k++) {
        auto coeff = std::vector<uint32_t>(VARS, 0);
        coeff[k] = 1;
        auto product = mul(mul(ux_k[k], u_k[k]), expr({{basisV(coeff), -1}}));
        bases.push_back(std::move(product));
    }
    // p
    {
        auto x1 = mul(u_k[0], u_k[0]);
        auto x2 = mul(u_k[1], u_k[1]);
        auto x3 = mul(u_k[2], u_k[2]);

        bases.push_back(add(x2, mul(x1, -1)));
        bases.push_back(add(x3, mul(x1, -1)));
        bases.push_back(add(x3, mul(x2, -1)));
    }
    return bases;
}

std::vector<expr> build_bases() {
    std::vector<expr> u_k;
    std::vector<expr> ux_k;
    constexpr int K = VARS / 2;
    for (int k = 0; k < K; k++) {
        expr u;
        expr ux;
        for (int i = 0; i < K; i++) {
            if (i == k)
                continue;
            auto coeffx = std::vector<uint32_t>(VARS, 0);
            coeffx[i] = 1;
            auto coeff1 = coeffx;
            coeff1[i + K] = 1;
            auto coeff2 = coeffx;
            coeff2[k + K] = 1;
            int factor = k > i ? 1 : -1;

            // auto test = coeffx;
            // test[i+K] = 1;

            ux.terms.insert({basisV(coeffx), factor});
            u.terms.insert({basisV(coeff1), -factor});
            u.terms.insert({basisV(coeff2), factor});
            // u.terms.insert({basisV(test), 1});
        }
        u_k.emplace_back(std::move(u));
        ux_k.emplace_back(std::move(ux));
    }
    std::vector<expr> bases;
    // m
    for (int k = 0; k < K; k++) {
        auto coeff = std::vector<uint32_t>(VARS, 0);
        coeff[k] = 1;
        auto product = mul(mul(ux_k[k], u_k[k]), expr({{basisV(coeff), -1}}));
        bases.push_back(std::move(product));
    }
    // x
    for (int k = 0; k < K; k++) {
        bases.push_back(mul(u_k[k], u_k[k]));
    }
    return bases;
}

int main() {
    // MatrixXd A(2, 2);
    // A << 2, -2,
    //      -2, 2;

    auto bases = build_bases();
    for (auto& base : bases) {
        print_expr(base);
    }
    // std::vector<expr> bases;
    // bases.emplace_back(expr({{basisV({1, 2, 0, 1}), 1}, {basisV({1, 2, 1, 0}), -1}}));
    // bases.emplace_back(expr({{basisV({2, 1, 0, 1}), -1}, {basisV({2, 1, 1, 0}), 1}}));
    // bases.emplace_back(expr({{basisV({0, 2, 2, 0}), 1}, {basisV({0, 2, 0, 2}), 1}, {basisV({0, 2, 1, 1}), -2}}));
    // bases.emplace_back(expr({{basisV({2, 0, 2, 0}), 1}, {basisV({2, 0, 0, 2}), 1}, {basisV({2, 0, 1, 1}), -2}}));

    std::vector<basisV> build_uses;
    int m_power = 6;
    int x_power = 4;
    std::vector<expr> exprs = build_exprs(bases, build_uses, m_power, x_power);
    for (auto& base : exprs) {
        // print_expr(base);
    }

    std::cout << "Exprs size: " << exprs.size() << std::endl;
    auto A = createA(exprs);
    // std::cout << "A:" << std::endl << A << std::endl;
    std::cout << "A - cols :" << A.cols() << std::endl << "A - rows :" << A.cols() << std::endl;

    // Compute the null space using Eigen's FullPivLU decomposition
    FullPivLU<MatrixXd> lu_decomp(A);
    MatrixXd null_space = lu_decomp.kernel();
    null_space = null_space.unaryExpr([](double x){return (abs(x)<1e-12)?0.:x;});

    // Display the null space matrix
    std::cout << "Null Space Matrix:" << std::endl << null_space << std::endl;
    auto constants = remove_duplicates(null_space, build_uses);
    // print_constants(null_space, {"m_1", "m_2", "x_1", "x_2"} );
    // print_constants(constants, {"m_1", "m_2", "m_3", "x_1", "x_2", "x_3"} );
    // print_python(constants, {"m1", "m2", "m3", "x1", "x2", "x3"} );
    // print_constants(constants, {"m_1", "m_2", "m_3", "x_1", "x_2", "x_3"} );
    print_python(constants, {"m1", "m2", "m3", "m4", "x1", "x2", "x3", "x4"} );
    // print_python(constants, {"m1", "m2", "x1", "x2"} );
    // print_constants(null_space, {"m_1", "m_2", "m_3", "\\varphi_1", "M", "\\varphi_3"} );

    return 0;
}
