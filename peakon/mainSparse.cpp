#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <cstdint>
#include <cassert>

using namespace Eigen;

template<size_t N>
struct basis {
    static_assert(N > 0, "Number of uint16_t values must be greater than 0.");
    static constexpr size_t uint16_per_uint64 = sizeof(uint64_t) / sizeof(uint16_t);

    std::vector<uint64_t> values;

    basis() : values((N + uint16_per_uint64 - 1) / uint16_per_uint64, 0) {}

    basis(std::vector<uint16_t> init) : basis() {
        assert(init.size() <= N);
        auto it = init.begin();
        for (size_t i = 0; i < N && it != init.end(); ++i, ++it) {
            (*this)(i) = *it;
        }
    }

    basis(std::initializer_list<uint16_t> init) : basis() {
        assert(init.size() <= N);
        auto it = init.begin();
        for (size_t i = 0; i < N && it != init.end(); ++i, ++it) {
            (*this)(i) = *it;
        }
    }

    template<int size>
    static basis mask2basis(int i) {
        static_assert(size > 0 && size <= 16, "Size must be between 1 and 16");
        std::vector<uint16_t> values;

        for (size_t j = 0; j < N; ++j) {
            int shiftAmount = j * size;
            uint16_t value = (i >> shiftAmount) & ((1 << size) - 1);
            values.push_back(value);
        }

        return basis(values);
    }

    uint16_t& operator()(size_t index) {
        assert(index < N);
        size_t vectorIndex = index / uint16_per_uint64;
        size_t offset = (index % uint16_per_uint64) * sizeof(uint16_t);
        uint8_t* bytePointer = reinterpret_cast<uint8_t*>(&values[vectorIndex]);
        return *reinterpret_cast<uint16_t*>(bytePointer + offset);
    }

    const uint16_t& operator()(size_t index) const {
        assert(index < N);
        size_t vectorIndex = index / uint16_per_uint64;
        size_t offset = (index % uint16_per_uint64) * sizeof(uint16_t);
        const uint8_t* bytePointer = reinterpret_cast<const uint8_t*>(&values[vectorIndex]);
        return *reinterpret_cast<const uint16_t*>(bytePointer + offset);
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
    const int columnWidth = 6;
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

const int VARS = 4;
using basisV = basis<VARS>;
using basis_hashV = basis_hash<VARS>;

struct expr {
    std::unordered_map<basisV, int, basis_hashV> terms;

    expr() = default;
    expr(std::unordered_map<basisV, int, basis_hashV> t)
        : terms{std::move(t)} {}
};

SparseMatrix<double> createA(std::vector<expr> exprs) {
    std::unordered_map<basisV, int, basis_hashV> basis_enumeration;
    std::vector<basisV> reverse_enumeration;
    for (auto& curr_expr : exprs) {
        for (auto& curr_term : curr_expr.terms) {
            if (basis_enumeration.count(curr_term.first) == 0) {
                basis_enumeration.insert({curr_term.first, basis_enumeration.size()});                
            }
        }
    }
    reverse_enumeration = std::vector<basisV>(basis_enumeration.size());
    for (auto& key_value : basis_enumeration) {
        reverse_enumeration[key_value.second] = key_value.first;
    }
    auto A = SparseMatrix<double>(basis_enumeration.size(), exprs.size());
    for (int i = 0; i < exprs.size(); i++) {
        for (auto& curr_term : exprs[i].terms) {
            int j = basis_enumeration[curr_term.first];
            A.insert(j, i) = curr_term.second;
        }
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
            }
        }
    }
    return res;
}

std::vector<basisV> build_uses;
std::vector<expr> build_exprs(std::vector<expr> bases) {
    std::vector<expr> res;
    for (int i = 1; i < (1 << (bases.size() * 1)); i++) {
        auto uses = basisV::mask2basis<1>(i);
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
    bases.emplace_back(expr({{basisV({1, 2, 0, 1}), 1}, {basisV({1, 2, 1, 0}), -1}}));
    bases.emplace_back(expr({{basisV({2, 1, 0, 1}), -1}, {basisV({2, 1, 1, 0}), 1}}));
    bases.emplace_back(expr({{basisV({0, 2, 2, 0}), 1}, {basisV({0, 2, 0, 2}), 1}, {basisV({0, 2, 1, 1}), -2}}));
    bases.emplace_back(expr({{basisV({2, 0, 2, 0}), 1}, {basisV({2, 0, 0, 2}), 1}, {basisV({2, 0, 1, 1}), -2}}));

    std::vector<expr> exprs = build_exprs(bases);

    auto A = createA(exprs);
    A.makeCompressed();
    std::cout << "A:" << std::endl << A << std::endl;

    // Compute the QR decomposition
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
    solver.compute(A);
    auto R = solver.matrixR().topLeftCorner(solver.rank(), solver.rank());
    std::cout << R << std::endl;

    int rank = solver.rank();
    int m = A.cols();

    for (int i = rank; i < m; ++i) {
        VectorXd nullVec = VectorXd::Zero(m);
        nullVec(i) = 1;

        // Back-solve for the other elements of nullVec
        for (int j = i - 1; j >= 0; --j) {
            double sum = 0;
            for (int k = j + 1; k < i; ++k) {
                sum += R.coeff(j, k) * nullVec(k);
            }
            nullVec(j) = -sum / R.coeff(j, j);
        }

        // Now, nullVec is in the null space of A
        // Store or process nullVec as needed
        std::cout << nullVec << std::endl;
    }

    // Check if the matrix is full rank
    if (solver.rank() == A.cols()) {
        std::cout << "Matrix is full rank, null space is empty." << std::endl;
    } else {
        std::cout << "Basis vectors of null space:" << std::endl;
        // Calculate the null space
        for (int i = solver.rank(); i < A.cols(); ++i) {
            VectorXd null_space_vector = VectorXd::Zero(A.cols());
            null_space_vector(i) = 1;
            std::cout << "Vector " << i+1 << ": " << null_space_vector.transpose() << std::endl;
        }
    }

    // Display the null space matrix
    // std::cout << "Null Space Matrix:" << std::endl << null_space << std::endl;
    // print_constants(null_space, {"m_1", "m_2", "x_1", "x_2"} );

    return 0;
}
