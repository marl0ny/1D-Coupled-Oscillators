/* 
This source file is a re-implementation of this script:

https://gist.github.com/marl0ny/db6aa96035046ba597d76dcf52f6921b

Coupled Fermions by constructing the Fermion
operators as matrices and then diagonalizing the subsequent Hamiltonian
to find the energies and energy eigenvectors.

The primary reference for this is Chapter 4 of 
Introduction to Many Body Physics by Piers Coleman, pg 71 - 77 on the
Jordan-Wigner transformation.

    Coleman P., "Simple examples of second quantization",
    <i>Introduction to Many Body Physics</i>,
    Cambridge University Press, 2015,  pg 71-94.

*/
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <iostream>
#include <vector>

using Eigen::ComplexEigenSolver;
typedef Eigen::Matrix2cd Matrix2x2;
typedef Eigen::VectorXcd Vector;
typedef Eigen::MatrixXcd Matrix;
typedef std::vector<Matrix> MatrixList;

static const double PI = 3.141592653589793;

/* Fermionic creation operator for a single particle. */
Matrix2x2 get_ladder_up() {
    Matrix2x2 ladder_up_mat;
    ladder_up_mat << 
        0.0, 1.0,
        0.0, 0.0;
    return ladder_up_mat;
}

/* The "strand" corresponds to the operator exp(-n i pi),
where n counts the occupation. */
Matrix2x2 get_strand() {
    Matrix2x2 strand_mat;
    strand_mat << 
        -1.0, 0.0,
        0.0, 1.0;
    return strand_mat;
}

static size_t powl(size_t base, size_t exponent) {
    size_t res = 1;
    for (int i = 1; i <= exponent; i++)
        res *= base;
    return res;
}

/* Get the state that corresponds to zero occupancy.*/
Vector get_vacant_state(size_t n_oscillators) {
    size_t n = powl(2, n_oscillators);
    Vector ground {n};
    for (int i = 0; i < n; i++)
        ground[i] = (i == n-1)? 1.0: 0.0;
    return ground;
}

/* Get the Kronecker product of two matrices.*/
Matrix kron(const Matrix &a, const Matrix &b) {
    Matrix res(a.rows()*b.rows(), a.cols()*b.cols());
    for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.cols(); j++) {
            res.block(i*b.rows(), j*b.cols(), b.rows(), b.cols())
                = a(i, j)*b;
        }
    }
    return res;
}

Matrix kron(const Matrix &a, const Matrix2x2 &b) {
    Matrix res(a.rows()*b.rows(), a.cols()*b.cols());
    for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.cols(); j++) {
            res.block(i*b.rows(), j*b.cols(), b.rows(), b.cols())
                = a(i, j)*b;
        }
    }
    return res;
}

Matrix kron(const Matrix2x2 &a, const Matrix &b) {
    Matrix res(a.rows()*b.rows(), a.cols()*b.cols());
    for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.cols(); j++) {
            res.block(i*b.rows(), j*b.cols(), b.rows(), b.cols())
                = a(i, j)*b;
        }
    }
    return res;
}

MatrixList fourier_transform(const MatrixList &x) {
    MatrixList f;
    size_t n = x.size();
    for (int k = 0; k < n; k++) {
        f.push_back(x[0]);
        for (int j = 1; j < n; j++) {
            f[k] += x[j]*std::exp(std::complex<double>(0.0, -2.0*PI*j*k/n));
        }
        f[k] *= 1.0/std::sqrt(n);
    }
    return f;
}

// Matrix kron_power(
//     Matrix base, size_t exponent) {
//     Matrix res = (exponent == 0)? 
//         Matrix::Identity(base.rows(), base.cols()): base;
//     for (int i = 0; i < (exponent - 1); i++)
//         res = kron(res, base);
//     return res;
// }

/*
Returns a list of the coupled Fermion creation operators for n particles.

The way this works is that one first starts with ladder spin operators
for a one-dimensional linear sequence of spins, where one now interprets
spin down as corresponding to the state of zero occupancy, and spin up as
the occupied state. The main problem with this approach is that the ladder
spin operators do not follow Fermion anticommuting relations, rather they 
commute. To fix this, use the Jordan-Wigner transformation, where each
ladder spin operator in the sequence is matrix multiplied by an operator
exp(-i n pi), where n counts the occupation of states for **only** the
preceding ladder spin operators in the sequence See the following:

    Coleman P., "Simple examples of second quantization",
    <i>Introduction to Many Body Physics</i>,
    Cambridge University Press, 2015,  pg 71-94.

*/
MatrixList get_creation_operators(size_t n) {
    MatrixList ladders_u;
    Matrix string (1, 1);
    string(0, 0) = 1;
    Matrix2x2 ladder_up = get_ladder_up();
    Matrix2x2 strand = get_strand();
    for (int i = 0; i < n; i++) {
        size_t id_size = powl(2, n-i-1);
        Matrix identity = Matrix::Identity(id_size, id_size);
        Matrix o = kron(kron(string, ladder_up),identity);
        ladders_u.push_back(o);
        string = kron(strand, string);
    }
    return ladders_u;
}

MatrixList get_annihilators_from_creators(const MatrixList &creators) {
    MatrixList annihilators;
    for (auto &creator: creators) {
        annihilators.push_back(creator.conjugate().transpose());
    }
    return annihilators;
}

Matrix make_density_vector_matrix(size_t n_oscillators) {
    size_t n_states = powl(2, n_oscillators);
    Matrix n_mat (n_oscillators, n_states);
    int val = 1;
    for (int i = 0; i < n_oscillators; i++) {
        for (int j = 0; j < n_states; j++) {
            double element = double(int((val&(~j)) > 0));
            n_mat(n_oscillators-i-1, j) = element;
        }
        val *= 2.0;
    }
    return n_mat;
}

Matrix get_hamiltonian(
    const Matrix &couplings, const MatrixList &annihilators) {
    size_t n_states = annihilators[0].rows();
    Matrix h = Matrix::Zero(n_states, n_states);
    for (int i = 0; i < couplings.rows(); i++)
        for (int j = 0; j < couplings.cols(); j++)
            if (couplings(i, j) != 0.0)
                h += 
                    annihilators[i].adjoint()*couplings(i, j)*annihilators[j];
    return h;
}

int main() {
    /* Compile with:
    clang++ --std=c++17 -lm -O3 -ffast-math -I./ main.cpp -o program; ./program
    */
    int n_oscillators = 8;
    Matrix couplings = Matrix::Zero(n_oscillators, n_oscillators);
    for (int i = 0; i < n_oscillators; i++) {
        if (i + 1 < n_oscillators)
            couplings(i, i+1) = -1.0;
        couplings(i, i) = 2.0;
        if (i - 1 >= 0)
            couplings(i, i-1) = -1.0;
    }
    couplings(0, n_oscillators-1) = -1.0;
    couplings(n_oscillators-1, 0) = -1.0;
    // std::cout << couplings << std::endl;

    Vector ground = get_vacant_state(n_oscillators);
    Matrix density_mat = make_density_vector_matrix(n_oscillators);
    MatrixList creators = get_creation_operators(n_oscillators);
    MatrixList annihilators = get_annihilators_from_creators(creators);
    MatrixList annihilators_k = fourier_transform(annihilators);
    Matrix h = get_hamiltonian(couplings, annihilators);
    ComplexEigenSolver<Matrix> solver(h);
    // std::cout << solver.eigenvalues() << std::endl;
    std::cout << 
        density_mat*(solver.eigenvectors().col(1)).cwiseAbs2() << std::endl;

    // for (auto &up: ups)
    //     std::cout << up << std::endl;
    // for (int i = 0; i < n_oscillators; i++) {
    //     Vector a = density_mat*creators_k[i]*ground;
    //     std::cout << a << std::endl << std::endl;
    // }
    return 0;
}