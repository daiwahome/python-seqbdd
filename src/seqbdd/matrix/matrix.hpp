#ifndef SEQBDD_MATRIX_HPP__
#define SEQBDD_MATRIX_HPP__

#include <cstddef>
#include <unordered_map>

namespace seqbdd { namespace matrix {
    class BLOSUM62 {
        static const int SIZE = 23;
        static std::unordered_map<char, size_t> INDEX;
        static const int MATRIX[][SIZE];

    public:
        static int at(char x, char y) {
            return MATRIX[INDEX[x]][INDEX[y]];
        }

        BLOSUM62() = delete;
        ~BLOSUM62() = delete;
        BLOSUM62(const BLOSUM62&) = delete;
        BLOSUM62& operator=(const BLOSUM62&) = delete;
        BLOSUM62(BLOSUM62&&) = delete;
        BLOSUM62& operator=(BLOSUM62&&) = delete;
    };
}};

#endif //SEQBDD_MATRIX_HPP__
