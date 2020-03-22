#include <vector>

template <typename Scalar>
class Matrix
{
private:
    std::vector<std::vector<Scalar>> mat;
    size_t rows_;
    size_t cols_;

public:
    Matrix(size_t _rows, size_t _cols, const Scalar &initial = Scalar(0));
    Matrix(const Matrix<Scalar> &rhs);
    ~Matrix();
    Matrix<Scalar> transpose();

    //Matrix Matrix operations
    Matrix<Scalar> &operator=(const Matrix<Scalar> &rhs);
    Matrix<Scalar> operator+(const Matrix<Scalar> &rhs);
    Matrix<Scalar> &operator+=(const Matrix<Scalar> &rhs);
    Matrix<Scalar> operator-(const Matrix<Scalar> &rhs);
    Matrix<Scalar> &operator-=(const Matrix<Scalar> &rhs);
    Matrix<Scalar> operator*(const Matrix<Scalar> &rhs);
    Matrix<Scalar> &operator*=(const Matrix<Scalar> &rhs);

    //Matrix Scalar operations
    Matrix<Scalar> operator+(const Scalar &rhs);
    Matrix<Scalar> &operator+=(const Scalar &rhs);
    Matrix<Scalar> operator-(const Scalar &rhs);
    Matrix<Scalar> &operator-=(const Scalar &rhs);
    Matrix<Scalar> operator*(const Scalar &rhs);
    Matrix<Scalar> &operator*=(const Scalar &rhs);
    Matrix<Scalar> operator/(const Scalar &rhs);
    Matrix<Scalar> &operator/=(const Scalar &rhs);

    Scalar &operator()(const size_t &row, const size_t &col);
    const Scalar &operator()(const size_t &row, const size_t &col) const;

    Matrix<Scalar> block(size_t startRow, size_t startCol, size_t blockRow, size_t blockCol);

    size_t get_rows() const;
    size_t get_cols() const;
    void print_info() const;
};
