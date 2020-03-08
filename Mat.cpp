#include "iostream"
#include "Mat.h"
#include <cassert>

template <typename Scalar>
Matrix<Scalar>::Matrix(size_t rows, size_t cols, const Scalar &initial)
{
    assert(rows >= 0);
    assert(cols >= 0);

    rows_ = rows;
    cols_ = cols;

    mat.resize(rows_);

    for (size_t i = 0; i < mat.size(); i++)
        mat[i].resize(cols_, initial);
}

template <typename Scalar>
Matrix<Scalar>::Matrix(const Matrix<Scalar> &rhs)
{
    rows_ = rhs.get_rows();
    cols_ = rhs.get_cols();
    mat = rhs.mat;
}

template <typename Scalar>
Matrix<Scalar>::~Matrix()
{
    rows_ = 0;
    cols_ = 0;
    mat.clear();
}

template <typename Scalar>
Scalar &Matrix<Scalar>::operator()(const size_t &row, const size_t &col)
{
    return this->mat[row][col];
}

template <typename Scalar>
const Scalar &Matrix<Scalar>::operator()(const size_t &row, const size_t &col) const
{
    return this->mat[row][col];
}

template <typename Scalar>
size_t Matrix<Scalar>::get_rows() const
{
    return this->rows_;
}

template <typename Scalar>
size_t Matrix<Scalar>::get_cols() const
{
    return this->cols_;
}


template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator=(const Matrix<Scalar> &rhs)
{
    size_t new_rows = rhs.get_rows();
    size_t new_cols = rhs.get_cols();

    mat.resize(new_rows);
    for (size_t i = 0; i < mat.size(); i++)
        mat[i].resize(new_cols);

    for (size_t i = 0; i < new_rows; i++)
        for (size_t j = 0; j < new_cols; j++)
            mat[i][j] = rhs(i, j);

    rows_ = new_rows;
    cols_ = new_cols;

    return *this;
}

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator+(const Matrix<Scalar> &rhs)
{
    assert(rows_ == rhs.get_rows());
    assert(cols_ == rhs.get_cols());

    Matrix result(rows_, cols_, Scalar(0));

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            result(i, j) = this->mat[i][j] + rhs(i, j);

    return result;
}

template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator+=(const Matrix<Scalar> &rhs)
{
    assert(rows_ == rhs.get_rows());
    assert(cols_ == rhs.get_cols());

    size_t rows_ = rhs.get_rows();
    size_t cols_ = rhs.get_cols();

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            this->mat[i][j] += rhs(i, j);

    return *this;
}

// Subtraction of this matrix and another
template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator-(const Matrix<Scalar> &rhs)
{
    assert(rows_ == rhs.get_rows());
    assert(cols_ == rhs.get_cols());

    size_t rows_ = rhs.get_rows();
    size_t cols_ = rhs.get_cols();
    Matrix result(rows_, cols_, Scalar(0));

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            result(i, j) = this->mat[i][j] - rhs(i, j);

    return result;
}

template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator-=(const Matrix<Scalar> &rhs)
{
    assert(rows_ == rhs.get_rows());
    assert(cols_ == rhs.get_cols());

    size_t rows_ = rhs.get_rows();
    size_t cols_ = rhs.get_cols();

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            this->mat[i][j] -= rhs(i, j);

    return *this;
}

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator*(const Matrix<Scalar> &rhs)
{
    assert(cols_ == rhs.get_rows());

    size_t rhs_rows = rhs.get_rows();
    size_t rhs_cols = rhs.get_cols();
    Matrix result(rows_, rhs_cols, Scalar(0));

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            for (size_t k = 0; k < rhs_cols; k++)
                result(i, k) = result(i, k) + ((this->mat[i][j]) * rhs(j, k));

    return result;
}

template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator*=(const Matrix<Scalar> &rhs)
{
    Matrix result = (*this) * rhs;
    (*this) = result;
    return *this;
}

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::transpose()
{
    Matrix result(cols_, rows_, Scalar(0));

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            result(j, i) = this->mat[i][j];

    return result;
}

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator+(const Scalar &rhs)
{
    Matrix result(rows_, cols_, Scalar(0));

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            result(i, j) = this->mat[i][j] + rhs;

    return result;
}

template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator+=(const Scalar &rhs)
{
    Matrix result = (*this) + rhs;
    (*this) = result;
    return *this;
}

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator-(const Scalar &rhs)
{
    return ((*this) + (-1. * rhs));
}


template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator-=(const Scalar &rhs)
{
    Matrix result = (*this) - rhs;
    (*this) = result;
    return *this;
}

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator*(const Scalar &rhs)
{
    Matrix result(rows_, cols_, Scalar(0));

    for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
            result(i, j) = this->mat[i][j] * rhs;

    return result;
}

template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator*=(const Scalar &rhs)
{
    Matrix result = (*this) * rhs;
    (*this) = result;
    return *this;
}


template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator/(const Scalar &rhs)
{

    return ((*this) * (1. / rhs));
}

template <typename Scalar>
Matrix<Scalar> &Matrix<Scalar>::operator/=(const Scalar &rhs)
{
    Matrix result = (*this) / rhs;
    (*this) = result;
    return *this;
}

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::block(size_t startRow, size_t startCol, size_t blockRow, size_t blockCol)
{
    assert(startRow >= 0 && startCol >= 0);
    assert(startRow + blockRow <= rows_ && startCol + blockCol <= cols_);
    Matrix result(blockRow, blockCol, Scalar(0));
    for (size_t i = 0; i < blockRow; i++)
        for (size_t j = 0; j < blockCol; j++)
            result(i, j) = this->mat[i + startRow][j + startCol];
    return result;
}


template <typename Scalar>
void Matrix<Scalar>::print_info() const
{
    std::cout << "Rows = " << this->get_rows() << " Cols = " << this->get_cols() << std::endl;
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
            std::cout << (*this)(i, j) << " ";
        std::cout << std::endl;
    }
}