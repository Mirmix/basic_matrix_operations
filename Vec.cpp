#include "Vec.h"
#include <cmath>

template <typename Scalar>
Vector<Scalar>::Vector(size_t rows, Scalar initial) : Matrix<Scalar>{rows, 1, initial}
{
}

template <typename Scalar>
Vector<Scalar>::Vector(const Matrix<Scalar> &rhs) : Matrix<Scalar>{rhs}
{
    assert(rhs.get_cols() == 1);
}

template <typename Scalar>
Scalar Vector<Scalar>::norm()
{
    Scalar norm = 0.;
    for (size_t i = 0; i < this->get_rows(); i++)
    {
        norm += (*this)(i) * (*this)(i);
    }
    norm = sqrt(norm);
    return norm;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::normalize()
{
    Vector<Scalar> result = (*this) / this->norm();
    (*this) = result;
    return (*this);
}

template <typename Scalar>
Scalar Vector<Scalar>::dotProduct(Vector<Scalar> &rhs)
{
    Matrix<Scalar> result = (this->transpose()) * (rhs);
    assert(result.get_cols() == 1);
    assert(result.get_rows() == 1);
    return result(0, 0);
}

template <typename Scalar>
Scalar &Vector<Scalar>::operator()(const size_t &row)
{
    return Matrix<Scalar>::operator()(row, 0);
}

template <typename Scalar>
const Scalar &Vector<Scalar>::operator()(const size_t &row) const
{
    return Matrix<Scalar>::operator()(row, 0);
}

template <typename Scalar>
Matrix<Scalar> Vector<Scalar>::hat()
{
    assert(this->get_cols() == 1);
    assert(this->get_rows() == 3);
    Matrix<Scalar> result(3, 3, 0);
    result(0, 1) = Scalar(-1) * ((*this)(2));
    result(0, 2) = ((*this)(1));
    result(1, 2) = Scalar(-1) * ((*this)(0));
    result(1, 0) = result(0, 1) * (-1);
    result(2, 0) = result(0, 2) * (-1);
    result(2, 1) = result(1, 2) * (-1);
    return result;
}

template <typename Scalar>
Vector<Scalar> Vector<Scalar>::crossProduct(Vector<Scalar> &rhs)
{
    assert(this->get_cols() == 1);
    assert(this->get_rows() == 3);
    assert(rhs.get_cols() == 1);
    assert(rhs.get_rows() == 3);

    return Vector<Scalar>((this->hat()) * rhs);
}