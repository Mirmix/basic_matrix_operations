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
        norm += (*this)(i, 0) * (*this)(i, 0);
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
    return result(0,0);
}
