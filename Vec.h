#include "Mat.cpp"

template <typename Scalar>
class Vector : public Matrix<Scalar>
{
public:
    Vector(size_t row, Scalar initial = Scalar(0));
    Vector(const Matrix<Scalar> &rhs);
    Scalar norm();
    Vector<Scalar> normalize();

    Scalar &operator()(const size_t &row);
    const Scalar &operator()(const size_t &row) const;

    Matrix<Scalar> hat();
    Scalar dotProduct(Vector<Scalar> &rhs);
    Vector<Scalar> crossProduct(Vector<Scalar> &rhs);
};