#include "Mat.cpp"


template <typename Scalar>
class Vector : public Matrix<Scalar>
{
public:
    Vector(size_t row, Scalar initial = Scalar(0));
    Vector(const Matrix<Scalar> &rhs);
    Scalar norm();
    Vector<Scalar> normalize();
    Scalar dotProduct(Vector<Scalar> &rhs);
};