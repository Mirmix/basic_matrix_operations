#include "iostream"
#include "Vec.cpp"

int main()
{
    //test
    Matrix<double> A(5, 5, 0);
    Matrix<double> B(4, 4, 0);
    Matrix<double> C(5, 5, 1);
    Vector<double> D(5, 0.);
    Vector<double> E(5, 1.);
    Vector<double> V3D(3, 1.);
    Vector<double> W3D(3, 0.);

    E.print_info();
    Vector<double>((E * 5)).normalize().print_info();
    std::cout << E.dotProduct(E) << std::endl;
    std::cout << Vector<double>((E * 5)).normalize().norm();
    C(2, 2) = 6;
    (C * E).print_info();
    std::cout << "(E*5)(3) = " << Vector<double>((E * 5))(3) << std::endl;
    W3D(2)=-1.;
    (V3D.crossProduct(W3D)).print_info();
}