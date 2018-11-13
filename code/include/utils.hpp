#ifndef ANPI_UTILS_HPP
#define ANPI_UTILS_HPP

#include <cstdlib>
#include <vector>

namespace anpi
{
template <typename T>
std::vector<T> linspace(T a, T b, size_t N)
{
    T h = (b - a) / static_cast<T>(N - 1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

/**
* Obtiene la temperatura de los bordes en caso de que no 
* se den los cuatro bordes
* @tparam T
* @param SizeVecBordes
* @param VecTemperaturas
* @param VecBordesValues
*/
template<typename T>
void obtainVecBordesValues(const int SizeVecBordes              ,
                          const std::vector<T>& VecTemperaturas,
                          std::vector<T>& VecBordesValues) {

    if (VecTemperaturas.size() >= 3) {
      spline<T>(SizeVecBordes, VecTemperaturas, VecBordesValues);
      std::cout << "using splines \n";
    }

    else if (VecTemperaturas.size() == 2) {
      linearIncrement<T>(SizeVecBordes, VecTemperaturas, VecBordesValues);
      std::cout << "lineal increment \n";
    }

    else if (VecTemperaturas.size() == 1) {
      for (int i = 0; i < SizeVecBordes; i++) {
        VecBordesValues.push_back(VecTemperaturas[0]);
      }
      std::cout << "constant values \n";
    }

    else {
      std::cout << "border isolated \n";
    }

  }
} // namespace anpi

#endif