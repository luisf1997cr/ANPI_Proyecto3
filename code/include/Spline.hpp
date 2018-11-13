#ifndef ANPI_SPLINE_HPP
#define ANPI_SPLINE_HPP

#include "Thomas.hpp"

namespace anpi
{

/** 
   * Encuentra los valores de X dado una vector de Y
   * @param SizeX
   * @param VecY
   * @param VecX
   */
template <typename T>
void VecXFiller(const int SizeX, std::vector<T> &VecX, const std::vector<T> &VecY)
{

  T spacing = T(SizeX) / T(VecY.size() - 1);

  T x = T(0);

  VecX.clear();

  while (VecX.size() < VecY.size())
  {

    VecX.push_back(x);
    if (x + spacing >= SizeX)
    {
      x = SizeX - 1;
    }
    else
    {
      x += spacing;
    }
  }
}

/** 
   * Completa los valores en los vectores x y y dados una funcion A = xB
   * @tparam T
   * @param VecX
   * @param VecY
   * @param A
   * @param x
   * @param b
   */
template <typename T>
void operandFiller(const std::vector<T> &VecX,
                   const std::vector<T> &VecY,
                   Matrix<T> &A,
                   std::vector<T> &x,
                   std::vector<T> &b)
{

  int SizeResult = VecY.size() - 2;

  A = Matrix<T>(SizeResult, SizeResult);

  if (SizeResult == 1)
  {
    A[0][0] = 2 * (VecX[2] - VecX[0]);
  }

  else
  {

    A[0][0] = 2 * (VecX[2] - VecX[0]);
    A[0][1] = VecX[2] - VecX[1];

    for (int i = 1; i < (SizeResult - 1); i++)
    {
      A[i][i - 1] = VecX[i + 1] - VecX[i];
      A[i][i] = 2 * (VecX[i + 2] - VecX[i]);
      A[i][i + 1] = VecX[i + 2] - VecX[i + 1];
    }

    A[SizeResult - 1][SizeResult - 2] = VecX[SizeResult + 1] - VecX[SizeResult];
    A[SizeResult - 1][SizeResult - 1] = 2 * (VecX[SizeResult + 2] - VecX[SizeResult]);
  }

  x.resize(SizeResult, T(1));

  b.clear();

  for (int i = 0; i < SizeResult; i++)
  {
    b.push_back(6 * ((VecY[i + 2] - VecY[i + 1]) / (VecX[i + 2] - VecX[i + 1])) -
                6 * ((VecY[i + 1] - VecY[i]) / (VecX[i + 1] - VecX[i])));
  }
}

/** 
   *Metodo utilizado para completar los valor de x dado un F(x) o su derivada
   * @tparam T
   * @param SizeX size of iterpolatedValues vector
   * @param VecX vector of x
   * @param VecY vector of f(x)
   * @param x     vector of f''(x)
   * @param interpolatedValues
   */
template <typename T>
void interpolate(const int SizeX,
                 const std::vector<T> &VecX,
                 const std::vector<T> &VecY,
                 const std::vector<T> &x,
                 std::vector<T> &interpolatedValues)
{

  int index = 1;

  for (int i = 1; i <= SizeX; i++)
  {

    if (i > VecX[index] && (index + 1) < int(VecX.size()))
    {
      index++;
    }

    interpolatedValues.push_back(x[index - 1] * (((i - 1) - VecX[index]) * ((i - 1) - VecX[index]) * ((i - 1) - VecX[index])) / (6 * (VecX[index - 1] - VecX[index])) +
                                 x[index] * (((i - 1) - VecX[index - 1]) * ((i - 1) - VecX[index - 1]) * ((i - 1) - VecX[index - 1])) / (6 * (VecX[index] - VecX[index - 1])) +
                                 (((VecY[index - 1]) / (VecX[index - 1] - VecX[index])) - ((x[index - 1] * (VecX[index - 1] - VecX[index])) / 6)) * (i - 1 - VecX[index]) +
                                 (((VecY[index]) / (VecX[index] - VecX[index - 1])) - ((x[index] * (VecX[index] - VecX[index - 1])) / 6)) * (i - 1 - VecX[index - 1]));
  }
}

/** Realiza A=xB cuando la ecuacion cuenta con 3 o mas bordes
     * @tparam T
     * @param Tamaño de vector X
     * @param Valores de Y
     * @return void
     */
template <typename T>
void spline(const int vecxSize, const std::vector<T> &VecY, std::vector<T> &interpolValue)
{

  std::vector<T> VecX;

  VecXFiller(vecxSize, VecX, VecY);

  Matrix<T> A;
  std::vector<T> x;
  std::vector<T> b;

  operandFiller(VecX, VecY, A, x, b);

  thomas(A, x, b);

  interpolate(SizeX, VecX, VecY, x, interpolatedValues);
}

} // namespace anpi

#endif