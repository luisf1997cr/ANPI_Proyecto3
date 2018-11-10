/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>

#include "Liebman.hpp"
#include "MatrixUtils.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
// #include <complex>

#include <functional>

#include <cmath>

namespace anpi
{
namespace test
{

template <typename T>
void simpleLiebTest()
{
  Matrix<double> m(2, 2);
  // m.fill(0);

  Edge top, bot, left, right;
  top.temperatures = {5, 5};
  bot.temperatures = {0, 0};
  right.temperatures = {1, 1};
  left.temperatures = {5, 5};

  int iter;

  LiebmnanSolver lieb(top, bot, right, left, m);

  iter = lieb.liebman();

  anpi::printMatrix(m);
  std::cout << "\nNumber of iterations it took: " << iter << std::endl;
  anpi::printMatrix(lieb.temperatureMatrix);

  // BOOST_CHECK(Ai == exAi);
}

template <typename T>
void prueba()
{
  std::cout << (9 / 8) << std::endl;
}

template <typename T>
void solverTest()
{

  // BOOST_CHECK();
}

} // namespace test
} // namespace anpi

BOOST_AUTO_TEST_SUITE(Liebman)

BOOST_AUTO_TEST_CASE(simpleLieb)
{
  anpi::test::simpleLiebTest<int>();
}

// BOOST_AUTO_TEST_CASE(Crout)
// {
// }

BOOST_AUTO_TEST_SUITE_END()
