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
#include "utils.hpp"
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

void simpleLiebTest()
{
  Matrix<double> m(2, 2);
  // m.fill(0);
  std::cout << "\nTest cuadrada " << std::endl;
  Edge top, bot, left, right;
  top.temperatures = {5, 5};
  bot.temperatures = {0, 0};
  right.temperatures = {1, 1};
  left.temperatures = {5, 5};

  int iter;

  LiebmnanSolver lieb(top, bot, right, left);

  iter = lieb.liebman(m, top, bot, right, left);

  anpi::printMatrix(m);

  lieb.fillBiggerTempMatrix(m, 3, 3);

  anpi::printMatrix(m);

  std::cout << "\nNumber of iterations it took: " << iter << std::endl;
  // anpi::printMatrix(lieb.temperatureMatrix);

  Matrix<double> n(1, 1);

  top.temperatures = {5};
  bot.temperatures = {0};
  right.temperatures = {1};
  left.temperatures = {5};

  LiebmnanSolver slieb(top, bot, right, left); //, n);
  std::cout
      << "\nLiebman Matriz tamaño 1" << std::endl;
  iter = slieb.liebman(n, top, bot, right, left);
  anpi::printMatrix(n);
  std::cout << "\nmatriz duplicada" << std::endl;
  slieb.duplicateTempMatrix(n);
  anpi::printMatrix(n);
  std::cout << "\nLAlto duplicado" << std::endl;
  slieb.duplicateHeightTempMatrix(n);
  anpi::printMatrix(n);
  std::cout << "\nLAanchoduplicado" << std::endl;
  slieb.duplicateWidthTempMatrix(n);

  anpi::printMatrix(n);
  std::cout << "\nNumber of iterations it took: " << iter << std::endl;

  // BOOST_CHECK(Ai == exAi);
}

void simpleLiebTestRect()
{
  Matrix<double> m(3, 4);
  // m.fill(0);

  Edge top, bot, left, right;
  top.temperatures = {5, 5, 5, 5};
  bot.temperatures = {0, 0, 0, 0};
  right.temperatures = {1, 1, 1};
  left.temperatures = {5, 5, 5};

  int iter;

  LiebmnanSolver lieb(top, bot, right, left); //, m);

  iter = lieb.liebman(m, top, bot, right, left);
  std::cout << "\nTest rectangular" << std::endl;
  anpi::printMatrix(m);
  std::cout << "\nNumber of iterations it took: " << iter << std::endl;
  // anpi::printMatrix(lieb.tempsMatrix);

  // BOOST_CHECK(Ai == exAi);
}

void liebmanSeriousTests()
{
  //temperature variables
  Edge top, bot, left, right;
  int numMatrices;

  //test 512x512
  std::cout << "\nTest 512x512" << std::endl;

  //create temperatures
  top.temperatures.resize(512, 50);
  bot.temperatures.resize(512, 0);
  right.temperatures.resize(512, 25);
  left.temperatures.resize(512, 10);

  LiebmnanSolver lieb(top, bot, right, left);

  numMatrices = lieb.lieb();

  // anpi::printMatrix(lieb.tempsMatrix);
  std::cout << "\nNumber of iterations it took: " << numMatrices << std::endl;

  //test 512x520
  std::cout << "\nTest 512x520" << std::endl;

  //create temperatures
  top.temperatures.resize(512, 50);
  bot.temperatures.resize(512, 0);
  right.temperatures.resize(520, 25);
  left.temperatures.resize(520, 10);

  LiebmnanSolver rlieb(top, bot, right, left);

  numMatrices = rlieb.lieb();
  // anpi::printMatrix(lieb.tempsMatrix);
  std::cout << "\nNumber of iterations it took: " << numMatrices << std::endl;

  //test 520x512
  std::cout << "\nTest 520x512" << std::endl;

  //create temperatures
  top.temperatures.resize(520, 50);
  bot.temperatures.resize(520, 0);
  right.temperatures.resize(512, 25);
  left.temperatures.resize(512, 10);

  LiebmnanSolver llieb(top, bot, right, left);

  numMatrices = llieb.lieb();
  // anpi::printMatrix(lieb.tempsMatrix);
  std::cout << "\nNumber of iterations it took: " << numMatrices << std::endl;

  //test 515x520
  std::cout << "\nTest 5152x520" << std::endl;

  //create temperatures
  top.temperatures.resize(515, 50);
  bot.temperatures.resize(515, 0);
  right.temperatures.resize(520, 25);
  left.temperatures.resize(520, 10);

  LiebmnanSolver tlieb(top, bot, right, left);

  numMatrices = tlieb.lieb();
  // anpi::printMatrix(lieb.tempsMatrix);
  std::cout << "\nNumber of iterations it took: " << numMatrices << std::endl;

  //test 515x520 with linspace temps
  std::cout << "\nTest 515x520 with linspace" << std::endl;

  //create temperatures
  top.temperatures = anpi::linspace<double>(0, 50, 515);
  bot.temperatures = anpi::linspace<double>(0, 70, 515);
  right.temperatures = anpi::linspace<double>(0, 100, 520);
  left.temperatures = anpi::linspace<double>(0, 17, 520);

  LiebmnanSolver linlieb(top, bot, right, left);

  numMatrices = linlieb.lieb();
  // anpi::printMatrix(lieb.tempsMatrix);
  std::cout << "\nNumber of iterations it took: " << numMatrices << std::endl;
}

struct ank
{
  int j = 0;
};
ank f()
{
  ank thisone;
  thisone.j = 23;
  return thisone;
}
void prueba()
{
  std::cout << (9 / 8) << std::endl;

  Matrix<int> a(2, 2, 1);
  Matrix<int> b(5, 6);

  ank test = f();
  std::cout << test.j << std::endl;

  std::vector<int> vec;
  vec.resize(3, 5);
  for (int i = 0; i < vec.size(); i++)
    std::cout << vec[i] << std::endl;

  // b.fill(a);

  // anpi::printMatrix(b);
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
  // anpi::test::simpleLiebTest();
  // anpi::test::simpleLiebTestRect();
  // anpi::test::prueba();
  anpi::test::liebmanSeriousTests();
}

// BOOST_AUTO_TEST_CASE(Crout)
// {
// }

BOOST_AUTO_TEST_SUITE_END()
