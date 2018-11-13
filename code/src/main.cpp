/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 24.03.2018
 */

#include <cstdlib>
#include <string>
#include <iostream>
#include <algorithm>
#include "Liebman.hpp"
#include "MatrixUtils.hpp"
#include "utils.hpp"
#include "PlotPy.hpp"
#include <chrono>
#include <fstream>
#include "Spline.hpp"

#include <boost/program_options.hpp>
#include <boost/type_traits/is_complex.hpp>

// #include <Exception.hpp>

namespace po = boost::program_options;

const int TEST_SIZE = 200;

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

  try
  {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")                                                                                //
        ("test", "produce a sample test with different options to show times")                                                        //                                                                                                                 //
        ("horizontal,h", po::value<int>()->default_value(512), "Indicates the number of horizontal pixels to use in the calculation") //
        ("vertical,v", po::value<int>()->default_value(512), "Indicates the number of vertical pixels to use in the calculation")     //

        ("top,t", po::value<std::vector<double>>(), "Indicates top border temperature in degrees Celsius")                                //
        ("bottom,b", po::value<std::vector<double>>(), "Indicates bottom border temperature in degrees Celsius")                          //
        ("left,l", po::value<std::vector<double>>(), "Indicates left border temperature in degrees Celsius")                              //
        ("right,r", po::value<std::vector<double>>(), "Indicates right border temperature in degrees Celsius")                            //
        ("error,e", po::value<double>()->default_value(0.001), "Desired error difference used for calculation")                           //
        ("lambda,a", po::value<double>()->default_value(1.3), "Desired relaxation step for Liebman. Must be between 1 and 2")             //
        ("isolate,i", po::value<std::string>(), "Which border to isolate, i.e top, bottom, left or right")                                //
        ("parallel_optimization,o", po::value<bool>()->default_value(true), "Option to use paralellism")                                  //
        ("simple_Liebman,s", po::value<bool>()->default_value(false), "Option to use the simple Liebman starting with a 0 filled Matrix") //
        (",q", "Deactivate visualization")                                                                                                //
        ("profile,p", po::value<std::string>(), "Indicates file to use for the termical profile [NOT WORKING]")                           //

        ("flux,f", "Show calculated heat flux vectors [NOT WORKING]") //
        ("grid,g", po::value<double>(), "Indicates the amount of pixels per grid cell to use in order to show the heat flux vectors [NOT WORKING]");

    po::variables_map vm;

    po::store(parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
      std::cout << desc << "\n\n";
      std::cout
          << "The program runs faster when sizes are a power of 2\n"
             "The program will output the matrices it used to do the calculations\n"
             "as well as the number of iterations it took on each of the matrices \n"
             "The program will also display the time it took to calcualte the solution\n"

          << std::endl;

      return EXIT_SUCCESS;
    }

    po::notify(vm);

    // Default values
    int h = 0, v = 0;
    bool parallel = true, simpleLieb = false;
    double error = 0.001, lambda = 1.3;
    anpi::Edge top, bot, left, right;

    //if no temp is given edges are isolated
    top.isolated = true;
    bot.isolated = true;
    right.isolated = true;
    left.isolated = true;

    if (vm.count("test"))
    {

      top.isolated = false;
      bot.isolated = false;
      right.isolated = false;
      left.isolated = false;
      top.temperatures.resize(TEST_SIZE, 50);
      bot.temperatures.resize(TEST_SIZE, 100);
      left.temperatures.resize(TEST_SIZE, 0);
      right.temperatures.resize(TEST_SIZE, 30);
      double timeSec;

      std::cout << "---------------------Simple Liebman (without parallelism, nor piramid calculation)---------------------------" << std::endl;
      anpi::LiebmnanSolver ls(top, bot, right, left, TEST_SIZE, TEST_SIZE, 0.001, 1.35, false, true);
      timeSec = ls.lieb();
      std::cout << "El método tardó: " << timeSec << " segundos\n\n"
                << std::endl;

      std::cout << "---------------------Simple Liebman with parallelism (no piramid calculation)---------------------------" << std::endl;
      anpi::LiebmnanSolver ls1(top, bot, right, left, TEST_SIZE, TEST_SIZE, 0.001, 1.35, true, true);
      timeSec = ls1.lieb();
      std::cout << "El método tardó: " << timeSec << " segundos\n\n"
                << std::endl;
      std::cout << "---------------------Simple Liebman with Piramid calculation (without parallelism)--------------------------" << std::endl;
      anpi::LiebmnanSolver ls2(top, bot, right, left, TEST_SIZE, TEST_SIZE, 0.001, 1.35, false, false);
      timeSec = ls2.lieb();
      std::cout << "El método tardó: " << timeSec << " segundos\n\n"
                << std::endl;
      std::cout << "---------------------Optimized Liebman with Piramid calculation & parallelism--------------------------" << std::endl;
      anpi::LiebmnanSolver ls3(top, bot, right, left, TEST_SIZE, TEST_SIZE, 0.001, 1.35, false, false);
      timeSec = ls3.lieb();
      std::cout << "El método tardó: " << timeSec << " segundos\n\n"
                << std::endl;

      return EXIT_SUCCESS;

      // anpi::LiebmnanSolver ls
    }

    if (vm.count("parallel_optimization"))
    {
      parallel = vm["parallel_optimization"].as<bool>();
    }
    if (vm.count("simple_Liebman"))
    {
      simpleLieb = vm["simple_Liebman"].as<bool>();
    }

    if (vm.count("horizontal"))
    {
      h = vm["horizontal"].as<int>();
      if (h <= 0)
      {
        std::cout << "Error in Horizontal: Matrix has to be bigger that 1.\n Provided value  " << h << std::endl;
        return EXIT_FAILURE;
      }
    }

    if (vm.count("vertical"))
    {
      v = vm["vertical"].as<int>();

      if (v <= 0)
      {
        std::cout << "Error in Vertical: Matrix has to be bigger that 1.\n Provided value:  " << v << std::endl;
        return EXIT_FAILURE;
      }
    }

    if (vm.count("lambda"))
    {
      lambda = vm["lambda"].as<double>();
      if (lambda <= 1 || lambda >= 2)
      {
        std::cout << "Error in Lamdba value: Value must be between 1 and 2.\n Provided value:  " << lambda << std::endl;
        return EXIT_FAILURE;
      }
    }

    if (vm.count("error"))
    {
      error = vm["error"].as<double>();
    }

    if (vm.count("profile"))
    {
      std::string filename = vm["profile"].as<std::string>();
      try
      {
        if (!anpi::readTempFile(filename, top, bot, right, left, v, h))
        {
          std::cout << "Error reading profile file: " << filename << "   Ignoring the file.  " << std::endl;
          // std::cout << "Received value cannot be made into double.  Received value: " << curr << std::endl;
        }
      }
      catch (anpi::Exception e)
      {
        std::cout << "Error reading profile file: " << filename << "   Ignoring the file.  " << std::endl;
        std::cerr << e.what() << '\n';
      }
    }

    if (vm.count("top"))
    {
      std::vector<double> temps = vm["top"].as<std::vector<double>>();
      if (temps.size() > 2)
      {
        //spline<T>(SizeVecBordes, VecTemperaturas, VecBordesValues);
        std::cout << " Error in top temperatures: Cannot have more than 2 values will use only first two" << std::endl;
      }
      if (temps.size() == 1)
      {
        double ttemp = temps[0];
        top.temperatures.resize(h, ttemp);
        top.isolated = false;
      }
      else
      {
        top.temperatures = anpi::linspace(temps[0], temps[1], h);
        top.isolated = false;
      }
    }

    if (vm.count("bottom"))
    {
      std::vector<double> temps = vm["bottom"].as<std::vector<double>>();
      if (temps.size() > 2)
      {
        std::cout << " Error in bottom temperatures: Cannot have more than 2 values will use only first two" << std::endl;
      }
      if (temps.size() == 1)
      {
        // double b = (double)vm["-b"].as<int>();
        double btemp = temps[0];
        bot.temperatures.resize(h, btemp);
        bot.isolated = false;
      }
      else
      {
        bot.temperatures = anpi::linspace(temps[0], temps[1], h);
        bot.isolated = false;
      }
    }

    if (vm.count("left"))
    {
      std::vector<double> temps = vm["left"].as<std::vector<double>>();
      if (temps.size() > 2)
      {
        std::cout << " Error in left temperatures: Cannot have more than 2 values will use only first two" << std::endl;
      }
      if (temps.size() == 1)
      {
        double ltemp = temps[0];
        left.temperatures.resize(h, ltemp);
        left.isolated = false;
      }
      else
      {
        left.temperatures = anpi::linspace(temps[0], temps[1], v);
        left.isolated = false;
      }
    }

    if (vm.count("right"))
    {
      std::vector<double> temps = vm["right"].as<std::vector<double>>();
      if (temps.size() > 2)
      {
        std::cout << " Error in right temperatures: Cannot have more than 2 values will use only first two" << std::endl;
      }
      if (temps.size() == 1)
      {

        double rtemp = temps[0];
        right.temperatures.resize(h, rtemp);
        right.isolated = false;
      }
      else
      {
        right.temperatures = anpi::linspace(temps[0], temps[1], v);
        right.isolated = false;
      }
    }

    if (vm.count("isolate"))
    {
      std::string isolateString = vm["isolate"].as<std::string>();
      int stringSize = isolateString.size();
      if (stringSize == 0)
      {
        std::cout << "Error in Isolate flag: if flag is present, at least one side must be present" << std::endl;
        return EXIT_FAILURE;
      }
      if (stringSize > 4)
      {
        std::cout << "Error in Isolate flag: there can't be more than 4 isolated sides" << std::endl;
        return EXIT_FAILURE;
      }

      for (int i = 0; i < stringSize; ++i)
      {
        char curr = isolateString.at(i);
        switch (curr)
        {
        case ('t'):
          top.isolated = true;
          break;
        case ('l'):
          left.isolated = true;
          break;
        case ('b'):
          bot.isolated = true;
          break;
        case ('r'):
          right.isolated = true;
          break;

        default:
          std::cout << "Error in Isolate flag: unknown character side added: " << curr << std::endl;
          std::cout << "Possible characters are t for top, b for bottom, r right and l for left' " << std::endl;

          return EXIT_FAILURE;
        }
      }
    }

    /////////////////////////////////////////////LIEBAMN CALL///////////////////////////////////////
    //do the liebman calculations
    anpi::LiebmnanSolver ls(top, bot, right, left, v, h, error, lambda, parallel, simpleLieb);

    std::cout << "\nEl método Liebmann:\n " << std::endl;

    double timeSec = ls.lieb();

    std::cout << "El método tardó: " << timeSec << " segundos" << std::endl;

    // std::cout << "\nEl método Liebman normal: \n " << std::endl;
    // anpi::Matrix<double> m(v, h);
    // t_start = std::chrono::high_resolution_clock::now();
    // auto it = ls.liebman(m, top, bot, right, left);
    // t_end = std::chrono::high_resolution_clock::now();
    // timeSec = ((std::chrono::duration<double, std::milli>(t_end - t_start).count()) / 1000);
    // std::cout << "iterations: " << it
    //           << std::endl;
    // // anpi::printMatrix(ls.tempsMatrix);
    // std::cout << "El método tardo: " << timeSec << " segundos" << std::endl;

    // ls.tempsMatrix.DumpToFile("matrizFinal.txt");

    //end calculations with no visualization

    if (vm.count("-q"))
    {

      return EXIT_SUCCESS;
    }

    //show color visuals

    //show visual grid
    anpi::Plot2d<double> plotter;
    plotter.initialize();
    std::string title = "Matriz de ";
    title.append(std::to_string(v) + std::string("x")).append(std::to_string(h));
    // plotter.imgshow(ls.tempsMatrix, title);

    if (vm.count("flux"))
    {
      std::cout << "se mete a FLUX\n";
      plotter.quiver(ls.tempsMatrix);
      plotter.show();
    }
    if (vm.count("grid"))
    {
    }

    plotter.show();

    return EXIT_SUCCESS;

    //show heat flux vectors

    // Dispatch with the proper types to call the real workers
  }
  catch (po::error &e)
  {
    std::cerr << "Error:\n  " << e.what() << std::endl
              << std::endl;
    return EXIT_FAILURE;
  }
  catch (std::exception &e)
  {
    std::cerr << "Unhandled Exception reached the top of main:\n  "
              << e.what() << "\nApplication will now exit" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
