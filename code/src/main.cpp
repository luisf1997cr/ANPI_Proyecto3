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

#include <boost/program_options.hpp>
#include <boost/type_traits/is_complex.hpp>

// #include <Exception.hpp>

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////
// int main(int argc, char *argv[])
// {
//   return 0;
// }
int main(int argc, char *argv[])
{

  try
  {
    po::options_description desc("Allowed options");
    desc.add_options()                                                                                                                //
        ("help", "produce help message")                                                                                              //
        ("top,t", po::value<std::vector<double>>(), "Indicates top border temperature in degrees Celsius")                            //
        ("bottom,b", po::value<std::vector<double>>(), "Indicates bottom border temperature in degrees Celsius")                      //
        ("left,l", po::value<std::vector<double>>(), "Indicates left border temperature in degrees Celsius")                          //
        ("right,r", po::value<std::vector<double>>(), "Indicates right border temperature in degrees Celsius")                        //
        ("error,e", po::value<double>()->default_value(0.001), "Desired error percentage used for calculation")                       //
        ("lambda,a", po::value<double>()->default_value(1.3), "Desired relaxation step for Liebman. Must be between 1 and 2")         //
        ("isolate,i", po::value<std::string>(), "Which border to isolate, i.e top, bottom, left or right")                            //
        ("profile,p", po::value<std::string>(), "Indicates file to use for the termical profile")                                     //
        ("horizontal,h", po::value<int>()->default_value(512), "Indicates the number of horizontal pixels to use in the calculation") //
        ("vertical,v", po::value<int>()->default_value(512), "Indicates the number of vertical pixels to use in the calculation")     //
        (",q", "Deactivate visualization")                                                                                            //
        ("flux,f", "Show calculated heat flux vectors")                                                                               //
        ("grid,g", po::value<double>(), "Indicates the amount of pixels per grid cell to use in order to show the heat flux vectors");

    po::variables_map vm;

    po::store(parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
      std::cout << desc << "\n\n";
      // std::cout
      //     << "For the coefficient and root types, use one of the following:\n"
      //        "  float    real numbers with single precision \n"
      //        "  double   real numbers with double precision \n"
      //        "  fcomplex complex numbers with float components\n"
      //        "  dcomplex complex numbers with double components\n\n"
      //     << "The polynomial formulae are composed by one or more \n"
      //        "polynomial terms:\n"
      //        "  polynomial  := term [ {'+'|'-'} term ]*\n"
      //        "  term        := coefficient 'x^' exponent\n"
      //        "  coefficient := double | complex\n"
      //        "  complex     := { '(' double ',' double ')' } | (double)\n"
      //        "  exponent    := unsigned integral\n\n"
      //     << "Examples of valid polynomials:\n"
      //        "  x^3 + 2x + 1\n"
      //        "  -5x^4 + 2.5x^2 + x\n"
      //        "  -3x + 5x^4 + 1 -2x^2\n"
      //        "  (0,1)x^4 + (5,2)x^2 + 1.5x + (1.5,2)\n"
      //        "The last example has some complex coefficients"
      //     << std::endl;

      return EXIT_SUCCESS;
    }

    po::notify(vm);

    // std::string poly = vm["poly"].as<std::string>();

    // Default values
    int h, v;
    double error, lambda;
    anpi::Edge top, bot, left, right;
    top.isolated = true;
    bot.isolated = true;
    right.isolated = true;
    left.isolated = true;

    if (vm.count("horizontal"))
    {
      h = vm["horizontal"].as<int>();
      std::cout << "agarro horizontal  " << vm["horizontal"].as<int>() << std::endl;
    }

    if (vm.count("vertical"))
    {
      v = vm["vertical"].as<int>();
      std::cout << "agarro vertical  " << vm["vertical"].as<int>() << std::endl;
    }

    if (vm.count("lambda"))
    {
      lambda = vm["lambda"].as<double>();
      // std::cout << "agarro horizontal  " << vm["horizontal"].as<int>() << std::endl;
    }
    if (vm.count("error"))
    {
      error = vm["error"].as<double>();
      // std::cout << "agarro horizontal  " << vm["horizontal"].as<int>() << std::endl;
    }

    if (vm.count("profile"))
    {
      //check file and parse
      std::cout << "agarro profile: " << vm["profile"].as<std::string>() << std::endl;
    }
    if (vm.count("top"))
    {
      std::vector<double> temps = vm["top"].as<std::vector<double>>();
      if (temps.size() > 2)
      {
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

      std::cout << "agarro isolate: " << vm["isolate"].as<std::string>() << std::endl;
    }

    //do the liebman calculations
    anpi::LiebmnanSolver ls(top, bot, right, left, v, h, error, lambda);
    ls.lieb();
    // anpi::printMatrix(ls.tempsMatrix);
    anpi::Plot2d<double> plotter;
    plotter.initialize();
    plotter.imgshow(ls.tempsMatrix);
    plotter.show();
    //end calculations with no visualization

    if (vm.count("-q"))
    {

      return EXIT_SUCCESS;
    }

    //show color visuals

    if (vm.count("flux"))
    {
    }
    if (vm.count("grid"))
    {
    }

    // return EXIT_SUCCESS;

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
