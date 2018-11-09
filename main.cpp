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

#include <boost/program_options.hpp>
#include <boost/type_traits/is_complex.hpp>

// #include <Exception.hpp>

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{

  try
  {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("top,t", po::value<double>(), "Indicates top border temperature in degrees Celsius") //
        ("bottom,b", po::value<double>(), "Indicates bottom border temperature in degrees Celsius")                                         //
        ("left,l", po::value<double>(), "Indicates left border temperature in degrees Celsius")                                             //
        ("right,r", po::value<double>(), "Indicates right border temperature in degrees Celsius")                                           //
        ("isolate,i", po::value<std::string>(), "Which border to isolate, i.e top, bottom, left or right")                                  //
        ("profile,p", po::value<std::string>(), "Indicates file to use for the termical profile")                                           //
        ("horizontal,h", po::value<double>(), "Indicates the number of horizontal pixels to use in the calculation")                        //
        ("vertical,v", po::value<double>(), "Indicates the number of vertical pixels to use in the calculation")                            //
        (",q", "Deactivate visualization")                                                                                                  //
        ("flux,f", "Show calculated heat flux vectors")                                                                                     //
        ("grid,g", po::value<double>(), "Indicates the amount of pixels per grid cell to use in order to show the heat flux vectors");

    po::variables_map vm;

    po::store(parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
      std::cout
          << desc << "\n\n";
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

    std::string poly = vm["poly"].as<std::string>();

    // Default values

    if (vm.count("profile"))
    {
    }
    if (vm.count("top"))
    {
    }

    if (vm.count("bottom"))
    {
    }

    if (vm.count("left"))
    {
    }

    if (vm.count("right"))
    {
    }

    if (vm.count("isolate"))
    {
    }

    if (vm.count("horizontal"))
    {
    }

    if (vm.count("vertical"))
    {
    }

    //do the liebman calculations

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
