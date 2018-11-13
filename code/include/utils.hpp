#ifndef ANPI_UTILS_HPP
#define ANPI_UTILS_HPP

#include <cstdlib>
#include <vector>
#include "Spline.hpp"

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

std::vector<double> makeEdgeTemps(std::vector<double> inputs, int size)
{
    int cantInputs = inputs.size();
    if (cantInputs > size)
    {

        throw anpi::Exception("Error making edge temps, there are more inputs than the required size");
    }

    std::vector<double> temps = anpi::linspace(inputs[0], inputs[1], size);

    return temps;
}

bool readTempFile(std::string filename, Edge &top, Edge &bottom, Edge &right, Edge &left, int vertical, int horizontal)
{
    std::ifstream inFile;

    //check file and parse
    // std::cout << "agarro profile: " << vm["profile"].as<std::string>() << std::endl;

    inFile.open(filename);
    if (!inFile)
    {
        std::cout << "Error al abrir el archivo " << filename << std::endl;
        return false;
    }
    else
    {
        std::string line;
        if (inFile.is_open())
        {
            int l = 0;
            while (getline(inFile, line)) //for each line on the file
            {
                ++l;

                std::stringstream ss(line);
                std::vector<double> temps;

                std::string curr;
                char currSide;
                getline(ss, curr, ' '); //grab the first word on the line, must be top, bottom, right or left

                if (curr.find("top") || curr.find("bottom") || curr.find("right") || curr.find("left"))
                {
                    currSide = curr.at(0); //get the first character of the word

                    getline(ss, curr, ' ');
                    if (curr.find('=')) //next character is equals
                    {
                        while (getline(ss, curr, ' '))
                        {
                            // if(strcmp(curr,""))

                            try
                            {
                                temps.push_back(std::stod(curr));
                            }
                            catch (const std::exception &e)
                            {
                                std::cout << "Error in profile file: " << filename << "  at line (" << l << ")" << std::endl;
                                std::cout << "Received value cannot be made into double.  Received value: " << curr << std::endl;
                                std::cerr << e.what() << '\n';
                            }
                        }
                    }

                    switch (currSide)
                    {
                    case ('t'):
                        if (temps.size() > 1)
                            top.gradient = true;
                        top.temperatures = anpi::makeEdgeTemps(temps, horizontal);
                        top.isolated = false;
                        break;
                    case ('b'):
                        if (temps.size() > 1)
                            bottom.gradient = true;
                        bottom.temperatures = anpi::makeEdgeTemps(temps, horizontal);
                        bottom.isolated = false;

                        break;
                    case ('r'):
                        if (temps.size() > 1)
                            right.gradient = true;
                        right.temperatures = anpi::makeEdgeTemps(temps, vertical);
                        right.isolated = false;
                        break;
                    case ('l'):
                        if (temps.size() > 1)
                            left.gradient = true;
                        left.temperatures = anpi::makeEdgeTemps(temps, vertical);
                        left.isolated = false;
                        break;
                    default:
                        return false;
                    }
                }

                else
                {
                    // std::cout << "Line bust be lowercase and start with top, bottom, left or right " << std::endl;
                    std::cout << "Error in profile file: " << filename << "  at line (" << l << ")" << std::endl;
                    std::cout << "Line bust be lowercase and start with top, bottom, left or right " << std::endl;

                    return false;
                }
            }
        }
        inFile.close();
        return true;
    }
}

/**
* Obtiene la temperatura de los bordes en caso de que no 
* se den los cuatro bordes
* @tparam T
* @param SizeVecBordes
* @param VecTemperaturas
* @param VecBordesValues
<<<<<<< HEAD
*/

=======
// */
// template <typename T>
// void obtainVecBordesValues(const int SizeVecBordes,
//                            const std::vector<T> &VecTemperaturas,
//                            std::vector<T> &VecBordesValues)
// {

//     if (VecTemperaturas.size() >= 3)
//     {
//         spline<T>(SizeVecBordes, VecTemperaturas, VecBordesValues);
//         std::cout << "using splines \n";
//     }

//     else if (VecTemperaturas.size() == 2)
//     {
//         anpi::linspace( VecTemperaturas, VecBordesValues);
//         std::cout << "lineal increment \n";
//     }

//     else if (VecTemperaturas.size() == 1)
//     {
//         for (int i = 0; i < SizeVecBordes; i++)
//         {
//             VecBordesValues.push_back(VecTemperaturas[0]);
//         }
//         std::cout << "constant values \n";
//     }

//     else
//     {
//         std::cout << "border isolated \n";
//     }
// }
>>>>>>> bb9bba5910d7b98728593a550b0bbbd77703c707
} // namespace anpi

#endif