/**
 * @file Liebman.hpp
 * @author Jorge Agüero
 * @brief Liebman Solver to calclute the horizontal and vertical temperature flux on a metal plaque g
 * @version 0.1
 * @date 2018-11-12
 * 
 * @copyright Copyright (c) 2018
 * 
 */

#ifndef LIEBMAN_HPP
#define LIEBMAN_HPP

#include "Matrix.hpp"
#include "MatrixUtils.hpp"
#include <cstdlib>
#include <iostream>
#include <chrono>

namespace anpi
{

//struct that holds the temperatures at the edges of the plaque
/**
 * @brief truct that holds the temperatures at the edges of the plaque
 * 
 */
struct Edge
{
    std::vector<double> temperatures;
    bool isolated, gradient;

    Edge()
    {
        gradient = false;
        isolated = true;
    }
};

/**
 * 
 * @brief CLass that solves and contains solution for the heat plaque
 * 
 * 
 */
class LiebmnanSolver
{

    const int MAX_ITER = 100;
    const int MAX_MATRIX_GROWTH = 10000;
    const double ERROR = 0.001;
    const double LAMBDA = 1.37;

    Edge top, bottom, right, left;
    double error;
    double lambda;
    int height, width;
    bool parallelized = true;
    bool simpleLiebman = false;

  public:
    Matrix<double> tempsMatrix;

    /**
     * @brief Empty constructor for a new Liebmnan Solver object
     * 
     */
    LiebmnanSolver();

    /**
 * @brief Constructs a solver given the different edge temperatures Uses the top and right side sizes to know
 *  what is the size of the matrix to create
 * 
 */
    LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp) //, Matrix<double> &temps)
    {
        top = toptemp;
        bottom = botttemp;
        right = rightemp;
        left = leftemp;
        height = right.temperatures.size();
        width = top.temperatures.size();
        // tempsMatrix = temps;
        error = ERROR;
        lambda = LAMBDA;
    }
    /**
 * @brief Constructs a solver given the different edge temperatures, and a final size of the matrix 
 * 
 */
    LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp, int vsize, int hsize) //, Matrix<double> &temps)
    {
        top = toptemp;
        bottom = botttemp;
        right = rightemp;
        left = leftemp;
        height = vsize;
        width = hsize;
        // tempsMatrix = temps;
        error = ERROR;
        lambda = LAMBDA;
    }

    /**
     * @brief Constructs a solver given the different edge temperatures, a final size of the matrix, and the error desired to use for calculation
     * 
     */
    LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp, int vsize, int hsize, double inerror) //, Matrix<double> &temps)
    {
        top = toptemp;
        bottom = botttemp;
        right = rightemp;
        left = leftemp;
        height = vsize;
        width = hsize;
        // tempsMatrix = temps;
        error = inerror;
        lambda = LAMBDA;
    }

    /**
     * @brief Constructs a solver given the different edge temperatures, a final size of the matrix, the error to use for calculation and the lambda used for relaxation (between 1 and 2)
     * 
     */
    LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp, int vsize, int hsize, double inerror, double inlambda) //, Matrix<double> &temps)
    {
        top = toptemp;
        bottom = botttemp;
        right = rightemp;
        left = leftemp;
        height = vsize;
        width = hsize;
        // tempsMatrix = temps;
        error = inerror;
        lambda = inlambda;
    }

    /**
     * @brief Constructs a solver given the different edge temperatures, a final size of the matrix, the error to use for calculation, the lambda used for relaxation (between 1 and 2), and wether the calculations should be parallelized
     * 
     */
    LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp, int vsize, int hsize, double inerror, double inlambda, bool parallel) //, Matrix<double> &temps)
    {
        top = toptemp;
        bottom = botttemp;
        right = rightemp;
        left = leftemp;
        height = vsize;
        width = hsize;
        // tempsMatrix = temps;
        error = inerror;
        lambda = inlambda;
        parallelized = parallel;
    }

    /**
     * @brief Constructs a solver given the different edge temperatures, a final size of the matrix, the error to use for calculation, the lambda used for relaxation (between 1 and 2), and wether the calculations should be parallelized
     * 
     */
    LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp, int vsize, int hsize, double inerror, double inlambda, bool parallel, bool simpleLieb) //, Matrix<double> &temps)
    {
        top = toptemp;
        bottom = botttemp;
        right = rightemp;
        left = leftemp;
        height = vsize;
        width = hsize;
        // tempsMatrix = temps;
        error = inerror;
        lambda = inlambda;
        parallelized = parallel;
        simpleLiebman = simpleLieb;
    }

    /**
 * @brief Calculates the Liebman method growing the matrix in size in order to reduce  the numbe of iterations required for big matrices
 * 
 * @return double The time in seconds that it took to complete
 */
    double lieb()
    {
        if (simpleLiebman)
        {
            Matrix<double> currMatrix(height, width);

            double timeSec;
            auto t_start = std::chrono::high_resolution_clock::now();
            if (parallelized)
                std::cout << "\nIterations used: " << liebmanOMP(currMatrix, top, bottom, right, left) << std::endl;
            else
                std::cout << "\nIterations used: " << liebman(currMatrix, top, bottom, right, left) << std::endl;

            auto t_end = std::chrono::high_resolution_clock::now();
            timeSec = ((std::chrono::duration<double, std::milli>(t_end - t_start).count()) / 1000);

            tempsMatrix = currMatrix;
            return timeSec;
        }

        int finalRows, finalCols, currRows, currCols,
            totalMatrices, iterPerMatrix;
        finalRows = height;
        finalCols = width;
        //Edge temperatures
        Edge currTop, currBot, currLeft, currRight;

        double timeSec;
        auto t_start = std::chrono::high_resolution_clock::now();
        //start with a size 1 matrix
        Matrix<double> currMatrix(1, 1);
        currCols = 1;
        currRows = 1;
        totalMatrices = 1;

        while (true)
        {

            //creates temperature vector the size of the matrix we are calculating
            currTop = averageEdgeTemp(top, currCols);
            currBot = averageEdgeTemp(bottom, currCols);
            currLeft = averageEdgeTemp(left, currRows);
            currRight = averageEdgeTemp(right, currRows);

            //calculate liebman for current matrix

            if (parallelized)
                iterPerMatrix = liebmanOMP(currMatrix, currTop, currBot, currRight, currLeft);
            else
                iterPerMatrix = liebman(currMatrix, currTop, currBot, currRight, currLeft);

            //print statistics
            std::cout << "\nMatrix number: " << totalMatrices << "  |  size:  " << currRows << "x" << currCols << "  |  used iterations: " << iterPerMatrix << std::endl;

            //calculated liebman for the final matrix
            if (currCols == finalCols && currRows == finalRows)
            {
                tempsMatrix = currMatrix;
                break;
            }

            //else we grow the temperature matrix and calculate liebman again

            // if we are able to duplicate the whole matrix
            if (currCols * 2 <= finalCols && currRows * 2 <= finalRows)
            {
                duplicateTempMatrix(currMatrix);
                currCols = currCols * 2;
                currRows = currRows * 2;
            }
            //if we are able to duplicate width
            else if (currCols * 2 <= finalCols)
            {
                duplicateWidthTempMatrix(currMatrix);
                currCols = currCols * 2;
            }
            //if we are able to duplicate height
            else if (currRows * 2 <= finalRows)
            {
                duplicateHeightTempMatrix(currMatrix);
                currRows = currRows * 2;
            }
            else
            {
                // if the size difference is bigger than a maximun widht, then grow in a step of 64
                // this value was obtained empirically by testing and seeing how the matrix behaved
                //this is because of the way we expand the matrix, as an average cross
                if (currCols + MAX_MATRIX_GROWTH < finalCols && currRows + MAX_MATRIX_GROWTH < finalRows)
                {
                    /* code */
                    fillBiggerTempMatrix(currMatrix, currRows + MAX_MATRIX_GROWTH, currCols + MAX_MATRIX_GROWTH);
                    currCols = currCols + MAX_MATRIX_GROWTH;
                    currRows = currRows + MAX_MATRIX_GROWTH;
                }
                else if (currCols + MAX_MATRIX_GROWTH < finalCols)
                {
                    fillBiggerTempMatrix(currMatrix, currRows, currCols + MAX_MATRIX_GROWTH);
                    currCols = currCols + MAX_MATRIX_GROWTH;
                    // currRows = currRows + MAX_MATRIX_GROWTH;
                }
                else if (currRows + MAX_MATRIX_GROWTH < finalRows)
                {
                    fillBiggerTempMatrix(currMatrix, currRows + MAX_MATRIX_GROWTH, currCols);
                    // currCols = currCols + MAX_MATRIX_GROWTH;
                    currRows = currRows + MAX_MATRIX_GROWTH;
                }
                else
                {
                    /* code */
                    fillBiggerTempMatrix(currMatrix, finalRows, finalCols);
                    currCols = finalCols;
                    currRows = finalRows;
                }
            }

            ++totalMatrices;
            // currMatrix.DumpToFile("currMatrix.txt");
        } //end of while

        auto t_end = std::chrono::high_resolution_clock::now();
        timeSec = ((std::chrono::duration<double, std::milli>(t_end - t_start).count()) / 1000);

        std::cout << "\nMatrices created: " << totalMatrices << std::endl;

        //initial step Matrix size 1
        return timeSec;
    }

    /**
 * @brief Calculates the temperature flow in a plaque using the Liebman algorithm given an initial temperature matrix and the temperature at the edges
 * 
 * @param top The top temperatures, given as an Edge vector of the same size as bottom and the matrix width
 * @param bottom The bottom temperatures, given as an Edge vector of the same size as top and the matrix width
 * @param right The right temperatures, given as an Edge vector of the same size as left and the matrix height
 * @param left The left temperatures, given as an Edge vector of the same size as right and the matrix height
 * @return int The amount of iterations it took to complete calculations
 */
    int liebman(Matrix<double> &temperatureMatrix, Edge top, Edge bottom, Edge right, Edge left)
    {
        int iterations = 0;
        int rows = temperatureMatrix.rows();
        int cols = temperatureMatrix.cols();

        //if Matrix of size 1
        if (rows == 1 && cols == 1)
        {
            bool tbIsolated = false;
            bool lrIsolated = false;
            double tt, bt, lt, rt;

            //if one of the sides is isolated
            if ((top.isolated || bottom.isolated || right.isolated || left.isolated))
            {
                //top isolated
                if (top.isolated)
                {
                    if (bottom.isolated)
                    {
                        // tt = 0;
                        // bt = 0;
                        tbIsolated = true;
                    }
                    else
                    {
                        tt = bottom.temperatures[0];
                        bt = bottom.temperatures[0];
                    }
                }
                //bottom isolated
                else if (bottom.isolated)
                {

                    tt = top.temperatures[0];
                    bt = top.temperatures[0];
                }
                else
                {
                    tt = top.temperatures[0];
                    bt = bottom.temperatures[0];
                }

                //right isolated
                if (right.isolated)
                {
                    if (left.isolated)
                    {
                        lrIsolated = true;
                    }
                    else
                    {
                        rt = left.temperatures[0];
                        lt = left.temperatures[0];
                    }
                }
                //left isolated
                else if (left.isolated)
                {
                    lt = right.temperatures[0];
                    rt = right.temperatures[0];
                }
                else
                {
                    rt = right.temperatures[0];
                    lt = left.temperatures[0];
                }
            }
            else
            {
                tt = bottom.temperatures[0];
                bt = bottom.temperatures[0];
                rt = right.temperatures[0];
                lt = left.temperatures[0];
            }

            //if the plaque is completely isolated
            if (lrIsolated && tbIsolated)
            {
                temperatureMatrix[0][0] = 0;
            }
            //both sides isolated
            else if (lrIsolated)
            {
                temperatureMatrix[0][0] = (tt + bt) / 2;
            }
            //top and bottom isolated
            else if (tbIsolated)
            {
                temperatureMatrix[0][0] = (rt + lt) / 2;
            }
            //calculate normally
            else
            {
                temperatureMatrix[0][0] = (rt + lt + tt + bt) / 4;
            }
            //no need to iterate, will always get the same result
            return 1;
        } //end of matrix size 1

        double bigerror = 0;

        double currTemp, calcTemp, calcError;
        do
        {
            //holds the value of the biggest difference (error) calculated
            bigerror = 0;

            //iterate through rows
            for (int i = 0; i < rows; ++i)
            {
                //iterate through columns
                for (int j = 0; j < cols; ++j)
                {

                    currTemp = temperatureMatrix[i][j];

                    //////////////////////////////////top check///////////////////////////////////////////
                    //if on the top
                    if (i == 0)
                    {

                        //if on the left corner
                        if (j == 0)
                        {
                            //if top is isolated
                            if (top.isolated)
                            {

                                //if left is also isolated
                                if (left.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i + 1][j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                }
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                }
                            }
                            //top not isolated
                            else
                            {
                                //if left is isolated
                                if (left.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                }
                                //neither is isolated
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                }
                            }
                        }

                        //if on the right corner
                        else if (j == cols - 1)
                        {
                            //if top is isolated
                            if (top.isolated)
                            {

                                //if right is also isolated
                                if (right.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i + 1][j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                }
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j - 1] + right.temperatures[i]) / 4;
                                }
                            }
                            //top not isolated
                            else
                            {
                                //if right is isolated
                                if (right.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                }
                                //neither is isolated
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j - 1] + right.temperatures[i]) / 4;
                                }
                            }
                        }
                        //other parts of the top
                        else
                        {
                            if (top.isolated)
                            {
                                //calculate temperature
                                calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                            }
                            else
                            {
                                //calculate temperature
                                calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                            }
                        }
                    }

                    //////////////////////////////////bottom check///////////////////////////////////////////
                    //if at the bottom
                    else if (i == rows - 1)
                    {
                        //if on the left corner
                        if (j == 0)
                        {
                            //if bottom is isolated
                            if (bottom.isolated)
                            {
                                //if left is also isolated
                                if (left.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i - 1][j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                }
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                }
                            }
                            //bottom not isolated
                            else
                            {
                                //if left is isolated
                                if (left.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                }
                                //neither is isolated
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                }
                            }
                        }
                        //if on the right corner
                        else if (j == cols - 1)
                        {
                            //if bottom is isolated
                            if (bottom.isolated)
                            {
                                //if right is also isolated
                                if (right.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i - 1][j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                }
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + right.temperatures[i]) / 4;
                                }
                            }
                            //bottom not isolated
                            else
                            {
                                //if right is isolated
                                if (right.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                }
                                //neither is isolated
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j + 1] + right.temperatures[i]) / 4;
                                }
                            }
                        }
                        //other parts of the bottom
                        else
                        {
                            if (bottom.isolated)
                            {
                                //calculate temperature
                                calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                            }
                            else
                            {
                                //calculate temperature
                                calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                            }
                        }
                    }

                    //corners were checked when calculating top and bottom

                    //////////////////////////////////left check///////////////////////////////////////////
                    //if on the left side non corner
                    else if (j == 0)
                    {
                        if (left.isolated)
                        {
                            //calculate temperature
                            calcTemp = (2 * temperatureMatrix[i][j + 1] + temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j]) / 4;
                        }
                        else
                        {
                            //calculate temperature
                            calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + left.temperatures[i] + temperatureMatrix[i][j + 1]) / 4;
                        }
                    }

                    //////////////////////////////////right check///////////////////////////////////////////
                    //if on the right side non corner
                    else if (j == cols - 1)
                    {
                        if (right.isolated)
                        {
                            //calculate temperature
                            calcTemp = (2 * temperatureMatrix[i][j - 1] + temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j]) / 4;
                        }
                        else
                        {
                            //calculate temperature
                            calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + right.temperatures[i] + temperatureMatrix[i][j - 1]) / 4;
                        }
                    }

                    //////////////////////////////////other cells///////////////////////////////////////////
                    //every other cell
                    else
                    {
                        //calculate temperature
                        calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + temperatureMatrix[i][j - 1]) / 4;
                    }

                    //apply relaxation
                    calcTemp = lambda * calcTemp + (1 - lambda) * currTemp;
                    //set the new temperture
                    temperatureMatrix[i][j] = calcTemp;

                    //calculate the error
                    // calcError = std::abs(std::abs(calcTemp) - std::abs(currTemp)) / std::abs(calcTemp) * 100;
                    calcError = std::abs(calcTemp - currTemp);

                    //check if the new calculated error is the biggest to compare
                    if (calcError > bigerror)
                    {
                        bigerror = calcError;
                    }

                } //end for column
            }     //end for row
            //////////////////////////////////end matrix iteration///////////////////////////////////////////

            ++iterations;

            //calculate error of the biggest temperature

        } while (bigerror > error);

        return iterations;
    }

    /**
 * @brief Calculates the temperature flow in a plaque using the Liebman algorithm, and using openMP directives to paralelize the computation
 * 
 * @param temperatureMatrix The intial temperature matrix used to calculate the new temperatures, the values of the matrix are updated with the new values
 * @param top The top temperatures, given as an Edge vector of the same size as bottom and the matrix width
 * @param bottom The bottom temperatures, given as an Edge vector of the same size as top and the matrix width
 * @param right The right temperatures, given as an Edge vector of the same size as left and the matrix height
 * @param left The left temperatures, given as an Edge vector of the same size as right and the matrix height
 * @return int The amount of iterations it took to complete calculations
 */
    int liebmanOMP(Matrix<double> &temperatureMatrix, Edge top, Edge bottom, Edge right, Edge left)
    {
        int iterations = 0;
        int rows = temperatureMatrix.rows();
        int cols = temperatureMatrix.cols();

        //if Matrix of size 1
        if (rows == 1 && cols == 1)
        {
            bool tbIsolated = false;
            bool lrIsolated = false;
            double tt, bt, lt, rt;

            //if one of the sides is isolated
            if ((top.isolated || bottom.isolated || right.isolated || left.isolated))
            {
                //top isolated
                if (top.isolated)
                {
                    if (bottom.isolated)
                    {
                        // tt = 0;
                        // bt = 0;
                        tbIsolated = true;
                    }
                    else
                    {
                        tt = bottom.temperatures[0];
                        bt = bottom.temperatures[0];
                    }
                }
                //bottom isolated
                else if (bottom.isolated)
                {

                    tt = top.temperatures[0];
                    bt = top.temperatures[0];
                }
                else
                {
                    tt = top.temperatures[0];
                    bt = bottom.temperatures[0];
                }

                //right isolated
                if (right.isolated)
                {
                    if (left.isolated)
                    {
                        lrIsolated = true;
                    }
                    else
                    {
                        rt = left.temperatures[0];
                        lt = left.temperatures[0];
                    }
                }
                //left isolated
                else if (left.isolated)
                {
                    lt = right.temperatures[0];
                    rt = right.temperatures[0];
                }
                else
                {
                    rt = right.temperatures[0];
                    lt = left.temperatures[0];
                }
            }
            else
            {
                tt = bottom.temperatures[0];
                bt = bottom.temperatures[0];
                rt = right.temperatures[0];
                lt = left.temperatures[0];
            }

            //if the plaque is completely isolated
            if (lrIsolated && tbIsolated)
            {
                temperatureMatrix[0][0] = 0;
            }
            //both sides isolated
            else if (lrIsolated)
            {
                temperatureMatrix[0][0] = (tt + bt) / 2;
            }
            //top and bottom isolated
            else if (tbIsolated)
            {
                temperatureMatrix[0][0] = (rt + lt) / 2;
            }
            //calculate normally
            else
            {
                temperatureMatrix[0][0] = (rt + lt + tt + bt) / 4;
            }
            //no need to iterate, will always get the same result
            return 1;
        } //end of matrix size 1
          //////////////////////////////////////////////////////////////////////////////////////////

        double parallelError, localError, bigError; //, oldbigError;

        double currTemp, calcTemp; //, calcError;
        do
        {

            parallelError = 0;

#pragma omp parallel private(calcTemp, currTemp, localError, bigError)
            {
                localError = 0;
                bigError = 0;
#pragma omp for collapse(2)
                //iterate through rows
                for (int i = 0; i < rows; ++i)
                {
                    //iterate through columns
                    for (int j = 0; j < cols; ++j)
                    {

                        currTemp = temperatureMatrix[i][j];

                        //////////////////////////////////top check///////////////////////////////////////////
                        //if on the top
                        if (i == 0)
                        {

                            //if on the left corner
                            if (j == 0)
                            {
                                //if top is isolated
                                if (top.isolated)
                                {

                                    //if left is also isolated
                                    if (left.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i + 1][j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                    }
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                    }
                                }
                                //top not isolated
                                else
                                {
                                    //if left is isolated
                                    if (left.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                    }
                                    //neither is isolated
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                    }
                                }
                            }

                            //if on the right corner
                            else if (j == cols - 1)
                            {
                                //if top is isolated
                                if (top.isolated)
                                {

                                    //if right is also isolated
                                    if (right.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i + 1][j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                    }
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j - 1] + right.temperatures[i]) / 4;
                                    }
                                }
                                //top not isolated
                                else
                                {
                                    //if right is isolated
                                    if (right.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                    }
                                    //neither is isolated
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j - 1] + right.temperatures[i]) / 4;
                                    }
                                }
                            }
                            //other parts of the top
                            else
                            {
                                if (top.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                                }
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                                }
                            }
                        }

                        //////////////////////////////////bottom check///////////////////////////////////////////
                        //if at the bottom
                        else if (i == rows - 1)
                        {
                            //if on the left corner
                            if (j == 0)
                            {
                                //if bottom is isolated
                                if (bottom.isolated)
                                {
                                    //if left is also isolated
                                    if (left.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i - 1][j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                    }
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                    }
                                }
                                //bottom not isolated
                                else
                                {
                                    //if left is isolated
                                    if (left.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + 2 * temperatureMatrix[i][j + 1]) / 4;
                                    }
                                    //neither is isolated
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
                                    }
                                }
                            }
                            //if on the right corner
                            else if (j == cols - 1)
                            {
                                //if bottom is isolated
                                if (bottom.isolated)
                                {
                                    //if right is also isolated
                                    if (right.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i - 1][j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                    }
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + right.temperatures[i]) / 4;
                                    }
                                }
                                //bottom not isolated
                                else
                                {
                                    //if right is isolated
                                    if (right.isolated)
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + 2 * temperatureMatrix[i][j - 1]) / 4;
                                    }
                                    //neither is isolated
                                    else
                                    {
                                        //calculate temperature
                                        calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j + 1] + right.temperatures[i]) / 4;
                                    }
                                }
                            }
                            //other parts of the bottom
                            else
                            {
                                if (bottom.isolated)
                                {
                                    //calculate temperature
                                    calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                                }
                                else
                                {
                                    //calculate temperature
                                    calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
                                }
                            }
                        }

                        //corners were checked when calculating top and bottom

                        //////////////////////////////////left check///////////////////////////////////////////
                        //if on the left side non corner
                        else if (j == 0)
                        {
                            if (left.isolated)
                            {
                                //calculate temperature
                                calcTemp = (2 * temperatureMatrix[i][j + 1] + temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j]) / 4;
                            }
                            else
                            {
                                //calculate temperature
                                calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + left.temperatures[i] + temperatureMatrix[i][j + 1]) / 4;
                            }
                        }

                        //////////////////////////////////right check///////////////////////////////////////////
                        //if on the right side non corner
                        else if (j == cols - 1)
                        {
                            if (right.isolated)
                            {
                                //calculate temperature
                                calcTemp = (2 * temperatureMatrix[i][j - 1] + temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j]) / 4;
                            }
                            else
                            {
                                //calculate temperature
                                calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + right.temperatures[i] + temperatureMatrix[i][j - 1]) / 4;
                            }
                        }

                        //////////////////////////////////other cells///////////////////////////////////////////
                        //every other cell
                        else
                        {
                            //calculate temperature
                            calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + temperatureMatrix[i][j - 1]) / 4;
                        }

                        //apply relaxation
                        calcTemp = lambda * calcTemp + (1 - lambda) * currTemp;
                        //set the new temperture
                        temperatureMatrix[i][j] = calcTemp;

                        localError = std::abs(calcTemp - currTemp);

                        //check if the new calculated temp is the biggest to compare for error
                        if (localError > bigError)
                        {
                            bigError = localError;
                        }

                    } //end for column
                }     //end for row
                      //////////////////////////////////end matrix iteration///////////////////////////////////////////

//runs for every thread, lloks for the highest error found by each thread
#pragma omp critical
                {
                    if (parallelError < bigError)
                    {
                        parallelError = bigError;
                    }
                }
            }
            //end PARALELL section

            ++iterations;

            //calculate error of the biggest temperature

            // calcError = (std::abs(parallelError)) * 100;

        } while (parallelError > error);

        return iterations;
    }

    /**
 * @brief Duplicates The size of the matrix and copies each value into 4 cells of the new matrix. Uses  openMP to paralelize the duplication
 * 
 * @param tempMatrix The matrix to duplicate
 */
    void duplicateTempMatrix(Matrix<double> &tempMatrix)
    {
        int scols = tempMatrix.cols();
        int srows = tempMatrix.rows();
        int bcols = scols * 2;
        int brows = srows * 2;
        int bigI, bigJ;
        Matrix<double> big(brows, bcols);
        bigI = 0;
        bigJ = 0;

#pragma omp parallel for collapse(2) private(bigI, bigJ)
        for (int i = 0; i < srows; ++i)
            for (int j = 0; j < scols; ++j)
            {
                bigI = i * 2;
                bigJ = j * 2;
                big[bigI][bigJ] = tempMatrix[i][j];
                big[bigI + 1][bigJ] = tempMatrix[i][j];
                big[bigI][bigJ + 1] = tempMatrix[i][j];
                big[bigI + 1][bigJ + 1] = tempMatrix[i][j];
            }

        tempMatrix = big;
    }

    /**
     * @brief Duplicates The width of a given matrix copies each value into 2 other horizontal values of the new matrix. Uses  openMP to paralelize the duplication
     * 
     * @param tempMatrix The matrix to duplicate
     */
    void duplicateWidthTempMatrix(Matrix<double> &tempMatrix)
    {
        int scols = tempMatrix.cols();
        int srows = tempMatrix.rows();
        int bcols = scols * 2;
        int brows = srows;
        Matrix<double> big(brows, bcols);

        int bigI, bigJ;
#pragma omp parallel for private(bigI, bigJ) collapse(2)
        for (int i = 0; i < srows; ++i)
            for (int j = 0; j < scols; ++j)
            {
                bigI = i;
                bigJ = j * 2;
                big[bigI][bigJ] = tempMatrix[i][j];
                big[bigI][bigJ + 1] = tempMatrix[i][j];
            }

        tempMatrix = big;
    }

    /**
     * @brief Duplicates the height of a given matrix copies each value into 2 other vertical values of the new matrix. Uses  openMP to paralelize the duplication
     * 
     * @param tempMatrix 
     */
    void duplicateHeightTempMatrix(Matrix<double> &tempMatrix)
    {
        int scols = tempMatrix.cols();
        int srows = tempMatrix.rows();
        int bcols = scols;
        int brows = srows * 2;

        Matrix<double> big(brows, bcols);
        int bigI, bigJ;

        // #pragma omp parallel private(i, j, bigI, bigJ)
        {
#pragma omp for collapse(2) private(bigI, bigJ)
            for (int i = 0; i < srows; ++i)
                for (int j = 0; j < scols; ++j)
                {
                    bigI = i * 2;
                    bigJ = j;
                    big[bigI][bigJ] = tempMatrix[i][j];
                    big[bigI + 1][bigJ] = tempMatrix[i][j];
                }
        }

        tempMatrix = big;
    }

    /**
     * @brief Fills a matrix up to the size specified. It grows the matrix by adding the average of the values of the middle temperatures to a centered cross section of the expanded matrix
     * 
     * @param small The matrix used to fill the new matrix, and also were this new matrix is returned
     * @param rows The row size of the new bigger matrix
     * @param cols The column size of the new bigger matrix
     */
    void
    fillBiggerTempMatrix(Matrix<double> &small, int rows, int cols)
    {
        int scols = small.cols();
        int srows = small.rows();
        int colDiff = cols - scols;
        int rowDiff = rows - srows;

        if (rows < srows || cols < scols)
        {
            throw anpi::Exception("new Matrix has to have bigger indices");
        }

        Matrix<double> big(rows, cols);

        //growing from the inside like a cross, using the average of the side cells as the new temp

        // #pragma omp parallel for
        for (int is = 0, ib = 0; is < srows; ++is, ++ib)
        {
            for (int js = 0, jb = 0; js < scols; ++js, ++jb)
            {
                //if on the middle horizontal band
                if (is == srows / 2)
                {

                    //if on the middle rectangle
                    if (js == scols / 2)
                    {
                        double temp = (small[is - 1][js - 1] + small[is - 1][js] + small[is][js - 1] + small[is][js]) / 4;
                        double rightSidetemp = (small[is][js] + small[is][js - 1]) / 2;
                        double downSidetemp = (small[is][js] + small[is - 1][js]) / 2;

                        //for each new item
                        for (int k = 0; k <= rowDiff; ++k)
                        {
                            for (int p = 0; p <= colDiff; p++)
                            {
                                big[ib + k][jb + p] = temp;

                                //add corner temp
                                if (p == colDiff && k == rowDiff)
                                {
                                    big[ib + k][jb + p] = small[is][js];
                                }
                                //add the right side temp
                                else if (p == colDiff)
                                {
                                    big[ib + k][jb + p] = rightSidetemp;
                                }
                                //add bottom temps
                                else if (k == rowDiff)
                                {
                                    big[ib + k][jb + p] = downSidetemp;
                                }

                                // ++jb;
                            }
                        }
                        jb += colDiff;
                    }
                    // if not on the rectangle, grow vertically
                    else
                    {
                        //temp is the average between the middle temps
                        double temp = (small[is][js] + small[is - 1][js]) / 2;
                        for (int k = 0; k < rowDiff; ++k)
                        {
                            big[ib + k][jb] = temp;
                        }
                        big[ib + rowDiff][jb] = small[is][js];
                    }
                }

                //if on the middle vertical band (horizontal growth)
                else if (js == scols / 2)
                {
                    double temp = (small[is][js] + small[is][js - 1]) / 2;

                    for (int p = 0; p < colDiff; ++p)
                    {
                        big[ib][jb + p] = temp;
                    }

                    jb += colDiff;

                    big[ib][jb] = small[is][js];
                }
                else
                {
                    big[ib][jb] = small[is][js];
                }
            }
            //if on the middle horizontal band
            if (is == srows / 2)
                ib += rowDiff;
        }

        // for (int i = 1; i <= srows; ++i)
        //     for (int j = 1; j <= scols; ++j)
        //     {
        //         big[i][j] = small[i - 1][j - 1];
        //     }
        // big.fill(small);
        // big.DumpToFile("UltimaCrecidaenUno.txt");
        small = big;

        // big.fill(small);

        // for (int i = 0; i < srows; ++i)
        //     for (int j = scols; j < cols; ++j)
        //     {
        //         big[i][j] = small[i][j - colDiff];
        //     }
        // for (int i = srows; i < rows; ++i)
        //     for (int j = 0; j < cols; ++j)
        //     {
        //         if (j < scols)
        //         {
        //             big[i][j] = small[i - rowDiff][j];
        //         }
        //         else
        //         {
        //             big[i][j] = small[i - rowDiff][j - colDiff];
        //         }
        //     }

        small = big;
    }

    /**
 * @brief Given an Edge vector of temps, returns another smaller Edge temperature of the specified size, by averaging the temperature values of the given Edge into equal sections to form the new one
 * 
 * @param temps The Edge temperature vector used to divide and calculate the temperatures of the new one
 * @param size The desired temperature size Edge
 * @return Edge 
 */
    Edge averageEdgeTemp(Edge &temps, int size)
    {
        //if it's isolated no need to work
        if (temps.isolated)
            return temps;

        int tempSize = temps.temperatures.size();

        //if the size is the same, just return the temps
        if (tempSize == size)
        {
            return temps;
        }

        Edge avgEdge;
        avgEdge.isolated = false;
        std::vector<double> avgTemps;

        //if it's all the same temperature
        if (!temps.gradient)
        {
            avgTemps.resize(size, temps.temperatures[0]);
            avgEdge.temperatures = avgTemps;
            return avgEdge;
        }
        else
        {
            avgEdge.gradient = true;
        }

        int chunkSize = tempSize / size;
        int overflow = tempSize % size;
        int end, j;
        double average;
#pragma omp parallel for private(average)
        for (int i = 0; i <= size; ++i)
        {
            average = 0;
            //if on the last chunk we also average the leftover temperatures
            if (i == size - 1)
            {
                j = i * chunkSize;
                end = j + chunkSize + overflow;
                for (; j < end; ++j)
                {
                    average += temps.temperatures[j];
                }
                average = average / (chunkSize + overflow);
            }
            else
            {
                j = i * chunkSize;
                end = j + chunkSize;
                for (; j < end; ++j)
                {
                    average += temps.temperatures[j];
                }
                average = average / chunkSize;
            }

            avgTemps.push_back(average);
        }

        avgEdge.temperatures = avgTemps;
        avgEdge.isolated = false;
        return avgEdge;
    }
}; // namespace anpi

} // namespace anpi

#endif // LIEBMAN_HPP