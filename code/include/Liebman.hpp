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
 * @brief CLass that solves and contains solution
 * 
 * 
 */
class LiebmnanSolver
{
    const int MAX_ITER = 100;
    const int MAX_MATRIX_GROWTH = 2;
    const double ERROR = 0.001;
    const double LAMBDA = 1.37;

    Edge top, bottom, right, left;
    double error;
    double lambda;
    int height, width;

  public:
    Matrix<double> tempsMatrix;

    /**
     * @brief Empty constructor for a new Liebmnan Solver object
     * 
     */
    LiebmnanSolver();

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

    int lieb()
    {
        int finalRows, finalCols, currRows, currCols,
            totalMatrices, iterPerMatrix;
        finalRows = height;
        finalCols = width;

        //start with a size 1 matrix
        Matrix<double> currMatrix(1, 1);
        currCols = 1;
        currRows = 1;
        totalMatrices = 1;

        Edge currTop, currBot, currLeft, currRight;

        while (true)
        {
            // currCols = currMatrix.cols();
            // currRows = currMatrix.rows();

            //creates temperature vector the size of the matrix we are calculating
            currTop = averageEdgeTemp(top, currCols);
            currBot = averageEdgeTemp(bottom, currCols);
            currLeft = averageEdgeTemp(left, currRows);
            currRight = averageEdgeTemp(right, currRows);

            //calculata liebman for current matrix
            iterPerMatrix = liebman(currMatrix, currTop, currBot, currRight, currLeft);

            // iterPerMatrix = liebmanOMP(currMatrix, currTop, currBot, currRight, currLeft);

            //print statistics
            std::cout << "\nMatrix number: " << totalMatrices << "  |  size:  " << currRows << "x" << currCols << "  |  used iterations: " << iterPerMatrix << std::endl;

            //calculated liebman for the final matrix
            if (currCols == finalCols && currRows == finalRows)
            {
                tempsMatrix = currMatrix;
                break;
            }

            // if (currCols + 2 > finalCols || currRows + 2 > finalRows)
            // {
            //     growOneMatrix(currMatrix);
            // }
            // else
            // {
            //     growMatrix(currMatrix);
            // }

            // currCols = currMatrix.cols();
            // currRows = currMatrix.rows();

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
                //if the size difference is bigger than 64, then grow in a step of 64
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

                // if (currCols + 2 > finalCols || currRows + 2 > finalRows)
                // {
                //     growOneMatrix(currMatrix);
                //     currCols += 1;
                //     currRows += 1;
                // }
                // else
                // {
                //     growMatrix(currMatrix);
                //     currCols += 2;
                //     currRows += 2;
                // }
            }

            ++totalMatrices;
            currMatrix.DumpToFile("currMatrix.txt");
        }

        std::cout << "\nMatrices created: " << totalMatrices << std::endl;
        //initial step Matrix size 1
        return totalMatrices;
    }

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
            // bigTemp = 0;
            // bigTempOld = 0;
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
                    calcError = std::abs(std::abs(calcTemp) - std::abs(currTemp)) / std::abs(calcTemp) * 100;

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

        double parallelError, localError, bigTemp; //, oldBigTemp;

        double currTemp, calcTemp; //, calcError;
        do
        {
            parallelError = 0;

#pragma omp parallel shared(parallelError) private(calcTemp, currTemp, localError, bigTemp)
            {
                localError = 0;
                bigTemp = 0;
#pragma omp for
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

                        //check if the new calculated temp is the biggest to compare for error
                        if (std::abs(calcTemp) > bigTemp)
                        {
                            bigTemp = std::abs(calcTemp);
                            // oldBigTemp = std::abs(currTemp);
                            localError = (std::abs(calcTemp) - std::abs(currTemp)) / std::abs(calcTemp) * 100;
                            // bigTempOld = std::abs(currTemp);
                        }

                    } //end for column
                }     //end for row
                      //////////////////////////////////end matrix iteration///////////////////////////////////////////

#pragma omp critical
                {
                    if (parallelError < localError)
                    {
                        parallelError = localError;
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

        // fill outer sides of matrix
        // for (int i = 0; i < srows / 2; ++i)
        // {
        //     for (int j = 0; j < scols / 2; ++j)
        //     {
        //         big[bigI][bigJ] = tempMatrix[i][j];
        //         ++bigJ;
        //     }
        //     ++bigI;
        //     bigJ = 0;
        // }
        // bigJ += srows;
        // bigJ = (scols / 2) + scols;

        // for (int i = srows / 2; i < srows; ++i)
        // {
        //     for (int j = scols / 2; j < scols; ++j)
        //     {
        //         big[bigI][bigJ] = tempMatrix[i][j];
        //         ++bigJ;
        //     }
        //     ++bigI;
        //     bigJ = (scols / 2) + scols;
        // }

        // duplicate with matrix at it's center

        // bigI = srows / 2;
        // bigJ = scols / 2;

        // for (int i = 0; i < srows; ++i)
        // {
        //     for (int j = 0; j < scols; ++j)
        //     {
        //         big[bigI][bigJ] = tempMatrix[i][j];
        //         ++bigJ;
        //     }
        //     ++bigI;
        //     bigJ = scols / 2;
        // }
        // // anpi::printMatrix(big);
        // big.DumpToFile("ultimaDuplicada.txt");

        // duplicate by copying
        // int bigI, bigJ;

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

    void growMatrix(Matrix<double> &small)
    {
        int scols = small.cols();
        int srows = small.rows();
        int bcols = scols + 2;
        int brows = srows + 2;

        Matrix<double> big(brows, bcols);

        // for (int i = 1; i <= srows; ++i)
        //     for (int j = 1; j <= scols; ++j)
        //     {
        //         big[i][j] = small[i - 1][j - 1];
        //     }

        //growing from the inside like a cross, using the average of the side cells as the new temp
        // #pragma omp parallel for
        for (int is = 0, ib = 0; is < srows; ++is, ++ib)
        {
            for (int js = 0, jb = 0; js < scols; ++js, ++jb)
            {
                //if on the middle horizontal band
                if (is == srows / 2)
                {
                    // if on the middle square
                    if (js == scols / 2)
                    {

                        // big[ib][jb] = small[is][js];

                        //all central square with same value
                        big[ib][jb] = big[ib][jb + 1] = big[ib + 1][jb] = big[ib + 1][jb + 1] = (small[is][js - 1] + small[is - 1][js] + small[is][js] + small[is - 1][js - 1]) / 4;

                        // big[ib][jb + 1] = (small[is][js] + small[is - 1][js]) / 2;

                        // big[ib + 1][jb] = (small[is][js - 1] + small[is + 1][js]) / 2;
                        // big[ib + 1][jb + 1] = (small[is][js] + small[is + 1][js]) / 2;

                        // set the lower squares
                        big[ib + 2][jb] = big[ib + 2][jb + 1] = (small[is][js] + small[is][js - 1]) / 2;
                        jb += 2;

                        big[ib][jb] = big[ib + 1][jb] = (small[is][js] + small[is - 1][js]) / 2;

                        big[ib + 2][jb] = small[is][js];
                    }

                    else
                    {

                        // big[ib][jb] = small[is][js];
                        big[ib][jb] = big[ib + 1][jb] = (small[is][js] + small[is - 1][js]) / 2;
                        big[ib + 2][jb] = small[is][js];
                        // ib += 2;
                    }
                }

                //if on the middle vertical band
                else if (js == scols / 2)
                {

                    // big[ib][jb] = small[is][js];
                    big[ib][jb] = big[ib][jb + 1] = (small[is][js] + small[is][js - 1]) / 2;
                    jb += 2;

                    big[ib][jb] = small[is][js];
                }
                else
                {
                    big[ib][jb] = small[is][js];
                }
                // big.DumpToFile("big.txt");
            }

            //if on the middle horizontal band
            if (is == srows / 2)
                ib += 2;
        }
        // big.fill(small);
        // small.DumpToFile("SMALL.txt");
        // big.DumpToFile("big.txt");
        small = big;
    }

    void growOneMatrix(Matrix<double> &small)
    {
        int scols = small.cols();
        int srows = small.rows();
        int bcols = scols + 1;
        int brows = srows + 1;

        Matrix<double> big(brows, bcols);

        //growing from the inside like a cross, using the average of the side cells as the new temp
        for (int is = 0, ib = 0; is < srows; ++is, ++ib)
        {
            for (int js = 0, jb = 0; js < scols; ++js, ++jb)
            {
                //if on the middle horizontal band
                if (is == srows / 2)
                {
                    big[ib + 1][jb] = small[is][js];
                    big[ib][jb] = (small[is][js] + small[is - 1][js]) / 2;

                    //if on the middle square
                    if (js == scols / 2)
                    {
                        big[ib][jb] = (small[is - 1][js - 1] + small[is - 1][js] + small[is][js - 1] + small[is][js]) / 4;
                        big[ib + 1][jb] = (small[is][js] + small[is][js - 1]) / 2;
                        big[ib][jb + 1] = (small[is][js] + small[is - 1][js]) / 2;

                        big[ib + 1][jb + 1] = small[is][js];

                        ++jb;
                    }
                }

                //if on the middle vertical band
                else if (js == scols / 2)
                {
                    big[ib][jb] = (small[is][js] + small[is][js - 1]) / 2;

                    ++jb;

                    big[ib][jb] = small[is][js];
                }
                else
                {
                    big[ib][jb] = small[is][js];
                }
            }
            //if on the middle horizontal band
            if (is == srows / 2)
                ++ib;
        }

        // for (int i = 1; i <= srows; ++i)
        //     for (int j = 1; j <= scols; ++j)
        //     {
        //         big[i][j] = small[i - 1][j - 1];
        //     }
        // big.fill(small);
        // big.DumpToFile("UltimaCrecidaenUno.txt");
        small = big;
    }

    /**
     * @brief Fills a matrix that is less than double the size of another matrix used to fill it. It copies 
     * the last rigthmost columns of the smaller matrix  into the rightmost columns of the big matirx, it
     * does similar to the bottom part, but the bottom part of the rigmtost section uses the values set when
     * filling the rightmost section
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

    void fillBiggerTempMatrix(Matrix<double> &small, Matrix<double> &big)
    {
        int scols = small.cols();
        int bcols = big.cols();
        int srows = small.rows();
        int brows = big.rows();

        //matrix size was duplicated each old cell corresponds to 4 cells of the new one
        if (scols * 2 == bcols && srows * 2 == brows)
        {
            int bigI, bigJ;
            for (int i = 0; i < srows; ++i)
                for (int j = 0; j < scols; ++j)
                {
                    bigI = i * 2;
                    bigJ = j * 2;
                    big[bigI][bigJ] = small[i][j];
                    big[bigI + 1][bigJ] = small[i][j];
                    big[bigI][bigJ + 1] = small[i][j];
                    big[bigI + 1][bigJ + 1] = small[i][j];
                }
        }
        else
        {
            //we just fill the bigger matrix with the small one
            big.fill(small);
        }
    }

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