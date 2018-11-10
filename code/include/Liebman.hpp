#include "Matrix.hpp"

namespace anpi
{

//struct that holds the temperatures at the edges of the plaque
struct Edge
{
    std::vector<double> temperatures;
    bool isolated;

    Edge()
    {
        isolated = false;
    }
};

class LiebmnanSolver
{

    Edge top, bottom, right, left;
    double error;
    double lambda;

  public:
    Matrix<double> temperatureMatrix;

    /**
     * @brief Empty constructor for a new Liebmnan Solver object
     * 
     */
    LiebmnanSolver();

    /**
 * @brief Construct a new Liebmnan Solver object
 * 
 * @param toptemp vector of temperatures on the top of the plaque
 * @param botttemp vector of temperatures on the bottom of the plaque
 * @param rightemp vector of temperatures on the right side of the plaque
 * @param leftemp vector of temperatures on the left side of the plaque
 * @param hsize the number of horizontal nodes to calculate, must be the same size
 *                  as toptemp an bottemp if they are not isolated.
 * @param vsize the number of vertical nodes to calculate, must be the same size
 *                  as lefttemp and righttemp if they are not isolated.
 */
    LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp, Matrix<double> &temps)
    {
        top = toptemp;
        bottom = botttemp;
        right = rightemp;
        left = leftemp;
        temperatureMatrix = temps;
        error = 0.0001;
        lambda = 1.5;
    }

    int liebman()
    {
        int iterations = 0;
        int rows = temperatureMatrix.rows();
        int cols = temperatureMatrix.cols();

        double bigTemp, bigTempOld;

        double currTemp, calcTemp, calcError;
        do
        {
            //iterate through rows
            for (int i = 0; i < rows; ++i)
            {
                //iterate through columns
                for (int j = 0; j < cols; ++j)
                {
                    bigTemp = 0;
                    bigTempOld = 0;
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

                    //aplyrelaxation
                    calcTemp = lambda * calcTemp + (1 - lambda) * currTemp;
                    //set the new temperture
                    temperatureMatrix[i][j] = calcTemp;

                    //check if the new calculated temp is the biggest to compare for error
                    if (calcTemp > bigTemp)
                    {
                        bigTemp = calcTemp;
                        bigTempOld = currTemp;
                    }

                } //end for column
            }     //end for row
            //////////////////////////////////end matrix iteration///////////////////////////////////////////

            ++iterations;

            //calculate error of the biggest temperature
            calcError = (abs(bigTemp - bigTempOld) / bigTemp) * 100;

        } while (calcError > error);

        return iterations;
    }
};

} // namespace anpi