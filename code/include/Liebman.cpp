// #include "Liebman.hpp"

// namespace anpi
// {
// using namespace anpi;

// LiebmnanSolver::LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp) //, Matrix<double> &temps)
// {
//     top = toptemp;
//     bottom = botttemp;
//     right = rightemp;
//     left = leftemp;
//     height = right.temperatures.size();
//     width = top.temperatures.size();
//     // tempsMatrix = temps;
//     error = ERROR;
//     lambda = LAMBDA;
// }
// LiebmnanSolver::LiebmnanSolver(Edge toptemp, Edge botttemp, Edge rightemp, Edge leftemp, int vsize, int hsize)
// {
//     right = rightemp;
//     left = leftemp;
//     height = vsize;
//     width = hsize;
//     // tempsMatrix = temps;
//     error = ERROR;
//     lambda = LAMBDA;
// }

// int LiebmnanSolver::lieb()
// {
//     int finalRows, finalCols, currRows, currCols,
//         totalMatrices, iterPerMatrix;
//     finalRows = height;
//     finalCols = width;

//     //start with a size 1 matrix
//     Matrix<double> currMatrix(1, 1);
//     currCols = 1;
//     currRows = 1;
//     totalMatrices = 1;

//     Edge currTop, currBot, currLeft, currRight;

//     while (true)
//     {
//         // currCols = currMatrix.cols();
//         // currRows = currMatrix.rows();

//         //creates temperature vector the size of the matrix we are calculating
//         currTop = averageEdgeTemp(top, currCols);
//         currBot = averageEdgeTemp(bottom, currCols);
//         currLeft = averageEdgeTemp(left, currRows);
//         currRight = averageEdgeTemp(right, currRows);

//         //calculata liebman for current matrix
//         iterPerMatrix = liebman(currMatrix, currTop, currBot, currRight, currLeft);

//         //print statistics
//         std::cout << "\nMatrix number: " << totalMatrices << "      used iterations: " << iterPerMatrix << std::endl;

//         //calcualted liebman for the final matrix
//         if (currCols == finalCols && currRows == finalRows)
//         {
//             tempsMatrix = currMatrix;
//             break;
//         }

//         //else we grow the temperature matrix and calculate liebman again
//         // if we are able to duplicate
//         if (currCols * 2 < finalCols && currRows * 2 < finalRows)
//         {
//             duplicateTempMatrix(currMatrix);
//             currCols = currCols * 2;
//             currRows = currRows * 2;
//         }
//         else if (currCols * 2 < finalCols)
//         {
//             duplicateWidthTempMatrix(currMatrix);
//             currCols = currCols * 2;
//         }
//         else if (currRows * 2 < finalRows)
//         {
//             duplicateHeightTempMatrix(currMatrix);
//             currRows = currRows * 2;
//         }
//         else
//         {
//             fillBiggerTempMatrix(currMatrix, finalRows, finalCols);
//             currCols = finalCols;
//             currRows = finalRows;
//         }

//         ++totalMatrices;
//     }

//     std::cout << "\nMatrices created: " << totalMatrices << std::endl;
//     //initial step Matrix size 1
//     return totalMatrices;
// }

// int LiebmnanSolver::liebman(Matrix<double> &temperatureMatrix, Edge top, Edge bottom, Edge right, Edge left)
// {
//     int iterations = 0;
//     int rows = temperatureMatrix.rows();
//     int cols = temperatureMatrix.cols();

//     //if Matrix of size 1
//     if (rows == 1 && cols == 1)
//     {
//         bool tbIsolated = false;
//         bool lrIsolated = false;
//         double tt, bt, lt, rt;

//         //if one of the sides is isolated
//         if ((top.isolated || bottom.isolated || right.isolated || left.isolated))
//         {
//             //top isolated
//             if (top.isolated)
//             {
//                 if (bottom.isolated)
//                 {
//                     // tt = 0;
//                     // bt = 0;
//                     tbIsolated = true;
//                 }
//                 else
//                 {
//                     tt = bottom.temperatures[0];
//                     bt = bottom.temperatures[0];
//                 }
//             }
//             //bottom isolated
//             else if (bottom.isolated)
//             {

//                 tt = top.temperatures[0];
//                 bt = top.temperatures[0];
//             }

//             //right isolated
//             if (right.isolated)
//             {
//                 if (left.isolated)
//                 {
//                     lrIsolated = true;
//                 }
//                 else
//                 {
//                     rt = left.temperatures[0];
//                     lt = left.temperatures[0];
//                 }
//             }
//             //left isolated
//             else if (left.isolated)
//             {
//                 tt = right.temperatures[0];
//                 bt = right.temperatures[0];
//             }
//         }
//         //no side is isolated
//         else
//         {
//             tt = top.temperatures[0];
//             bt = bottom.temperatures[0];
//             lt = left.temperatures[0];
//             rt = right.temperatures[0];
//         }

//         //if the plaque is completely isolated
//         if (lrIsolated && tbIsolated)
//         {
//             temperatureMatrix[0][0] = 0;
//         }
//         //both sides isolated
//         else if (lrIsolated)
//         {
//             temperatureMatrix[0][0] = (tt + bt) / 2;
//         }
//         //top and bottom isolated
//         else if (tbIsolated)
//         {
//             temperatureMatrix[0][0] = (rt + lt) / 2;
//         }
//         //calculate normally
//         else
//         {
//             temperatureMatrix[0][0] = (rt + lt + tt + bt) / 4;
//         }
//         //no need to iterate, will always get the same result
//         return 1;
//     } //end of matrix size 1

//     double bigTemp, bigTempOld;

//     double currTemp, calcTemp, calcError;
//     do
//     {
//         bigTemp = 0;
//         bigTempOld = 0;
//         //iterate through rows
//         for (int i = 0; i < rows; ++i)
//         {
//             //iterate through columns
//             for (int j = 0; j < cols; ++j)
//             {

//                 currTemp = temperatureMatrix[i][j];

//                 //////////////////////////////////top check///////////////////////////////////////////
//                 //if on the top
//                 if (i == 0)
//                 {

//                     //if on the left corner
//                     if (j == 0)
//                     {
//                         //if top is isolated
//                         if (top.isolated)
//                         {

//                             //if left is also isolated
//                             if (left.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i + 1][j] + 2 * temperatureMatrix[i][j + 1]) / 4;
//                             }
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
//                             }
//                         }
//                         //top not isolated
//                         else
//                         {
//                             //if left is isolated
//                             if (left.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + 2 * temperatureMatrix[i][j + 1]) / 4;
//                             }
//                             //neither is isolated
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
//                             }
//                         }
//                     }
//                     //if on the right corner
//                     else if (j == cols - 1)
//                     {
//                         //if top is isolated
//                         if (top.isolated)
//                         {

//                             //if right is also isolated
//                             if (right.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i + 1][j] + 2 * temperatureMatrix[i][j - 1]) / 4;
//                             }
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j - 1] + right.temperatures[i]) / 4;
//                             }
//                         }
//                         //top not isolated
//                         else
//                         {
//                             //if right is isolated
//                             if (right.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + 2 * temperatureMatrix[i][j - 1]) / 4;
//                             }
//                             //neither is isolated
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j - 1] + right.temperatures[i]) / 4;
//                             }
//                         }
//                     }
//                     //other parts of the top
//                     else
//                     {
//                         if (top.isolated)
//                         {
//                             //calculate temperature
//                             calcTemp = (2 * temperatureMatrix[i + 1][j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
//                         }
//                         else
//                         {
//                             //calculate temperature
//                             calcTemp = (temperatureMatrix[i + 1][j] + top.temperatures[j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
//                         }
//                     }
//                 }
//                 //////////////////////////////////bottom check///////////////////////////////////////////
//                 //if at the bottom
//                 else if (i == rows - 1)
//                 {
//                     //if on the left corner
//                     if (j == 0)
//                     {
//                         //if bottom is isolated
//                         if (bottom.isolated)
//                         {
//                             //if left is also isolated
//                             if (left.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i - 1][j] + 2 * temperatureMatrix[i][j + 1]) / 4;
//                             }
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
//                             }
//                         }
//                         //bottom not isolated
//                         else
//                         {
//                             //if left is isolated
//                             if (left.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + 2 * temperatureMatrix[i][j + 1]) / 4;
//                             }
//                             //neither is isolated
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j + 1] + left.temperatures[i]) / 4;
//                             }
//                         }
//                     }
//                     //if on the right corner
//                     else if (j == cols - 1)
//                     {
//                         //if bottom is isolated
//                         if (bottom.isolated)
//                         {
//                             //if right is also isolated
//                             if (right.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i - 1][j] + 2 * temperatureMatrix[i][j - 1]) / 4;
//                             }
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + right.temperatures[i]) / 4;
//                             }
//                         }
//                         //bottom not isolated
//                         else
//                         {
//                             //if right is isolated
//                             if (right.isolated)
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + 2 * temperatureMatrix[i][j - 1]) / 4;
//                             }
//                             //neither is isolated
//                             else
//                             {
//                                 //calculate temperature
//                                 calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j + 1] + right.temperatures[i]) / 4;
//                             }
//                         }
//                     }
//                     //other parts of the bottom
//                     else
//                     {
//                         if (bottom.isolated)
//                         {
//                             //calculate temperature
//                             calcTemp = (2 * temperatureMatrix[i - 1][j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
//                         }
//                         else
//                         {
//                             //calculate temperature
//                             calcTemp = (temperatureMatrix[i - 1][j] + bottom.temperatures[j] + temperatureMatrix[i][j - 1] + temperatureMatrix[i][j + 1]) / 4;
//                         }
//                     }
//                 }
//                 //corners were checked when calculating top and bottom

//                 //////////////////////////////////left check///////////////////////////////////////////
//                 //if on the left side non corner
//                 else if (j == 0)
//                 {
//                     if (left.isolated)
//                     {
//                         //calculate temperature
//                         calcTemp = (2 * temperatureMatrix[i][j + 1] + temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j]) / 4;
//                     }
//                     else
//                     {
//                         //calculate temperature
//                         calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + left.temperatures[i] + temperatureMatrix[i][j + 1]) / 4;
//                     }
//                 }

//                 //////////////////////////////////right check///////////////////////////////////////////
//                 //if on the right side non corner
//                 else if (j == cols - 1)
//                 {
//                     if (right.isolated)
//                     {
//                         //calculate temperature
//                         calcTemp = (2 * temperatureMatrix[i][j - 1] + temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j]) / 4;
//                     }
//                     else
//                     {
//                         //calculate temperature
//                         calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + right.temperatures[i] + temperatureMatrix[i][j - 1]) / 4;
//                     }
//                 }

//                 //////////////////////////////////other cells///////////////////////////////////////////
//                 //every other cell
//                 else
//                 {
//                     //calculate temperature
//                     calcTemp = (temperatureMatrix[i + 1][j] + temperatureMatrix[i - 1][j] + temperatureMatrix[i][j + 1] + temperatureMatrix[i][j - 1]) / 4;
//                 }

//                 //aplyrelaxation
//                 calcTemp = lambda * calcTemp + (1 - lambda) * currTemp;
//                 //set the new temperture
//                 temperatureMatrix[i][j] = calcTemp;

//                 //check if the new calculated temp is the biggest to compare for error
//                 if (calcTemp > bigTemp)
//                 {
//                     bigTemp = calcTemp;
//                     bigTempOld = currTemp;
//                 }

//             } //end for column
//         }     //end for row
//         //////////////////////////////////end matrix iteration///////////////////////////////////////////

//         ++iterations;

//         //calculate error of the biggest temperature
//         calcError = (abs(bigTemp - bigTempOld) / bigTemp) * 100;

//     } while (calcError > error);

//     return iterations;
// }

// void LiebmnanSolver::duplicateTempMatrix(Matrix<double> &tempMatrix)
// {
//     int scols = tempMatrix.cols();
//     int srows = tempMatrix.rows();
//     int bcols = scols * 2;
//     int brows = srows * 2;
//     Matrix<double> big(brows, bcols);

//     int bigI, bigJ;
//     for (int i = 0; i < srows; ++i)
//         for (int j = 0; j < scols; ++j)
//         {
//             bigI = i * 2;
//             bigJ = j * 2;
//             big[bigI][bigJ] = tempMatrix[i][j];
//             big[bigI + 1][bigJ] = tempMatrix[i][j];
//             big[bigI][bigJ + 1] = tempMatrix[i][j];
//             big[bigI + 1][bigJ + 1] = tempMatrix[i][j];
//         }

//     tempMatrix = big;
// }
// void LiebmnanSolver::duplicateWidthTempMatrix(Matrix<double> &tempMatrix)
// {
//     int scols = tempMatrix.cols();
//     int srows = tempMatrix.rows();
//     int bcols = scols * 2;
//     int brows = srows;
//     Matrix<double> big(brows, bcols);

//     int bigI, bigJ;
//     for (int i = 0; i < srows; ++i)
//         for (int j = 0; j < scols; ++j)
//         {
//             bigI = i;
//             bigJ = j * 2;
//             big[bigI][bigJ] = tempMatrix[i][j];
//             big[bigI][bigJ + 1] = tempMatrix[i][j];
//         }
//     tempMatrix = big;
// }
// void LiebmnanSolver::duplicateHeightTempMatrix(Matrix<double> &tempMatrix)
// {
//     int scols = tempMatrix.cols();
//     int srows = tempMatrix.rows();
//     int bcols = scols;
//     int brows = srows * 2;

//     Matrix<double> big(brows, bcols);

//     int bigI, bigJ;
//     for (int i = 0; i < srows; ++i)
//         for (int j = 0; j < scols; ++j)
//         {
//             bigI = i * 2;
//             bigJ = j;
//             big[bigI][bigJ] = tempMatrix[i][j];
//             big[bigI + 1][bigJ] = tempMatrix[i][j];
//         }
//     tempMatrix = big;
// }

// /**
//      * @brief Fills a matrix that is less than double the size of another matrix used to fill it. It copies
//      * the last rigthmost columns of the smaller matrix  into the rightmost columns of the big matirx, it
//      * does similar to the bottom part, but the bottom part of the rigmtost section uses the values set when
//      * filling the rightmost section
//      *
//      * @param small The matrix used to fill the new matrix, and also were this new matrix is returned
//      * @param rows The row size of the new bigger matrix
//      * @param cols The column size of the new bigger matrix
//      */
// void LiebmnanSolver::fillBiggerTempMatrix(Matrix<double> &small, int rows, int cols)
// {
//     int scols = small.cols();
//     int srows = small.rows();
//     int colDiff = cols - scols;
//     int rowDiff = rows - srows;

//     if (rows < srows || cols < scols)
//     {
//         throw anpi::Exception("new Matrix has to have bigger indices");
//     }

//     Matrix<double> big(rows, cols);

//     big.fill(small);

//     for (int i = 0; i < srows; ++i)
//         for (int j = scols; j < cols; ++j)
//         {
//             big[i][j] = small[i][j - colDiff];
//         }
//     for (int i = srows; i < rows; ++i)
//         for (int j = 0; j < cols; ++j)
//         {
//             if (j < scols)
//             {
//                 big[i][j] = small[i - rowDiff][j];
//             }
//             else
//             {
//                 big[i][j] = small[i - rowDiff][j - colDiff];
//             }
//         }

//     small = big;
// }

// void LiebmnanSolver::fillBiggerTempMatrix(Matrix<double> &small, Matrix<double> &big)
// {
//     int scols = small.cols();
//     int bcols = big.cols();
//     int srows = small.rows();
//     int brows = big.rows();

//     //matrix size was duplicated each old cell corresponds to 4 cells of the new one
//     if (scols * 2 == bcols && srows * 2 == brows)
//     {
//         int bigI, bigJ;
//         for (int i = 0; i < srows; ++i)
//             for (int j = 0; j < scols; ++j)
//             {
//                 bigI = i * 2;
//                 bigJ = j * 2;
//                 big[bigI][bigJ] = small[i][j];
//                 big[bigI + 1][bigJ] = small[i][j];
//                 big[bigI][bigJ + 1] = small[i][j];
//                 big[bigI + 1][bigJ + 1] = small[i][j];
//             }
//     }
//     else
//     {
//         //we just fill the bigger matrix with the small one
//         big.fill(small);
//     }
// }

// Edge LiebmnanSolver::averageEdgeTemp(Edge &temps, int size)
// {
//     //if it's isolated no need to work
//     if (temps.isolated)
//         return temps;

//     int tempSize = temps.temperatures.size();

//     //if the size is the same, just return the temps
//     if (tempSize == size)
//     {
//         return temps;
//     }

//     Edge avgEdge;
//     avgEdge.isolated = false;
//     std::vector<double> avgTemps;

//     //if it's all the same temperature
//     if (!temps.gradient)
//     {
//         avgTemps.resize(size, temps.temperatures[0]);
//         avgEdge.temperatures = avgTemps;
//         return avgEdge;
//     }
//     else
//     {
//         avgEdge.gradient = true;
//     }

//     int chunkSize = tempSize / size;
//     int overflow = tempSize % size;
//     int end, j;
//     double average;

//     for (int i = 0; i <= size; ++i)
//     {
//         average = 0;
//         //if on the last chunk we also average the leftover temperatures
//         if (i == size - 1)
//         {
//             j = i * chunkSize;
//             end = j + chunkSize + overflow;
//             for (; j < end; ++j)
//             {
//                 average += temps.temperatures[j];
//             }
//             average = average / (chunkSize + overflow);
//         }
//         else
//         {
//             j = i * chunkSize;
//             end = j + chunkSize;
//             for (; j < end; ++j)
//             {
//                 average += temps.temperatures[j];
//             }
//             average = average / chunkSize;
//         }

//         avgTemps.push_back(average);
//     }

//     avgEdge.temperatures = avgTemps;
//     avgEdge.isolated = false;
//     return avgEdge;
// }
// } // namespace anpi

