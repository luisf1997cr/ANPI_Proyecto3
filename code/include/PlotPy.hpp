/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * Plotting methods, wrapping Matplotlib (Python)
 * @author: David Badilla S.
 * @email: davbs94@gmail.com
 * @date:   15.07.2017
 */

#ifndef ANPI_PLOTPY_HPP
#define ANPI_PLOTPY_HPP

#include <python2.7/Python.h>
#include <stdlib.h>
#include <string>
#include <vector>

namespace anpi
{

/**
   * Two-dimensional plots
   *
   * Given a set of x coordinates and a corresponding set of y values,
   * plot a line between each point (x,y).
   *
   * You give a pair of vectors with the x and y values with the
   * plot() method to plot a curve.  You may overlay as many curves as
   * you need simply by calling plot() as many times as you need to.
   *
   * Finally, you call show() to display the window with all plotted curves.
   */
template <typename T>
class Plot2d
{
private:
  //Titulo de la grafica.
  std::string _title;
  //Estiqueta del eje x
  std::string _xlabel;
  //Etiqueta del eje y
  std::string _ylabel;
  //Tamano de la cuadricula
  T _sizeGrid;

public:
  /// Constructors
  //@{
  Plot2d();
  ~Plot2d();
  //@}

  /**
     * Initialize a plot window.
     *
     * Each id is associated with a different plot window.
     */
  void initialize();

  /// Set plot title
  void setTitle(const std::string &title);

  /// Set label for the X axis
  void setXLabel(const std::string &label);
  /// Set label for the Y axis
  void setYLabel(const std::string &label);

  /// Set the grid size
  void setGridSize(const T sizegrid);

  /// Set initial and final limits of the X axis
  void setXRange(const T xi, const T xs);

  /// Set initial and final limits of the Y axis
  void setYRange(const T yi, const T ys);

  /**
     * Plot a curve by drawing line segments from
     * the sequence of points (datax[i],datay[i]).  The
     * curve will have the given legend 
     */
  void plot(const std::vector<T> &datax,
            const std::vector<T> &datay,
            const std::string &legend,
            const std::string &color = "");

  /**
     * Plot an area range between the min and max values
     * and the average data inbetween.
     * @param datax values of x
     * @param averagey average values of y, corresponding to each x
     * @param miny minimum values of y, corresponding to each x
     * @param maxy maximum values of y, corresponding to each x
     */
  void plot(const std::vector<T> &datax,
            const std::vector<T> &averagey,
            const std::vector<T> &miny,
            const std::vector<T> &maxy,
            const std::string &legend,
            const std::string &color = "r");

  /**
     * Plot an area using quiver to print arrows
     * @param datax values of x
     * @param datay values of y 
     * @param datax values of u 
     * @param datav values of v
     */
  void quiver(std::vector<T> &datax,
              std::vector<T> &datay,
              std::vector<T> &datau,
              std::vector<T> &datav);

  void quiver(anpi::Matrix<T> &image);
  /**
     * Plot a image using interpolation
     * the image comes as a matrix
     * @param Matrix of element
     */

  void imgshow(anpi::Matrix<T> &);

  void imgshow(anpi::Matrix<T> &image, std::string title);

  /**
     * Show all curves plotted so far.
     */
  void show();

}; //class Plot2d

} // namespace anpi

#include "PlotPy.tpp"

#endif // PLOTPY_H
