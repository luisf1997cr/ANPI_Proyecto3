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

namespace anpi
{

template <typename T>
Plot2d<T>::Plot2d() {}

template <typename T>
Plot2d<T>::~Plot2d() {}

template <typename T>
void Plot2d<T>::initialize()
{
  Py_Initialize();
  PyRun_SimpleString("import matplotlib.pyplot as plt");
  PyRun_SimpleString("import numpy as np");
  PyRun_SimpleString("from mpl_toolkits.axes_grid1 import make_axes_locatable");
}

template <typename T>
void Plot2d<T>::setTitle(const std::string &title)
{
  _title = title;
  std::string strtitle = "plt.title('" + _title + "')";
  PyRun_SimpleString(strtitle.c_str());
}

template <typename T>
void Plot2d<T>::quiver(anpi::Matrix<T> &image)
{

  // std::string FluxFun = "def Flux(Matriz):\n";
  // auxFun.append("	k = 1\n");
  // auxFun.append("	delX = 1\n");
  // auxFun.append("	delY = 1\n");
  // auxFun.append("	if((len(Z[0])+len(Z))//(2)>=10):\n");
  // auxFun.append("		n=(len(Z[0])+len(Z))//(2*10)\n");
  // auxFun.append("	else:\n");
  // auxFun.append("		n=1\n");
  // auxFun.append("	dx = np.linspace(0,(len(Z[0])-n),len(Z[0])//n)	\n");
  // auxFun.append("	dy = np.linspace(0,(len(Z)-n),len(Z)//n)	\n");
  // auxFun.append("	X,Y = np.meshgrid(dx,dy)\n");
  // auxFun.append("	u = np.zeros((len(X),len(X[0])))\n");
  // auxFun.append("	v = np.zeros((len(X),len(X[0])))\n");
  // auxFun.append("	a,b=0,0\n");
  // auxFun.append("	for x in range(len(X)):\n");
  // auxFun.append("		for y in range(len(X[0])):\n");
  // auxFun.append("			if(a!=0 and a!=len(Z)-1 and b!=0 and b!=len(Z[0])-1):\n");
  // auxFun.append("				v[x][y]= k*(Z[a+1][b]-Z[a-1][b])/(2*delX)\n");
  // auxFun.append("				u[x][y]= -k*(Z[a][b+1]-Z[a][b-1])/(2*delY)	\n");
  // auxFun.append("			else:\n");
  // auxFun.append("				v[x][y]= 0\n");
  // auxFun.append("				u[x][y]= 0\n");
  // auxFun.append("			b+=n\n");
  // auxFun.append("		a+=n\n");
  // auxFun.append("		b=0\n");
  // auxFun.append("	return [X,Y,u,v]\n");
  // PyRun_SimpleString(FluxFun.c_str());
  // PyRun_SimpleString("FluxMatriz = Flux(image)");
  // PyRun_SimpleString("ax.quiver(FluxMatriz[0],FluxMatriz[1],FluxMatriz[2],FluxMatriz[3])");

  // Convert the vectors of data into Python strings
  //std::string xstr  = "X = [";
  //std::string ystr  = "Y = [";
  //std::string ustr  = "U = [";
  //std::string vstr  = "V = [";

  //char c=',';
  //for(size_t i = 0; i < datax.size(); i++) {
  //  if (i == datax.size()-1) {
  //    c=']';
  //  }
  //  xstr.append(std::to_string(datax[i])   + c);
  //  ystr.append(std::to_string(datay[i])   + c);

  //}
  //c=',';
  //for(size_t i = 0; i < datav.size(); i++) {
  //  if (i == datav.size()-1) {
  //    c=']';
  //  }
  //  ustr.append(std::to_string(datau[i])   + c);
  //  vstr.append(std::to_string(datav[i])   + c);
  //}
  //c=',';

  //PyRun_SimpleString(xstr.c_str());
  //PyRun_SimpleString(ystr.c_str());
  //PyRun_SimpleString(ustr.c_str());
  //PyRun_SimpleString(vstr.c_str());
  //PyRun_SimpleString("ax.quiver(X, Y, U, V)");
}

template <typename T>
void Plot2d<T>::setXLabel(const std::string &xlabel)
{
  _xlabel = xlabel;
  std::string strxlabel = "plt.xlabel('" + _xlabel + "')";
  PyRun_SimpleString(strxlabel.c_str());
}

template <typename T>
void Plot2d<T>::setYLabel(const std::string &ylabel)
{
  _ylabel = ylabel;
  std::string strylabel = "plt.ylabel('" + _ylabel + "')";
  PyRun_SimpleString(strylabel.c_str());
}

template <typename T>
void Plot2d<T>::setGridSize(const T sizegrid)
{
  _sizeGrid = sizegrid;
  std::string strgrid = "plt.grid(" + std::to_string(_sizeGrid) + ")";
  PyRun_SimpleString(strgrid.c_str());
}

template <typename T>
void Plot2d<T>::setXRange(const T xi, const T xs)
{
  std::string strxlim = "plt.xlim(" +
                        std::to_string(xi) + "," + std::to_string(xs) + ")";
  PyRun_SimpleString(strxlim.c_str());
}

template <typename T>
void Plot2d<T>::setYRange(const T yi, const T ys)
{
  std::string strylim = "plt.ylim(" +
                        std::to_string(yi) + "," + std::to_string(ys) + ")";
  PyRun_SimpleString(strylim.c_str());
}

template <typename T>
void Plot2d<T>::plot(const std::vector<T> &datax,
                     const std::vector<T> &datay,
                     const std::string &label,
                     const std::string &color)
{
  std::string xstr = "datax = [";
  std::string ystr = "datay = [";
  std::string pltcmd = "plt.plot(datax,datay";
  if (!label.empty())
  {
    pltcmd += ",label='" + label + "'";
  }
  if (!color.empty())
  {
    pltcmd += ",color='" + color + "'";
  }
  pltcmd += ")";

  for (size_t i = 0; i < datax.size(); i++)
  {
    if (i == datax.size() - 1)
    {
      xstr.append(std::to_string(datax[i]) + "]");
      ystr.append(std::to_string(datay[i]) + "]");
    }
    else
    {
      xstr.append(std::to_string(datax[i]) + ",");
      ystr.append(std::to_string(datay[i]) + ",");
    }
  }

  PyRun_SimpleString(xstr.c_str());
  PyRun_SimpleString(ystr.c_str());
  PyRun_SimpleString(pltcmd.c_str());
  PyRun_SimpleString("plt.legend()");
}

template <typename T>
void Plot2d<T>::imgshow(anpi::Matrix<T> &image)
{
  std::string xstr = "x = [[";
  std::string c = ",";
  for (size_t i = 0; i < image.rows(); i++)
  {
    for (size_t j = 0; j < image.cols(); j++)
    {
      if (j == image.cols() - 1 && i != image.rows() - 1)
      {
        c = "],[";
      }
      if (j == image.cols() - 1 && i == image.rows() - 1)
      {
        c = "]]";
      }
      xstr.append(std::to_string(image[i][j]) + c);
      c = ",";
    }
  }

  PyRun_SimpleString(xstr.c_str());
  PyRun_SimpleString("fig, ax = plt.subplots(nrows=1, sharex=True, figsize=(10, 10))");
  PyRun_SimpleString("ax.set_title('placa')");
  PyRun_SimpleString("im = ax.imshow(x, origin='upper', interpolation='bilinear',cmap='plasma')");

  PyRun_SimpleString("divider = make_axes_locatable(ax)");
  PyRun_SimpleString("cax = divider.append_axes('right', size='5%', pad=0.05)");
  PyRun_SimpleString("plt.colorbar(im, cax=cax)");
}

template <typename T>
void Plot2d<T>::plot(const std::vector<T> &datax,
                     const std::vector<T> &averagey,
                     const std::vector<T> &miny,
                     const std::vector<T> &maxy,
                     const std::string &legend,
                     const std::string &color)
{

  // Convert the vectors of data into Python strings
  std::string xstr = "datax = [";
  std::string avgystr = "avgy = [";
  std::string minystr = "miny = [";
  std::string maxystr = "maxy = [";

  char c = ',';
  for (size_t i = 0; i < datax.size(); i++)
  {
    if (i == datax.size() - 1)
    {
      c = ']';
    }
    xstr.append(std::to_string(datax[i]) + c);
    avgystr.append(std::to_string(averagey[i]) + c);
    minystr.append(std::to_string(miny[i]) + c);
    maxystr.append(std::to_string(maxy[i]) + c);
  }

  std::string lstr = legend.empty()
                         ? ""
                         : std::string(",label='") + legend + "'";

  std::string cstr = color.empty()
                         ? ""
                         : std::string(",color='") + color + "'";

  std::string pltcmd = "plt.plot(datax,avgy" + lstr + cstr + ",lw=2)";

  std::string fillcmd = "plt.fill_between(datax,miny,maxy" + cstr + ",alpha=0.1)";

  // Python lines with the data
  PyRun_SimpleString(xstr.c_str());
  PyRun_SimpleString(avgystr.c_str());
  PyRun_SimpleString(minystr.c_str());
  PyRun_SimpleString(maxystr.c_str());

  // Plot the lines
  PyRun_SimpleString(pltcmd.c_str());
  PyRun_SimpleString(fillcmd.c_str());
  PyRun_SimpleString("plt.legend()");
}

template <typename T>
void Plot2d<T>::show()
{
  PyRun_SimpleString("plt.show()");
}

} // namespace anpi
