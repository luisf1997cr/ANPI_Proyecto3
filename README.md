# ANPI_Proyecto3

Proyecto 3 de ANPI
Jorge Agüero Zamora
Luis Fernando Murillo 


*****************************************************************************************************************************
**********************************************  Proyecto 3 ***************e**************************************************
*****************************************************************************************************************************

README

In this project we calculate the temperature flow in a metal plaque usign the Liebman
method to solve the Prial differential equations of the problem

CONTACT

If you have problems or comments with this program you
can contact please contact us.

This project can also be found at GitHub in:
https://github.com/luisf1997cr/ANPI_Proyecto3

-----------------------------------------------------------------------------------------------------------
------------------------------------------- BUILD INSTRUCTIONS --------------------------------------------
-----------------------------------------------------------------------------------------------------------

Unzip the project, open a terminal and change your working directory to the unzipped folder

Create a directory build:

> mkdir build;

An ls should look something like
> ls
build  code  README.md

Go into the build directory

> cd build;

You can choose to build a release version with:

> cmake ../code -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../code -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make


----------------------------------------------------------------------------------------------------------
------------------------------------------- RUN INSTRUCTIONS ---------------------------------------------
----------------------------------------------------------------------------------------------------------

The executables will be stored at build/bin.

To execute the program go to the build/bin directory. 

Suposing you are currently on the build directory

> cd bin

To run the program you can use

>./proyecto3

But that doesn do much, go ahead and try

>./proyecto3 --help

That will give you a list of the available options you can use.

The project provide a flag that can be used to compare the different options by using the --test flag

>./proyecto3 --test

Have fun! 
