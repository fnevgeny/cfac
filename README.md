# What's in a name?

cFAC was started around 2010 as a fork of [Flexible Atomic Code (FAC)](https://github.com/flexible-atomic-code/fac), written by Ming Feng Gu at the Space Science Laboratory of Berkeley.

The initial focus was on providing large volumes of data as required, e.g., for
collisional-radiative plasma modeling, and eliminating reliance upon third-party
Fortran numerical libraries with their C equivalents (hence the change in the
package name).

# Licensing

FAC, and hence cFAC, are released under the [GPL (version 3 or higher) license](https://www.gnu.org/licenses/gpl-3.0.en.html).

However, some bits of FAC which are still used in cFAC, were published in
[Computer Physics Communications](https://www.sciencedirect.com/journal/computer-physics-communications), and as such, [are licensed for
non-profit or academic use only](https://www.elsevier.com/about/policies-and-standards/open-access-licenses/elsevier-user/cpc).

In order to compile in these CPC-licensed modules, you need to pass the
"\-\-with-cpc-modules" configure flag and explicitly agree to the CPC licensing
terms. Please note that as a result, the "sfac" executable will *not* be redistributable.

# Installation

## Prerequisites

  - C and FORTRAN compilers (such as gcc and gfortran).

  - The [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/),
   version 1.15 or above.

  - The [SQLite software library](http://www.sqlite.org/).

      - If you install these software using a package manager of your OS, make sure you select the development packages (\*-dev\*).
        For example, under a Debian(-based) Linux, do

        ```
        sudo apt install gcc gfortran libgsl0-dev libsqlite3-dev
        ```

  - In addition, if compiling a snapshot from GitHub (strongly **not** recommended for "production" calculations!), GNU Automake, help2man, and a fairly complete LaTeX installation will be required.

## Compilation

  1. Get the latest _release_ from <https://github.com/fnevgeny/cfac/releases>.
   The filename should be like cfac-1.7.0.tar.gz.
   Unpack it and change to the newly created directory.

  2. ./configure

        - Specify \-\-prefix=my/dir, if the default /usr/local is not what you want.

        - Use the "\-\-with-extrainc" and "\-\-with-extralib" options if the third-party libraries are not installed in the default system places.

        - Run ./configure \-\-help for more options.

        - Note in particular the "\-\-with-cpc-modules" option (see [above](#licensing)).

  3. make

  4. make check

      - This will run a test suite of demo scripts and compare their output against "reference" results.
      **Please report any failure**.

  5. make install

      - This installs the "sfac" executable and manual, as well as the CFACDB library, C header file, and utility.

# Usage
The "sfac" executable accepts an input file as its only argument. See the "demo"
directory for examples. For more details, please read the manual, and **make sure to
read the FAQ** section of it.

Evgeny Stambulchik  
Faculty of Physics  
Weizmann Institute of Science  
Rehovot 7610001  
Israel
