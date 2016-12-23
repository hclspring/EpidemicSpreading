#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([src/dynamic_network.cc])
AC_CONFIG_HEADERS([config.h])

CFLAGS=-O3
CXXFLAGS="-O3 -std=c++11"

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [main])
# FIXME: Replace `main' with a function in `-ligraph':
AC_CHECK_LIB([igraph], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
