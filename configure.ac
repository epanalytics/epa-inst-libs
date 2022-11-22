#
# Process this file with autoconf to produce a configure script.
#/* 
# * This file is part of the pebil project.
# * 
# * Copyright (c) 2010, University of California Regents
# * All rights reserved.
# * 
# * This program is free software: you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation, either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program.  If not, see <http://www.gnu.org/licenses/>.
# */

AC_PREREQ([2.63])
AC_INIT(epa-inst-libs, 1.0, allysonc@epanalytics.com)
AC_CONFIG_SRCDIR([InstrumentationCommon.hpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

# Get the following flags from a tool configuration
AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(MPI_FLAGS)
AC_SUBST(MPI_CXXFLAGS)
AC_SUBST(PEBIL_DEBUG)
AC_SUBST(PEBIL_EPA_TOOLS_AVAILABLE)
AC_SUBST(PEBIL_SUBDIR_EPA_TOOLS_PATH)
AC_SUBST(PEBIL_PAPI_BUILD)
AC_ARG_VAR(PAPI_PATH,[path to PAPI folder])
AC_SUBST(CPUFREQ_FLAGS)


AC_CHECK_PROG(OPENMPI_IS_NEW, mpifort, "yes")
if test x"$OPENMPI_IS_NEW" == x"yes"; then
  INST_LIBS_MPIFC_OPENMPI_LIB="-lmpi_mpifh"
else
  INST_LIBS_MPIFC_OPENMPI_LIB="-lmpi_f77 -lmpi_f90"
fi
AC_SUBST(INST_LIBS_MPIFC_OPENMPI_LIB)


AC_ARG_WITH(throttler_headers,
            [AS_HELP_STRING([--with-throttler-headers=<path>],[specify path to throttler headers])],
            [throttler_headers=$withval
             AC_CHECK_HEADER([$throttler_headers/throttler.h],
                             [THROTTLER_FLAGS="-DHAVE_THROTTLER_H -I$throttler_headers"],
                             [AC_MSG_ERROR(throttler.h not found using $throttler_headers)])],
            [throttler_headers=no])
AC_SUBST(THROTTLER_FLAGS)




# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h strings.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([bzero memset])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([ReuseDistance])
AC_OUTPUT

exit 0
