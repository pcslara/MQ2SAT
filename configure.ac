#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.68])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])

AC_PROG_CC
	
# Checks for header files.
AC_FUNC_ALLOCA
                     
AC_CHECK_HEADERS([malloc.h math.h stdio.h stdlib.h string.h unistd.h])
AC_C_INLINE
AC_TYPE_INT32_T


AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([fopen fgets printf fprintf memset strcpy strcat strlen])

AC_CHECK_PROG(ANF2CNF,anf2cnf,yes)
if test x"$ANF2CNF" != x"yes" ; then
    AC_MSG_ERROR([Please install anf2cnf before installing.])
fi

AC_CHECK_PROG(MINISAT,minisat,yes)
if test x"$MINISAT" != x"yes" ; then
    AC_MSG_ERROR([Please install minisat before installing.])
fi


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
