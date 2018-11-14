# adapted from configure.ac used in rgeos package (Roger Bivand)
define([pkgversion], esyscmd([sh -c "grep Version: DESCRIPTION | cut -d' ' -f2 | tr -d '\n'"]))

AC_INIT(exactextractr, [pkgversion], dbaston@isciences.com)
AC_MSG_NOTICE([${PACKAGE_NAME}: ${PACKAGE_VERSION}])

AC_CONFIG_SRCDIR(src/exact_extract.cpp)

# find R home and set correct compiler + flags
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([cannot determine R_HOME. Make sure you use R CMD INSTALL!])
fi
RBIN="${R_HOME}/bin/R"

# pick all flags for testing from R
: ${CXX=`"${RBIN}" CMD config CXX14`}
: ${CXXFLAGS=`"${RBIN}" CMD config CXX14FLAGS` `"${RBIN}" CMD config CXX14STD`}
: ${LDFLAGS=`"${RBIN}" CMD config LDFLAGS`}

GEOS_CONFIG="geos-config"

GEOS_CONFIG_SET="no"

AC_ARG_WITH([geos-config],
  AC_HELP_STRING([--with-geos-config=GEOS_CONFIG],
                 [the location of geos-config]),
                 [geos_config=$withval])

  if test [ -n "$geos_config" ] ; then
    GEOS_CONFIG_SET="yes"
    AC_SUBST([GEOS_CONFIG],["${geos_config}"])
    AC_MSG_NOTICE(geos-config set to $GEOS_CONFIG)
fi

if test ["$GEOS_CONFIG_SET" = "no"] ; then
  AC_PATH_PROG([GEOS_CONFIG], ["$GEOS_CONFIG"], ["no"])
  if test ["$GEOS_CONFIG" = "no"] ; then
    AC_MSG_ERROR([geos-config not found or not executable])
  fi
else
  AC_MSG_CHECKING(geos-config exists)
  if test -r "${GEOS_CONFIG}"; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([geos-config not found - configure argument error.])
  fi
  AC_MSG_CHECKING(geos-config executable)
  if test -x "${GEOS_CONFIG}"; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([geos-config not executable.])
  fi
fi

AC_MSG_CHECKING(geos-config usability)
if test `${GEOS_CONFIG} --version`;
then
  GEOS_VER=`${GEOS_CONFIG} --version`
  GEOS_VER_DOT=`${GEOS_CONFIG} --version | sed 's/[[^0-9]]*//g'`
  GEOS_CXXFLAGS=`${GEOS_CONFIG} --cflags`
  GEOS_LDFLAGS=`${GEOS_CONFIG} --clibs`
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([${GEOS_CONFIG} not usable])
fi

AC_MSG_NOTICE([GEOS version: ${GEOS_VER}])
AC_MSG_CHECKING([geos version at least 3.6.0])
if test ${GEOS_VER_DOT} -lt 360 ; then
  AC_MSG_RESULT(no)
  AC_MSG_RESULT([Upgrade GEOS to version 3.6.0 or greater.])
else
  AC_MSG_RESULT(yes)
fi

CXXFLAGS="${CXXFLAGS} ${GEOS_CXXFLAGS}"
LDFLAGS="${LDFLAGS} ${GEOS_LDFLAGS}"

AC_MSG_CHECKING(compiling and building against geos_c)
[cat > geos_test.cpp << _EOCONF
#include <geos_c.h>
#include <iostream>

int main() {
  auto handle = initGEOS_r(nullptr, nullptr);
  finishGEOS_r(handle);
  return 0;
}
_EOCONF]

${CXX} ${CXXFLAGS} -o geos_test geos_test.cpp ${LDFLAGS} 2> errors.txt
if test `echo $?` -ne 0 ; then
  AC_MSG_RESULT(no)
  cat errors.txt
  AC_MSG_NOTICE([Compilation and/or linkage problems.])
else
  AC_MSG_RESULT(yes)
fi

rm -f geos_test errors.txt geos_test.cpp

AC_SUBST([PKG_CXX], ["${CXX}"])
AC_SUBST([PKG_CXXFLAGS], ["${CXXFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LDFLAGS}"])

AC_MSG_NOTICE([PKG_CXX: ${PKG_CXX}])
AC_MSG_NOTICE([PKG_CXXFLAGS: ${PKG_CXXFLAGS}])
AC_MSG_NOTICE([PKG_LIBS: ${PKG_LIBS}])

AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT

