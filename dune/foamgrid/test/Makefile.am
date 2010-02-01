
TESTPROGS = foamgrid-test

# which tests to run
TESTS = $(TESTPROGS)

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)

AM_LDFLAGS = @AM_LDFLAGS@
AM_CXXFLAGS = -Wall -g

# define the programs
foamgrid_test_SOURCES = foamgrid-test.cc

include $(top_srcdir)/am/global-rules
