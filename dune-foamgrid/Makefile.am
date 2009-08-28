
SUBDIRS = foamgrid test

# the standard debug streams are put into the libdune
#noinst_LTLIBRARIES = libcommon.la

#libcommon_la_SOURCES = stdstreams.cc configparser.cc

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/..

subdirincludedir = $(includedir)/dune/foamgrid
subdirinclude_HEADERS = foamgrid.hh

include $(top_srcdir)/am/global-rules
