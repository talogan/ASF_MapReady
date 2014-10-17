# Makefile for meta2ddr/ddr2meta
# Author Patrick Denny

include ../../make_support/system_rules

INCLUDES = -I../../include/libsgp4

LIBS  = \
	$(LIBDIR)/libasf_sar.a \
	$(LIBDIR)/libasf_raster.a \
	$(LIBDIR)/asf_meta.a \
	$(GSL_LIBS) \
	$(LIBDIR)/libasf_proj.a \
	$(LIBDIR)/asf.a \
	$(PROJ_LIBS) \
	$(XML_LIBS) \
	$(LIBDIR)/libfftw3f.a \
	$(LIBDIR)/libsgp4.a \
	-lm -lz

OBJS = ksat_ingest.o

all: ksat_ingest

runtest.o: runtest.cpp
	c++ -c runtest.cpp -I../../include/libsgp4

ksat_ingest: $(OBJS)
	cc $(CFLAGS) $(LFLAGS) -o ksat_ingest $(OBJS) $(LIBS) $(LDFLAGS) $(INCLUDES)
	mv ksat_ingest$(BIN_POSTFIX) $(BINDIR)

clean:
	rm -f $(OBJS)
