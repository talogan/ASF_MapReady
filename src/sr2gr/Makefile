# Makefile for          :  sr2gr
# Module Author         :  Tom Logan 
# Module Version        :  2.0

include ../../make_support/system_rules

CFLAGS += $(GLIB_CFLAGS)

LIBS  = $(LIBDIR)/libasf_sar.a \
	$(LIBDIR)/libasf_raster.a \
	$(LIBDIR)/asf_meta.a \
	$(GSL_LIBS) \
	$(LIBDIR)/libasf_proj.a \
	$(LIBDIR)/asf.a \
	$(PROJ_LIBS) \
	$(XML_LIBS) \
	$(GLIB_LIBS) \
	$(ZLIB_LIBS) \
	-lm

OBJS  = sr2gr.o help.o

all: sr2gr clean

sr2gr : $(OBJS)
	$(CC) $(CFLAGS) -o sr2gr $(OBJS) $(LIBS) $(LDFLAGS)
	mv sr2gr$(BIN_POSTFIX) $(BINDIR)
	cp sr2gr.1 ../../man/cat1

clean:
	rm -f *.o
