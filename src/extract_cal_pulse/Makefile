

include ../../make_support/system_rules

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


OBJS = create_roi_in.o \
	dop.o \
	spectra.o \
	fix_vectors.o \
	get_peg.o \
	tle_code.o \
	runtest.o \
	esa_nodes.o

all: simple_amp_pulse extract_cal_pulse simple_cal_pulse extract_sinc

simple_amp_pulse: spectra.o fft_oversamp.o match_pulse.o simple_amp_pulse.o
	cc $(CFLAGS) $(LFLAGS) -o simple_amp_pulse simple_amp_pulse.o $(LIBS) $(LDFLAGS) $(INCLUDES)
	mv simple_amp_pulse $(BINDIR)

extract_cal_pulse: spectra.o fft_oversamp.o match_pulse.o extract_cal_pulse.o
	cc $(CFLAGS) $(LFLAGS) -o extract_cal_pulse spectra.o fft_oversamp.o match_pulse.o extract_cal_pulse.o $(LIBS) $(LDFLAGS) $(INCLUDES)
	mv extract_cal_pulse $(BINDIR)

simple_cal_pulse: spectra.o fft_oversamp.o match_pulse.o simple_cal_pulse.o
	cc $(CFLAGS) $(LFLAGS) -o simple_cal_pulse simple_cal_pulse.o $(LIBS) $(LDFLAGS) $(INCLUDES)
	mv simple_cal_pulse $(BINDIR)

extract_sinc: match_sinc.o extract_sinc.o
	cc $(CFLAGS) $(LFLAGS) -o extract_sinc match_sinc.o extract_sinc.o $(LIBS) $(LDFLAGS) $(INCLUDES)
	mv extract_sinc $(BINDIR)

fft_ov_stub: fft_oversamp.o fft_ov_stub.o
	cc $(CFLAGS) $(LFLAGS) -o fft_ov_stub fft_oversamp.o fft_ov_stub.o $(LIBS) $(LDFLAGS) $(INCLUDES)
	mv fft_ov_stub $(BINDIR)
clean:
	rm -f *.o
