# Makefile generated by imake - do not edit!

# ----------------------------------------------------------------------
# Makefile generated from "Imake.tmpl" and <Imakefile>
# $Xorg: Imake.tmpl,v 1.4 2000/08/17 19:41:46 cpqbld Exp $
# $XdotOrg: xc/config/cf/Imake.tmpl,v 1.9 2005/01/24 06:37:31 daniels Exp $
#
#
#
#
# $XFree86: xc/config/cf/Imake.tmpl,v 3.155 2003/12/24 18:58:41 dickey Exp $
# ----------------------------------------------------------------------

all::

.SUFFIXES: .i

# $XdotOrg: xc/config/cf/Imake.cf,v 1.7 2005/03/02 11:20:29 gisburn Exp $
# $Xorg: Imake.cf,v 1.4 2000/08/17 19:41:45 cpqbld Exp $

# $XFree86: xc/config/cf/Imake.cf,v 3.88 2003/12/16 21:30:21 herrb Exp $

# -----------------------------------------------------------------------
# site-specific configuration parameters that need to come before
# the platform-specific parameters - edit site.def to change

# site:  $TOG: site.sample /main/r64_final/1 1998/02/05 16:28:49 kaleb $

# site:  $XFree86: xc/config/cf/site.def,v 3.24 2000/06/25 20:17:29 dawes Exp $

# $XFree86: xc/config/cf/xf86site.def,v 3.186 2003/06/25 18:06:22 eich Exp $

# ----------------------------------------------------------------------
# platform-specific configuration parameters - edit darwin.cf to change

# $XdotOrg: xc/config/cf/darwin.cf,v 1.7 2005/03/17 01:12:07 torrey Exp $
# platform:  $XFree86: xc/config/cf/darwin.cf,v 1.50 2004/01/20 23:54:50 torrey Exp $

# operating system:  Darwin (12.6.0)

#   PROJ_BUILD = xcodebuild

#   XPLIBDIR = /usr/lib
#   XPINCDIR = $(SERVERSRC)/hw/darwin/quartz/xpr

# $XFree86: xc/config/cf/darwinLib.rules,v 1.8 2003/10/09 22:43:18 herrb Exp $

# $XdotOrg: xc/config/cf/xorg.cf,v 1.44 2005/01/27 03:50:46 ajax Exp $

# $Xorg: xfree86.cf,v 1.4 2000/08/17 19:41:49 cpqbld Exp $

# ----------------------------------------------------------------------
# site-specific configuration parameters that go after
# the platform-specific parameters - edit site.def to change

# site:  $TOG: site.sample /main/r64_final/1 1998/02/05 16:28:49 kaleb $

# site:  $XFree86: xc/config/cf/site.def,v 3.24 2000/06/25 20:17:29 dawes Exp $

# ---------------------------------------------------------------------
# Imake rules for building libraries, programs, scripts, and data files
# rules:  $Xorg: Imake.rules,v 1.3 2000/08/17 19:41:46 cpqbld Exp $
# rules:  $XdotOrg: xc/config/cf/Imake.rules,v 1.8 2005/02/01 22:27:00 ajax Exp $
#
#
#
#
# rules:  $XFree86: xc/config/cf/Imake.rules,v 3.128 2003/11/15 03:25:17 dawes Exp $

          PATHSEP = /
            SHELL = /bin/sh

              TOP = .
      CURRENT_DIR = .

            IMAKE = imake
           DEPEND = gccmakedep
        MKDIRHIER = mkdir -p
          REVPATH = revpath
    EXPORTLISTGEN =
             RMAN = RmanCmd
     RMANBASENAME = rman
      RMANOPTIONS = RmanOptions
        CONFIGSRC = $(TOP)/config
         IMAKESRC = $(CONFIGSRC)/imake
        DEPENDSRC = $(CONFIGSRC)/util

          INCROOT = /usr/local/include
        USRLIBDIR = /usr/local/lib
           VARDIR = /var
        VARLIBDIR = $(VARDIR)/lib
  SYSTEMUSRLIBDIR = /usr/lib
  SYSTEMUSRINCDIR = /usr/include
         SHLIBDIR = /usr/local/lib
     LOGDIRECTORY = $(VARDIR)/log

               AR = ar clq

.SUFFIXES: .cc

         COMPRESS = compress
          GZIPCMD = gzip

               LD = ld
               LN = ln -s
             MAKE = make
               MV = mv -f
               CP = cp

           RANLIB = ranlib

  RANLIBINSTFLAGS =

               RM = rm -f
             PERL = perl
         PERLOPTS =
     PERLENVSETUP = env LC_ALL=C

   COMPRESSMANCMD = gzip -n

            TROFF = groff -Tps
            NROFF = nroff

         MSMACROS = -ms
              TBL = tbl
              EQN = eqn
             NEQN = neqn
              COL = col
         COLFLAGS = -b


            MODLD = ld

       MODLDFLAGS =
MODLDCOMBINEFLAGS = -X -r
            MODAR = ar clq

        MODRANLIB = ranlib

     STD_INCLUDES =
#  STD_CPP_OPTIONS = -traditional -D__GNUC__
#  STD_CPP_DEFINES = -traditional -D__GNUC__ -traditional -D__GNUC__ -D__x86_64__ -D__DARWIN__                                 -DNO_ALLOCA -DCSRG_BASED  $(PROJECT_DEFINES)
#      STD_DEFINES = -D__x86_64__ -D__DARWIN__                                 -DNO_ALLOCA -DCSRG_BASED  $(PROJECT_DEFINES)
 EXTRA_LOAD_FLAGS =
  EXTRA_LDOPTIONS =
  EXTRA_LIBRARIES =
             TAGS = ctags

   PARALLELMFLAGS =

    SHAREDCODEDEF =
         SHLIBDEF =

     SHLIBLDFLAGS = -dynamiclib $(SHLIBGLOBALSFLAGS)

         NOSTDLIB = -nostdlib
     POSTNOSTDLIB = -Wl,-Bstatic -lgcc -Wl,-Bdynamic

         PICFLAGS = -dynamic

    PROTO_DEFINES =

     INSTPGMFLAGS =

     INSTBINFLAGS = -m 0755
     INSTUIDFLAGS = -m 4711
     INSTLIBFLAGS = -m 0644
     INSTINCFLAGS = -m 0444
     INSTMANFLAGS = -m 0444
     INSTDATFLAGS = -m 0444
    INSTKMEMFLAGS = -m 4711

      PROJECTROOT = /usr/local

      CDEBUGFLAGS = -Os
        CCOPTIONS = -Wall -Wpointer-arith -no-cpp-precomp

      ALLINCLUDES = $(INCLUDES) $(EXTRA_INCLUDES) $(TOP_INCLUDES) $(INSTALLED_INCLUDES) $(STD_INCLUDES)
       ALLDEFINES = $(ALLINCLUDES) $(STD_DEFINES) $(PROTO_DEFINES) $(THREADS_DEFINES) $(MODULE_DEFINES) $(DEFINES) $(EXTRA_DEFINES)
           CFLAGS = $(CDEBUGFLAGS) $(CCOPTIONS) $(THREADS_CFLAGS) $(MODULE_CFLAGS) $(ALLDEFINES)
         LDPRELIB = -L$(USRLIBDIR) $(INSTALLED_LIBS)
        LDPOSTLIB =
        LDOPTIONS = $(CDEBUGFLAGS) $(CCOPTIONS)  $(EXTRA_LDOPTIONS) $(THREADS_LDFLAGS) $(LOCAL_LDFLAGS) $(LDPRELIBS)

           LDLIBS = $(LDPOSTLIBS) $(THREADS_LIBS) $(SYS_LIBRARIES) $(EXTRA_LIBRARIES)

           CCLINK = $(CC)

     LDSTRIPFLAGS = -x
   LDCOMBINEFLAGS = -X -r
      DEPENDFLAGS =
   DEPEND_DEFINES =

        MACROFILE = darwin.cf
           RM_CMD = $(RM)

    IMAKE_DEFINES =
   IMAKE_WARNINGS = -Wundef

         IRULESRC = $(CONFIGDIR)
        IMAKE_CMD = $(IMAKE) -DUseInstalled -I$(IRULESRC) $(IMAKE_DEFINES) 		    $(IMAKE_WARNINGS)

     ICONFIGFILES = $(IRULESRC)/Imake.tmpl $(IRULESRC)/X11.tmpl 			$(IRULESRC)/site.def $(IRULESRC)/$(MACROFILE) 			$(IRULESRC)/xfree86.cf $(IRULESRC)/xf86.rules $(IRULESRC)/xorgsite.def $(IRULESRC)/host.def $(EXTRA_ICONFIGFILES)

# $Xorg: X11.rules,v 1.4 2000/08/17 19:41:46 cpqbld Exp $

# $XFree86: xc/config/cf/X11.rules,v 1.6 2001/01/17 16:22:31 dawes Exp $

# ----------------------------------------------------------------------
# X Window System Build Parameters and Rules
# $XdotOrg: xc/config/cf/X11.tmpl,v 1.52 2005/09/22 23:54:18 alanc Exp $
# $Xorg: X11.tmpl,v 1.6 2000/08/17 19:41:46 cpqbld Exp $
#
#
#
#
# $XFree86: xc/config/cf/X11.tmpl,v 1.248 2004/02/16 04:07:37 dawes Exp $

XORGRELSTRING = `echo 7 6 0 | sed -e 's/ /./g' -e 's/^/Version\\\ /'`
  XORGMANNAME = X Version 11

STICKY_DEFINES = -DHAS_STICKY_DIR_BIT

FCHOWN_DEFINES = -DHAS_FCHOWN

# -----------------------------------------------------------------------
# X Window System make variables; these need to be coordinated with rules

             XTOP = $(TOP)
           BINDIR = /usr/local/bin
     BUILDINCROOT = $(TOP)/exports
      BUILDINCDIR = $(BUILDINCROOT)/include
      BUILDINCTOP = ../..
      BUILDLIBDIR = $(TOP)/exports/lib
      BUILDLIBTOP = ../..
      BUILDBINDIR = $(TOP)/exports/bin
      BUILDBINTOP = ../..
   BUILDMODULEDIR = $(BUILDLIBDIR)/modules
     BUILDI18NDIR = $(BUILDLIBDIR)/locale
   BUILDMODULETOP = $(BUILDLIBTOP)/..
    XBUILDINCROOT = $(XTOP)/exports
     XBUILDINCDIR = $(XBUILDINCROOT)/include/X11
     XBUILDINCTOP = ../../..
     XBUILDBINDIR = $(XBUILDINCROOT)/bin
           INCDIR = $(INCROOT)
           ADMDIR = /usr/adm
           LIBDIR = /usr/X11//lib
         SHAREDIR = /usr/X11/share
       LIBEXECDIR = /usr/local/libexec
        MODULEDIR = $(USRLIBDIR)/modules
          XBINDIR = $(PROJECTROOT)/bin

       INSTSRCDIR = /usr/local/src

        ETCX11DIR = /etc/X11

          CONFDIR = $(ETCX11DIR)

#           DOCDIR = $(LIBDIR)/doc
#         DOCPSDIR = $(DOCDIR)/PostScript
#        DOCPDFDIR = $(DOCDIR)/PDF
#          FONTDIR = /opt/X11/share/fonts
#      FONTROOTDIR = /opt/X11/share/fonts
#     ENCODINGSDIR = ${FONTROOTDIR}/encodings
#         XINITDIR = $(LIBDIR)/xinit
#           XDMDIR = $(LIBDIR)/xdm
#        XDMVARDIR = $(VARLIBDIR)/xdm
#           TWMDIR = $(LIBDIR)/twm
#           XSMDIR = $(LIBDIR)/xsm
#           NLSDIR = $(LIBDIR)/nls
#       XLOCALEDIR = $(LIBDIR)/locale
#      LBXPROXYDIR = $(LIBDIR)/lbxproxy
#  PROXYMANAGERDIR = $(LIBDIR)/proxymngr
#        XPRINTDIR = $(LIBDIR)/xserver
#      XAPPLOADDIR = $(LIBDIR)/app-defaults
#       FONTCFLAGS = -t

     INSTAPPFLAGS = $(INSTDATFLAGS)

              RGB = $(XBINDIR)/rgb
            FONTC = $(XBINDIR)/bdftopcf
      MKFONTSCALE = $(XBINDIR)/mkfontscale
        MKFONTDIR = $(XBINDIR)/mkfontdir
          UCS2ANY = $(XBINDIR)/ucs2any
      BDFTRUNCATE = $(XBINDIR)/bdftruncate
     UCSMAPPREFIX = $(FONTDIR)/util/map-
       XCURSORGEN = $(XBINDIR)/xcursorgen

       DOCUTILSRC = $(XTOP)/doc/util
        CLIENTSRC = $(TOP)/clients
          DEMOSRC = $(TOP)/demos
       XDOCMACROS = $(DOCUTILSRC)/macros.t
       XIDXMACROS = $(DOCUTILSRC)/indexmacros.t
       PROGRAMSRC = $(TOP)/programs
           LIBSRC = $(XTOP)/lib
          FONTSRC = $(XTOP)/fonts
     ENCODINGSSRC = $(FONTSRC)/encodings
       INCLUDESRC = $(BUILDINCROOT)/include
      XINCLUDESRC = $(INCLUDESRC)/X11
        SERVERSRC = $(XTOP)/programs/Xserver
       CONTRIBSRC = $(XTOP)/../contrib
   UNSUPPORTEDSRC = $(XTOP)/unsupported
           DOCSRC = $(XTOP)/doc
           RGBSRC = $(XTOP)/programs/rgb
      BDFTOPCFSRC = $(PROGRAMSRC)/bdftopcf
     MKFONTDIRSRC = $(PROGRAMSRC)/mkfontdir
    FONTSERVERSRC = $(PROGRAMSRC)/xfs
       FONTINCSRC = $(XTOP)/include/fonts
        EXTINCSRC = $(XTOP)/include/extensions
      FTSOURCEDIR = FreeTypeSrcDir
        DRMSRCDIR = $(TOP)/extras/drm
       MESASRCDIR = $(TOP)/extras/Mesa
  OGLSAMPLESRCDIR = $(TOP)/extras/ogl-sample
        PSWRAPSRC = $(XTOP)/config/pswrap
     TRANSCOMMSRC = $(LIBSRC)/xtrans
   TRANS_INCLUDES = -I$(TRANSCOMMSRC)
 CONNECTION_FLAGS = -DTCPCONN -DUNIXCONN  $(STICKY_DEFINES) $(FCHOWN_DEFINES) -DIPv6
XTRANS_FAILDEFINES =

 VENDORWEBSUPPORT = http://wiki.X.Org

VENDORSUPPORTDEFS = -D__VENDORDWEBSUPPORT__='"$(VENDORWEBSUPPORT)"'

       XENVLIBDIR = $(USRLIBDIR)
   CLIENTENVSETUP = DYLD_LIBRARY_PATH=$(XENVLIBDIR) XLOCALEDIR=$(BUILDLIBDIR)/locale

# $XFree86: xc/config/cf/darwinLib.tmpl,v 1.19 2003/11/04 00:24:36 torrey Exp $

          XLIBSRC = $(LIBSRC)/X11

SOXLIBREV = 6.2
DEPXONLYLIB =
XONLYLIB =  -lX11

#LINTXONLY = $(LINTLIBDIR)/llib-lX11.ln

      DEPXLIBONLY = $(DEPXONLYLIB)
         XLIBONLY = $(XONLYLIB)
#     LINTXONLYLIB = $(LINTXONLY)

#LINTEXTENSION = $(LINTLIBDIR)/llib-lXext.ln

#LINTEXTENSIONLIB = $(LINTEXTENSION)
          DEPXLIB = $(DEPEXTENSIONLIB) $(DEPXONLYLIB)
             XLIB = $(EXTENSIONLIB) $(XONLYLIB)
#         LINTXLIB = $(LINTXONLYLIB)

FONTCONFIGINCDIR=$(INCDIR)

FONTCONFIGINCLUDES = -I$(FONTCONFIGINCDIR)

FCCACHE = $(BINDIR)/fc-cache

FONTCONFIGDEFINES = -DFONTCONFIG

XFTINCLUDES=$(FONTCONFIGINCLUDES) $(FREETYPE2INCLUDES)

LIBPNGINCDIR = /usr/include

LIBPNGINC=

LIBPNGDIR = /usr
LIBPNGLIBDIR = /usr/lib
LIBPNGINCDIR = /usr/include

LIBPNGLIB = -lpng

EXTRAXAWREQS =

EXTRAXAWCLIENTLIBS =

        CONFIGDIR = $(LIBDIR)/config

    USRLIBDIRPATH = $(USRLIBDIR)
        LDPRELIBS = -L$(USRLIBDIR)  $(INSTALLED_LIBS)
       LDPOSTLIBS =
#     TOP_INCLUDES = -I$(INCROOT) $(TOP_X_INCLUDES)
  PROJECT_DEFINES =
   VENDOR_DEFINES = -DXVENDORNAME='"$(VENDORNAME)"' -DXVENDORNAMESHORT='"$(VENDORNAMESHORT)"'

# ----------------------------------------------------------------------
# start of Imakefile

PGPLOT_LIBDIR = /opt/pgplot

PGPLOT_INCDIR = /opt/pgplot

CFITSIO_LIBDIR = /opt/local/lib

CFITSIO_INCDIR = /opt/local/include

GSL_LIBDIR = /opt/local/lib

GSL_INCDIR = /opt/local/include

LOCAL_LIBDIR = ./lib

LOCAL_INCDIR = ./include

LOCAL_LIBRARIES = $(XLIB) -L$(LIBDIR) 	-L$(PGPLOT_LIBDIR) -lcpgplot -lpgplot 	-L$(CFITSIO_LIBDIR) -lcfitsio 	-L$(GSL_LIBDIR) -lgsl -lgslcblas 	-L$(LOCAL_LIBDIR) -lastrotools -lpgtools -lmathtools -lm -lpng -lz

CCOPTIONS = -O3 -Wuninitialized -g
#CCOPTIONS = -O3 -ansi -Wuninitialized -g

CDEBUGFLAGS =

F77 = gfortran

CC = gcc

CCLINK = $(F77)

TOP_INCLUDES =  -I./include 		-I$(LOCAL_INCDIR) -I$(PGPLOT_INCDIR) -I$(CFITSIO_INCDIR)                 	-I$(GSL_INCDIR)

EXTRA_LDOPTIONS = -I./include 		-I$(LOCAL_INCDIR) -I$(PGPLOT_INCDIR) -I$(CFITSIO_INCDIR)                 	-I$(GSL_INCDIR)

SRCS =  aris.c 	baseline_number.c 	free_memory_block.c 	free_memory_block_char.c 	free_memory_block_float.c 	free_memory_block_double.c 	obs_param_file_io.c 	obs_param_set.c 	obs_param_input.c 	in__src_proc.c 	out_src_proc.c 	ch_time_set.c 	antenna_selection.c 	baseline_check.c 	number_char_cut.c 	srt_info_disp.c 	err_parameter_set.c 		 antenna_base_phase.c	 trk_priority_check.c 	tracking_init.c 	tracking_button_disp.c 	srtatt_button_disp.c 	sun_angle_check.c 	source_info_disp.c 	menu_config.c 	on_source_disp.c 	TV_menu_hatch.c 	fits_data_select.c 	text_bottom.c 	tv_button_disp.c 	comment_init.c 	comment_disp.c 	tv_get_param.c 	_off_button.c 	_on_button.c 	_toggle_button.c 	_button_chk.c 	off_button.c 	on_button.c 	toggle_button.c 	button_chk.c 	char_copy.c 	char_ncopy.c 	source_position.c 	transformation_matrices.c 	atmospheric_fluctuation.c 	Gaussian_Noise1D.c 	Gaussian_noise_check.c 	gauss_dev.c 	turbulent_phase_screen.c 	phase_screen_check.c 	antenna_visibility.c 	GRT_TEC_TID.c 	TEC_fluctuation_amp.c 	ref_pos_shift.c 	array_config.c 	position_on_screen.c 	azel_disp.c 	SSF_disp.c 	airmass.c 	nmf20.c 	nmf20_coeffi.c 	fitsidi_save.c 	brightness_disp.c 	source_model.c 	source_model3.c 	ADAF_disk.c 	VSOP2_logo.c 	jet_cjet.c 	PA_rotate.c 	ordinate_axis.c 	earth_eclipse.c 	fft2d.c 	fringe_disp.c 	fringe_fitting.c 	wfft2d.c 	get_vis.c 	visibility_calc.c 	bhs_model.c 	ccm_read.c 	mcm_read.c 	vis_smoothing.c 	input_star_position.c 	output_star_position.c 	str_init.c 	orbit_disp.c 	draw_earth.c 	data_link_schedule.c 	link_TRK_sim.c 	link_SLR_sim.c 	orbit_info_print.c 	peak_normalize.c 	spherical_geometry.c 	pg_color_bar.c 	pg_color_map.c 	pgcont_map.c 	position2uvw.c 	EPL_disp.c 	spectl_disp.c 	allanv_disp.c 	coherence_disp.c 	station_select.c 	wave_select.c 	random_val0.c 	random_val1.c 	diff.c 	imager.c 	pixel_calc.c 	phase_reference.c 	qlook_imager.c 	seed_random.c 	set_color.c 	spacecraft_position.c 	tracking_status.c 	tracking_station_name_read.c 	slr_config.c 	get_srt_link.c 	snrcal.c 	thermal_noise_cal.c 	parabo.c 	uvw_calc.c 	az_rotation.c 	az_adjustment.c 	weight_flag_add.c 	weight_flag_chk.c 	attitude_Q.c 	slew_time.c 	SEFD.c 	tracking_condition.c 	freq_std_instability.c 	asd_model.c 	modified_asd_model.c 	coherence_factor_calc.c 	uv_display.c	fit_power_two_number.c	ant_pos_err_disp.c

OBJS =  $(SRCS:.c=.o)

        PROGRAM = aris

all:: aris

aris: libmathtools libastrotools libpgtools $(OBJS)
	$(CCLINK) -o $@ $(LDOPTIONS) $(OBJS) $(LOCAL_LIBRARIES) $(LDLIBS)  $(EXTRA_LOAD_FLAGS)

libmathtools:
	cd src/mathtools; make all

libastrotools:
	cd src/astrotools; make all

libpgtools:
	cd src/pgtools; make all

install:: aris
	@if [ -d $(DESTDIR)$(BINDIR) ]; then \
		set +x; \
	else \
		if [ -h $(DESTDIR)$(BINDIR) ]; then \
			(set -x; rm -f $(DESTDIR)$(BINDIR)); \
		fi; \
		(set -x; $(MKDIRHIER) $(DESTDIR)$(BINDIR)); \
	fi
	$(INSTALL) $(INSTALLFLAGS) $(INSTPGMFLAGS)  aris $(DESTDIR)$(BINDIR)/aris

cleandir::
	cd lib; $(RM) -f *

cleandir::
	cd src/mathtools; make clean

cleandir::
	cd src/astrotools; make clean

cleandir::
	cd src/pgtools; make clean

depend::
	$(DEPEND) $(DEPENDFLAGS) -- $(ALLDEFINES) $(DEPEND_DEFINES) -- $(SRCS)

cleandir::
	$(RM) aris

# ----------------------------------------------------------------------
# common rules for all Makefiles - do not edit

.c.i:
	$(RM) $@
	 	$(CC) -E $(CFLAGS) $(_NOOP_) $*.c > $@

.SUFFIXES: .ii

.cc.ii:
	$(RM) $@
	 	$(CC) -E $(CFLAGS) $(_NOOP_) $*.cc > $@

.SUFFIXES: .s

.c.s:
	$(RM) $@
	 	$(CC) -S $(CFLAGS) $(_NOOP_) $*.c

.cc.s:
	$(RM) $@
	 	$(CC) -S $(CFLAGS) $(_NOOP_) $*.cc

emptyrule::

cleandir::
	$(RM) *.CKP *.ln *.BAK *.bak *.o core errs ,* *~ *.a .emacs_* tags TAGS make.log MakeOut   "#"*

Makefile::
	-@if [ -f Makefile ]; then set -x; \
	$(RM) Makefile.bak; $(MV) Makefile Makefile.bak; \
	else exit 0; fi
	$(IMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR)

tags::
	$(TAGS) -w *.[ch]
	$(TAGS) -xw *.[ch] > TAGS

man_keywords::

html_index::

clean:: cleandir

distclean:: cleandir

# ----------------------------------------------------------------------
# empty rules for directories that do not have SUBDIRS - do not edit

install::
	@echo "install in $(CURRENT_DIR) done"

#install.man::
#	@echo "install.man in $(CURRENT_DIR) done"

#install.sdk::
#	@echo "install.sdk in $(CURRENT_DIR) done"

Makefiles::

includes::

depend::

distclean::
	$(RM) Makefile Makefile.dep

# ----------------------------------------------------------------------
# dependencies generated by makedepend

