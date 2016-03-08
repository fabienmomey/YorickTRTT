# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/home/momey/apps/easy-yorick/arch/linux-x86_64/yorick
Y_EXE=/home/momey/apps/easy-yorick/arch/linux-x86_64/yorick/bin/yorick
Y_EXE_PKGS=
Y_EXE_HOME=/home/momey/apps/easy-yorick/arch/linux-x86_64/yorick
Y_EXE_SITE=/home/momey/apps/easy-yorick/arch/linux-x86_64/yorick
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

PKG_NAME=Ytrtt
PKG_I=Ytrtt_plugin.i

OBJS=Ygauss.o trtt_common.o trtt_xform3d.o trtt_2D_projectors.o trtt_2D_projectors_cubic.o trtt_3D_projectors.o trtt_3D_projectors_cubic.o Ytrtt.o

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS=-g
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c
Ygauss.o: gauss.c
trtt_common.o: trtt_common.c trtt_common.h
trtt_xform3d.o: trtt_xform3d.c trtt_xform3d.h trtt_common.h gauss.c
trtt_2D_projectors.o: trtt_2D_projectors.c trtt_2D_projectors.h trtt_xform3d.h trtt_common.h
trtt_3D_projectors.o: trtt_3D_projectors.c trtt_3D_projectors.h trtt_xform3d.h trtt_common.h
trtt_2D_projectors_cubic.o: trtt_2D_projectors_cubic.c trtt_2D_projectors.h trtt_xform3d.h trtt_common.h
trtt_3D_projectors_cubic.o: trtt_3D_projectors_cubic.c trtt_3D_projectors_cubic.h trtt_xform3d.h trtt_common.h
Ytrtt.o: Ytrtt.c trtt_xform3d.h trtt_2D_projectors.h trtt_3D_projectors.h trtt_common.h
# -------------------------------------------------------- end of Makefile
