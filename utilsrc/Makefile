# DEBUG = y
MAKERULESDIR = ..
DESTDIR = ../bin.$(UNAME)
INCLUDES = -I../include 
# -I /home/debarati/softwares/arpack++/include
LIBDIR = -L../lib.$(UNAME) -L /usr/lib/

include $(MAKERULESDIR)/Makerules
# include /home/debarati/softwares/arpack++/make.inc

#VIEWSOURCES =	mesh_perceptual.cc \

VIEWSOURCES =   dynamic_mesh.cc \
#OTHERSOURCES =	mesh_info.cc \

VIEWOFILES = $(addprefix $(OBJDIR)/,$(VIEWSOURCES:.cc=.o))
OTHEROFILES = $(addprefix $(OBJDIR)/,$(OTHERSOURCES:.cc=.o))
OFILES = $(VIEWOFILES) $(OTHEROFILES)

#VIEWPROG = $(DESTDIR)/mesh_perceptual$(EXE)

VIEWPROG = $(DESTDIR)/dynamic_mesh$(EXE)
OTHERPROGS = $(addsuffix $(EXE), $(addprefix $(DESTDIR)/, $(OTHERSOURCES:.cc=)))
PROGS = $(OTHERPROGS) $(VIEWPROG)

default: $(PROGS)


LIBS += -ltrimesh 
$(VIEWPROG) : $(VIEWOFILES)
	$(LINK) $(GLLIBS)

$(OTHERPROGS) : $(DESTDIR)/%$(EXE) : $(OBJDIR)/%.o
	$(LINK)

$(PROGS) : ../lib.$(UNAME)/libtrimesh.a

clean :
	-rm -f $(OFILES) $(OBJDIR)/Makedepend $(OBJDIR)/*.d
	-rm -rf $(OBJDIR)/ii_files
	-rm -rf $(OBJDIR)/

spotless : clean
	-rm -f $(PROGS)
	-rmdir $(DESTDIR)

