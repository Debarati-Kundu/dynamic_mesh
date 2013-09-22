include ../../../lapack-3.4.2/make.inc

# DEBUG = y
LOADER = gfortran
MAKERULESDIR = ..
DESTDIR = ../bin.$(UNAME)
INCLUDES = -I../include -I../../../lapack-3.4.2/lapacke/include/ 
#-I /home/debarati/softwares/arpack++/include
LIBDIR = -L../lib.$(UNAME) -L /usr/lib/

include $(MAKERULESDIR)/Makerules
# include /home/debarati/softwares/arpack++/make.inc

#VIEWSOURCES =	mesh_perceptual.cc \

VIEWSOURCES =  mesh_view.cc \
#dynamic_mesh.cc \
#OTHERSOURCES =	mesh_view.cc \

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

svd :   basis_SVD.o ../../../lapack-3.4.2/lib/liblapack.a ../../../lapack-3.4.2/liblapacke.a
	$(LOADER) basis_SVD.o \
	../../../lapack-3.4.2/liblapacke.a ../../../lapack-3.4.2/lib/liblapack.a ../../../lapack-3.4.2/lib/libblas.a -o $@

corsvd :  cor_SVD.o ../../../lapack-3.4.2/lib/liblapack.a ../../../lapack-3.4.2/liblapacke.a
	$(LOADER) cor_SVD.o \
	../../../lapack-3.4.2/liblapacke.a ../../../lapack-3.4.2/lib/liblapack.a ../../../lapack-3.4.2/lib/libblas.a -o $@

clean :
	-rm -f $(OFILES) $(OBJDIR)/Makedepend $(OBJDIR)/*.d
	-rm -rf $(OBJDIR)/ii_files
	-rm -rf $(OBJDIR)/
	-rm *.o
	-rm svd

spotless : clean
	-rm -f $(PROGS)
	-rmdir $(DESTDIR)

