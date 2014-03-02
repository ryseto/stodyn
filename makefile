objsmain = main.o system.o grid.o contable.o bond.o particle.o sbond.o sparticle.o demsystem.o
objsryuon = XA-slip.o YH.o ft.o XA.o YM-slip.o fts.o XC-slip.o YM.o gmres.o \
XC.o ZM-slip.o libiter.o XG-slip.o ZM.o lub-matrix.o XG.o bench.o\
lub.o XM-slip.o bi-cgstab.o  XM.o  bicgstab.o matrix.o\
XP.o bonds.o minv-poly.o XQ.o cg.o myblas.o YA-slip.o cg_.o non-ewald.o\
YA.o cgs.o orthomin.o YB-slip.o dgetri_c.o steepest.o YB.o\
ewald-3ft-matrix.o stokes.o YC-slip.o ewald-3ft.o \
YC.o ewald-3fts-matrix.o two-body-res.o YG-slip.o ewald-3fts.o\
twobody-slip.o YG.o ewald.o twobody.o YH-slip.o f.o 
objs = $(objsmain) $(objsryuon)

CC = g++
LFLAGS = -L/usr/lib/atlas -llapack -lf77blas -lcblas -latlas -lgfortran
CXXFLAGS = -O3 

all : stodyn 

stodyn: $(objsmain) $(objsryuon)
	$(CC) -o /home/seto/bin/stodyn $(CXXFLAGS) $(objs) /usr/lib/lapack/lapack_DellHPC.a  $(LFLAGS)

clean: 
	@rm -rf $(objsmain) /home/seto/bin/stodyn
cleanall: 
	@rm -rf $(objs) /home/seto/bin/stodyn

main.o : main.cpp
	$(CC) -c $(CXXFLAGS) $<
system.o : system.cpp 
	$(CC) -c $(CXXFLAGS) $<
grid.o : grid.cpp
	$(CC) -c $(CXXFLAGS) $<
contable.o : contable.cpp 
	$(CC) -c $(CXXFLAGS) $<
bond.o : bond.cpp
	$(CC) -c $(CXXFLAGS) $<
particle.o : particle.cpp
	$(CC) -c $(CXXFLAGS) $<
sbond.o : sbond.cpp
	$(CC) -c  $(CXXFLAGS) $<
sparticle.o : sparticle.cpp
	$(CC) -c $(CXXFLAGS) $<
bonds.o : bonds.cpp
	$(CC) -c  $(CXXFLAGS) $<
demsystem.o : demsystem.cpp
	$(CC) -c  $(CXXFLAGS) $<
XA-slip.o : XA-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
YG.o : YG.cpp
	$(CC) -c $(CXXFLAGS) $<
YH.o : YH.cpp
	$(CC) -c $(CXXFLAGS) $<
ft.o : ft.cpp
	$(CC) -c $(CXXFLAGS) $<
XA.o : XA.cpp
	$(CC) -c $(CXXFLAGS) $<
YM-slip.o : YM-slip.cpp 
	$(CC) -c $(CXXFLAGS) $<
fts.o : fts.cpp 
	$(CC) -c $(CXXFLAGS) $<
XC-slip.o : XC-slip.cpp 
	$(CC) -c $(CXXFLAGS) $<
YM.o : YM.cpp 
	$(CC) -c $(CXXFLAGS) $<
gmres.o : gmres.cpp 
	$(CC) -c $(CXXFLAGS) $< 
XC.o : XC.cpp 
	$(CC) -c $(CXXFLAGS) $<
ZM-slip.o : ZM-slip.cpp 
	$(CC) -c $(CXXFLAGS) $<
libiter.o : libiter.cpp 
	$(CC) -c $(CXXFLAGS) $<
XG-slip.o : XG-slip.cpp 
	$(CC) -c $(CXXFLAGS) $<
ZM.o : ZM.cpp
	$(CC) -c $(CXXFLAGS) $<
lub-matrix.o : lub-matrix.cpp
	$(CC) -c $(CXXFLAGS) $<
XG.o : XG.cpp
	$(CC) -c $(CXXFLAGS) $<
bench.o : bench.cpp
	$(CC) -c $(CXXFLAGS) $<
lub.o : lub.cpp
	$(CC) -c $(CXXFLAGS) $<
XM-slip.o : XM-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
bi-cgstab.o : bi-cgstab.cpp
	$(CC) -c $(CXXFLAGS) $<
XM.o : XM.cpp
	$(CC) -c $(CXXFLAGS) $<
bicgstab.o : bicgstab.cpp
	$(CC) -c $(CXXFLAGS) $<
matrix.o : matrix.cpp
	$(CC) -c $(CXXFLAGS) $<
XP.o : XP.cpp
	$(CC) -c $(CXXFLAGS) $<
minv-poly.o : minv-poly.cpp
	$(CC) -c $(CXXFLAGS) $<
XQ.o : XQ.cpp 
	$(CC) -c $(CXXFLAGS) $<
cg.o : cg.cpp 
	$(CC) -c $(CXXFLAGS) $<
myblas.o : myblas.cpp
	$(CC) -c $(CXXFLAGS) $<
YA-slip.o : YA-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
cg_.o : cg_.cpp
	$(CC) -c $(CXXFLAGS) $<
non-ewald.o : non-ewald.cpp
	$(CC) -c $(CXXFLAGS) $<
YA.o : YA.cpp 
	$(CC) -c $(CXXFLAGS) $<
cgs.o : cgs.cpp
	$(CC) -c $(CXXFLAGS) $<
orthomin.o : orthomin.cpp
	$(CC) -c $(CXXFLAGS) $<
YB-slip.o : YB-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
dgetri_c.o : dgetri_c.cpp
	$(CC) -c $(CXXFLAGS) $<
steepest.o : steepest.cpp
	$(CC) -c $(CXXFLAGS) $<
YB.o : YB.cpp
	$(CC) -c $(CXXFLAGS) $<
ewald-3ft-matrix.o : ewald-3ft-matrix.cpp
	$(CC) -c $(CXXFLAGS) $<
stokes.o : stokes.cpp
	$(CC) -c $(CXXFLAGS) $<
YC-slip.o : YC-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
ewald-3ft.o : ewald-3ft.cpp
	$(CC) -c $(CXXFLAGS) $<
YC.o : YC.cpp
	$(CC) -c $(CXXFLAGS) $<
ewald-3fts-matrix.o : ewald-3fts-matrix.cpp
	$(CC) -c $(CXXFLAGS) $<
two-body-res.o : two-body-res.cpp 
	$(CC) -c $(CXXFLAGS) $<
YG-slip.o : YG-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
ewald-3fts.o : ewald-3fts.cpp
	$(CC) -c $(CXXFLAGS) $<
twobody-slip.o : twobody-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
ewald.o : ewald.cpp
	$(CC) -c $(CXXFLAGS) $<
twobody.o : twobody.cpp 
	$(CC) -c $(CXXFLAGS) $<
YH-slip.o : YH-slip.cpp
	$(CC) -c $(CXXFLAGS) $<
f.o : f.cpp 
	$(CC) -c $(CXXFLAGS) $<
