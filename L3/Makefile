FC = gfortran
FF = -g -O3  -fbounds-check -fbacktrace

NCDF_ROOT = /usr

LIB = -L$(NCDF_ROOT)/lib -lnetcdf -lnetcdff
OBJ = mod_write_nc.o
NAME = st_venant.f90
MODNAME = src/mod_write_nc.f90

All: st_venant.x

st_venant.x: src/$(NAME)
	$(FC) $(FF) -o st_venant.x ${MODNAME} src/$(NAME) -I $(NCDF_ROOT)/include $(LIB)
hej:
	@echo hej!
clean:
	rm -f *.x *.bin src/*.o src/*.mod *~ *\#
