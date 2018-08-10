TARGET = miev
FSOURCES = miev0_py.f ErrPack.f
SKIP = skip: LPCOEF LPCO1T LPCO2T BIGA MIPRNT SMALL1 SMALL2 TESTMI ErrMsg
$(TARGET).so:
	f2py2 -c -m $(TARGET)for $(FSOURCES) $(SKIP)
	f2py -c -m $(TARGET)for $(FSOURCES) $(SKIP)
$(TARGET):
	gfortran RDI1MACHmodule.F90 MVTstNew.f $(FSOURCES) -o $(TARGET)
install: 
	cp $(TARGET)for.*.so $(QSOFT)/python_modules
	cp $(TARGET)for.so $(QSOFT)/python_modules
	cp $(TARGET).py $(QSOFT)/python_modules
