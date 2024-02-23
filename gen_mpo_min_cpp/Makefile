# 1. Put this file in the same folder as your 'driver' code 
#    (the code containing the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it)
LIBRARY_DIR = /home/mankai/Documents/2024_documents/research/fuzzy_sphere_dmrg/itensor

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP = gen_mpo#measure_OES#fourLayer/Excited4#Excited##defect_correlator#RES_GENERAL#RES_QHFM#Laughlin#MUTUAL_QHFM#RES_IQH#measure_ope#measure_H2

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS = interaction/cg.h site/fermionBilayer.h interaction/operator.h
APPHEADERS = $(HEADERS)
#APPHEADERS += $(HEADERS)

# 5. For any additional .cc (source) files making up your project,
#    add their full filenames here.
CCFILES  = interaction/operator.cc

SINGLEFILES = interaction/cg.cc 

#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h 
CCFLAGS += -I/usr/local/include/eigen3 -I/${BOOST_INCLUDE}
CCGFLAGS += -I/usr/local/include/eigen3 -I/${BOOST_INCLUDE}

#Mappings --------------

OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %.cc,%-g.o, $(CCFILES))
#GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS) $(SINGLEOBJECTS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

#.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
%-g.o: %.cc $(HEADERS) $(TENSOR_HEADERS) $(SINGLEGOBJECTS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<


#Mappings --------------

SINGLEOBJECTS=$(patsubst %.cc,%.o, $(SINGLEFILES))
SINGLEGOBJECTS=$(patsubst %.cc,%-g.o, $(SINGLEOBJECTS))
#SINGLEGOBJECTS=$(patsubst %,.debug_objs/%, $(SINGLEOBJECTS))

#Rules ------------------

%.o: %.cc %.h
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

%-g.o: %.cc %.h
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

OBJECTS += $(SINGLEOBJECTS) $(APP).o#envop/localmpoFQHE.o fqhe/fqhe_sl.o $(APP).o
GOBJECTS+= $(SINGLEGOBJECTS) $(APP)-g.o#envop/localmpoFQHE-g.o fqhe/fqhe_sl-g.o $(APP)-g.o

build: bin/$(APP).exe
debug: bin/$(APP)-g.exe

bin/$(APP).exe: makebin $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o bin/$(APP).exe $(LIBFLAGS)

bin/$(APP)-g.exe: makebin mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o bin/$(APP)-g.exe $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o *.exe bin/$(APP).exe bin/$(APP)-g
	rm -fr $(OBJECTS) $(GOBJECTS)
	rm -fr bin

mkdebugdir:
	mkdir -p .debug_objs

makebin:
	mkdir -p bin bin/fourLayer

#APP -------------------
envop/localmpoFQHE.o: envop/localmpoFQHE.cc $(APPHEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<
envop/localmpoFQHE-g.o: envop/localmpoFQHE.cc $(APPHEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<
	
$(APP).o:  $(APP).cc $(APPHEADERS)  $(SINGLEHEADERS)#solver/dmrgLabour.h solver/linearSolver.h res/res.h
	$(CCCOM) -c $(CCFLAGS) -o $@ $<
$(APP)-g.o:  $(APP).cc $(APPHEADERS) $(SINGLEHEADERS)#solver/dmrgLabour.h solver/linearSolver.h res/res.h
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

fqhe/fqhe_sl.o: fqhe/fqhe_sl.cc $(APPHEADERS) #solver/dmrgLabour.h solver/linearSolver.h envop/TransferMatrix.h
	$(CCCOM) -c $(CCFLAGS) -o $@ $< 
fqhe/fqhe_sl-g.o: fqhe/fqhe_sl.cc $(APPHEADERS) #solver/dmrgLabour.h solver/linearSolver.h envop/TransferMatrix.h
	$(CCCOM) -c $(CCGFLAGS) -o $@ $< 

#test/test: test/test.cc $(ITENSOR_LIBS) 
#	$(CCCOM) $(CCFLAGS) test/test.cc -o test/test $(LIBFLAGS)
#test/test-g: test/test.cc $(ITENSOR_LIBS) 
#	$(CCCOM) $(CCFLAGS) test/test.cc -o test/test $(LIBFLAGS)
