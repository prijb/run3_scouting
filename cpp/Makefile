CC=g++
INCLUDE=-I. -I$${CMSSW_BASE}/src -I$${CMSSW_BASE}/../../../external/boost/1.75.0/include
LIBRARIES=-L$${CMSSW_BASE}/lib/$${SCRAM_ARCH} -lDataFormatsFWLite -lFWCoreUtilities -lDataFormatsCandidate -lDataFormatsCommon -lDataFormatsHepMCCandidate -lDataFormatsProvenance

SOURCES=$(wildcard tools/*.cc)
OUT=main.exe 
OBJS=main.cc rootDict.o $(SOURCES:.cc=.o)

all : $(OUT)

$(OUT): $(OBJS) run3ScoutingLooper.C
	$(CC) $(OBJS) $(INCLUDE) $(LIBRARIES) -Wno-attributes -o $(OUT) `root-config --cflags --glibs`

rootDict : LinkDef.h
	rootcint -f rootDict.cc -c headerDict.h LinkDef.h

%.o : tools/%.cc tools/%.h
	$(CC) -c $< -o $@

rootDict.o : rootDict.cc
	$(CC) -I$../include -c $< -o $@ `root-config --cflags --glibs`
	
clean:
	rm main.exe *.o tools/*.o
