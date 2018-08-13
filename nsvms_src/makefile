include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/test
Compiler=mpiicpc
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX}  -std=c++11 -qopenmp -Wall
#-D_GLIBCXX_USE_CXX11_ABI=0
all: nsvms

nsvms: BasicDataStructure.o UserSetting.o NS_3Dsteady.o main.o 
	$(Compiler) BasicDataStructure.o UserSetting.o NS_3Dsteady.o main.o -o nsvms ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
UserSetting.o: UserSetting.cpp
	$(Compiler)	-c UserSetting.cpp $(CFLAGS)

NS_3Dsteady.o: NS_3Dsteady.cpp
	$(Compiler)	-c NS_3Dsteady.cpp $(CFLAGS)

main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)

clean:
	rm *o nsvms