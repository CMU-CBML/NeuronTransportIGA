# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
Compiler=icpc
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=  -c -std=c++11 -Wall 
#-DMKL_ILP64
FLAGS_LINK = -qopenmp -lpthread -lm -ldl


all: spline

spline: BasicDataStructure.o kernel.o KnotInsertion.o  main.o 
	$(Compiler) BasicDataStructure.o kernel.o KnotInsertion.o  main.o -o spline $(FLAGS_LINK)

BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	$(CFLAGS) BasicDataStructure.cpp $(FLAGS_LINK)
	
kernel.o: kernel.cpp
	$(Compiler)	$(CFLAGS) kernel.cpp $(FLAGS_LINK)

KnotInsertion.o: KnotInsertion.cpp
	$(Compiler)	$(CFLAGS) KnotInsertion.cpp $(FLAGS_LINK)

main.o: main.cpp
	$(Compiler)	$(CFLAGS) main.cpp $(FLAGS_LINK)

clean:
	rm *o

