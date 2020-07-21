PWD=$(shell pwd)

CXX?=g++
CXXFLAGS =-O3  #-O0 -g -fsanitize=address

RELEASE_DIR=$(PWD)
SQUIC_DIR=$(shell echo $$HOME)

all: SQUIC_CMD
release: $(TARGET)

SQUIC_CMD: 
	@echo -e "\033[92m>> Compiling SQUIC_CMD (-lSQUIC expected in ${SQUIC_DIR} ) \033[0m"
	$(shell mkdir -p $(RELEASE_DIR)) # make the folder
	$(CXX) $(CXXFLAGS) -o $(RELEASE_DIR)/$@ SQUIC_CMD.cpp  -static-libgfortran -static-libstdc++ -static-libgcc  -L${SQUIC_DIR} -lSQUIC  -Wl,-rpath=${SQUIC_DIR} #compile & link
clean: 
	@echo -e "\033[92m>> Cleaning \033[0m"	
	rm -rf  SQUIC_CMD