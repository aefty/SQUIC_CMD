PWD=$(shell pwd)

CXX:=g++-9
CXXFLAGS = -O0 -g -fsanitize=address


RELEASE_DIR=$(PWD)

all: SQUIC_CMD
release: $(TARGET)

SQUIC_CMD: 
	@echo "\033[92m>> Compiling SQUIC_CMD (-lSQUIC expected in /usr/local/lib ) \033[0m"
	$(shell mkdir -p $(RELEASE_DIR)) # make the folder
	$(CXX) $(CXXFLAGS) -o $(RELEASE_DIR)/$@ SQUIC_CMD.cpp -lSQUIC          # compile & link

clean: 
	@echo "\033[92m>> Cleaning \033[0m"	
	rm -rf  SQUIC_CMD