FLAGS=-g -std=c++11

suffix_tree_test: main.cc
	$(CXX) $(FLAGS) $? -o $@

default: all

all: suffix_tree_test

clean:
	rm suffix_tree_test

