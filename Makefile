sampleobjects = buffer_manager.o file_manager.o r_tree.o

r_tree : $(sampleobjects)
	     g++ -std=c++11 -o r_tree $(sampleobjects)

r_tree.o : r_tree.cpp
	g++ -std=c++11 -c r_tree.cpp

buffer_manager.o : buffer_manager.cpp
	g++ -std=c++11 -c buffer_manager.cpp

file_manager.o : file_manager.cpp
	g++ -std=c++11 -c file_manager.cpp

clean :
	rm -f *.o
	rm -f r_tree
