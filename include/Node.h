#pragma once

#include <vector>
#include <string>

// NODE CLASS
class Node
{
public:

	long id;
	static long num_nodes;
	std::vector<double> coords;
	long left;
	long right;
	int boundary;
	std::vector<long> elems;

public:
	Node(const std::vector<std::string>&);
	Node(double, double, double);
	~Node();

	void print_coords(void);
	void print_elements(void);
	
private:

};
