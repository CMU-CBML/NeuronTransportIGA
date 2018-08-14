#include <iostream>
#include "kernel.h"
#include <sstream>
#include <string>
#include <stdio.h>
using namespace std;
int main(int argc, char **argv)
{
	if (argc == 2) {
		stringstream ss;
		string file_in;
		ss << argv[1];
		ss >> file_in;
		kernel app;
		app.run(file_in);
		cout << "DONE!\n";
		//getchar();
	}
	else if (argc > 2) {
		cout << "Too many arguments.\n";
	}
	else {
		cout << "One argument expected.\n";
	}
	return 0;
}