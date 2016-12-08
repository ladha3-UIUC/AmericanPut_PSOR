#include <iostream>
#include <fstream>
#include "AmericanPut.h"

int main(void) {
	American P(81, 200, 0.02, 0.2, 20, 50, 20, 1.3);

	double ** result = NULL;
	result = P.putPrice();

	// write result to a .csv file
	std::ofstream output_file;
	output_file.open("aPut.csv");
	for (int i = 0; i < P.M + 1; i++) {
		for (int j = 0; j < P.N - 1; j++) {
			output_file << result[i][j] << ',';
		}
		output_file << result[i][P.N - 1] << "\n";
	}

	std::cout << "Done!" << std::endl;
}
