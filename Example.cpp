// Example of using the class Pso

// Here we are optimising dimensions of a box so that neither
// of the side lenghes exceeds some limit


//    ____________
//   /           /|
//  /___________/ |
//  |           | |  Volume = a*b*c
//  |           | |  a < limit, b < limit, c < limit
//  |           | |
//  |           | |
//  |___________|/


#include"partswarmopt.h"
#include <iostream>
#include <fstream>

string createInputFile();

int main()
{
	string inpFileName = createInputFile();
	
	vector<int> position; // for debugging

	ofstream fout("output.txt", ios::out);
	const unsigned int n_part = 10; // number of particles
	const unsigned int n_gen = 200; // number of generations

	Pso box(inpFileName, n_part);
	box.initialize();
	for (unsigned int i = 0; i < n_gen - 1; i++)
	{
		++box;
		fout << box;

		// The code below is usefull for debugging purposes only
		for (unsigned int i = 1; i <= n_part; i++)
		{
			position = box.currentPosition(i);
			fout << "Position of particle " << i << ": ";
			for (int j = 0; j < position.size(); j++)
			{
				fout << position[j] << " ";
			}
			fout << endl;
		}
		fout << "Current best position: ";
		position = box.bestPosision();
		for (unsigned int i = 0; i < 3; i++)
		{
			fout << position[i] << " ";
		}
		fout << endl << endl << endl;
	}
	// system("pause");
	return 0;
}


double Udm::Objective(const vector<double>& parameterValues)
{ // Volume of a box = a*b*c
	double ans = 1;
	for (unsigned int i = 0; i < parameterValues.size(); i++)
	{
		ans *= parameterValues[i];
	}
	return -ans;
}

vector<double> Udm::Constraints(const vector<double>& parameterValues)
{ // Box sides' lengthes are limited
	double limit = 12;
	vector<double> ans;
	double penalty;
	if (parameterValues[0] > limit)
	{
		penalty = (parameterValues[0] - limit) / (parameterValues[0] + limit);
		ans.push_back(penalty);
	}

	if (parameterValues[1] > limit)
	{
		penalty = (parameterValues[1] - limit) / (parameterValues[1] + limit);
		ans.push_back(penalty);
	}

	if (parameterValues[2] > limit)
	{
		penalty = (parameterValues[2] - limit) / (parameterValues[2] + limit);
		ans.push_back(penalty);
	}

	return ans;
}

string createInputFile()
{ // Creates file with possible side lengthes
	string filename = "boxsides.txt";
	ofstream file(filename, ios::out);

	for (int j = 1; j <= 3; j++)
	{
		file << 20 << endl;
		for (int i = 0; i < 20; i++)
		{
			file << (i + 1)*2 << " ";
		}
		file << endl;
	}
	return filename;
}