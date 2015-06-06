// Example of using the class Pso

// Here we are optimising dimensions of a box so that neither
// of the side lengths exceeds some limit while the volume
// remains greatest possible.

//    _____________
//   /            /|
//  /____________/ |
//  |            | |  Volume = a*b*c
//  |            | |  a < limit, b < limit, c < limit
//  |            | |
//  |____________|/

#include<iostream>
#include<fstream>
#include<string>
#include"pso.h"

// Creates input file with all possible side lengthes
std::string createInputFile()
{
	std::string filename = "boxsides.txt";
	std::ofstream file(filename, std::ios::out);

	for (int j = 1; j <= 3; j++)
	{
		for (int i = 0; i < 20; i++)
			file << (i + 1) * 2 << " ";
		if (j != 3)
			file << '\n';
	}
	file.close();
	return filename;
}

// Specifies what to optimize, user-defined
class Fitness :public IFitness
{
public:
	Fitness(double limit):m_limit(limit)
	{}
	// Volume of the box
	double objective(const std::vector<double>& param)
	{
		double volume = 1;
		for (std::vector<double>::const_iterator it = param.begin(); it != param.end(); ++it)
			volume *= *it;
		return (-1)*volume;
	}
	// Each side must not exceed LIMIT
	std::set<double> constraint(const std::vector<double>& param)
	{
		std::set<double> constr;
		for (std::vector<double>::const_iterator it = param.begin(); it != param.end(); ++it)
			constr.insert(*it - m_limit);
		return constr;
	}
private:
	double m_limit;
};

template <class T> std::ostream& operator << (std::ostream& out, std::vector<T> vec)
{	// Just for convenience
	for (std::vector<T>::const_iterator it = vec.begin(); it != vec.end(); ++it)
		out << *it << ' ';
	return out;
}

int main()
{
	std::string inpFileName = createInputFile();

	Fitness fit(10);
	Pso* ppso = Pso::discrete(inpFileName.c_str(), &fit, 50);
	
	ppso->run(50, Pso::Term::fit_dev); //-1728 for LIMIT = 12, -512 for LIMIT = 8
	ppso->print_log(std::cout);

	std::cout << "\nResulting volume: " << ppso->get_best().first
		<< "\nParameters: " << ppso->get_best().second;
	
	delete ppso;
	system("pause");
	return 0;
}
