// Example of using the class Pso

#include"pso.h"
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<vector>
#include<set>
#include<cmath>
using namespace Opt;
using std::pow;

// Benchmark problems 1 and 2
class Benchmark :public IFitness
{
public:
	explicit Benchmark(int n) :m_what(n)
	{
		std::stringstream ss;
		ss << "benchmark" << m_what << "_output.txt";
		outstream = std::ofstream(ss.str());
		outstream << "Results obtained using Benchmark problem " << m_what << "\n\n";
	}
	~Benchmark()
	{
		outstream.close();
	}
	double objective(const std::vector<double>& p)
	{
		if (m_what == 1)
		{
			return pow(p.at(0) - 10, 2) + 5 * pow(p.at(1) - 12, 2) + pow(p.at(2), 4) +
				3 * pow(p.at(3) - 11, 2) + 10 * pow(p.at(4), 6) + 7 * pow(p.at(5), 2) +
				pow(p.at(6), 4) - 4 * p.at(5)*p.at(6) - 10 * p.at(5) - 8 * p.at(6);
		}
		else if (m_what == 2)
		{
			return -(-2 * p.at(0)*p.at(0) + 2 * p.at(0)*p.at(1) - 2 * p.at(1)*p.at(1) + 
					 4 * p.at(0) + 6 * p.at(1));
		}
	}
	std::set<double> constraint(const std::vector<double>& p)
	{
		std::set<double>temp;
		if (m_what == 1)
		{
			double c[4];
			c[0] = -127 + 2 * p.at(0)*p.at(0) + 3 * pow(p.at(1), 4) + p.at(2) +
				4 * p.at(3)*p.at(3) + 5 * p.at(4);
			c[1] = -282 + 7 * p.at(0) + 3 * p.at(1) + 10 * p.at(2)*p.at(2) +
				p.at(3) - p.at(4);
			c[2] = -196 + 23 * p.at(0) + p.at(1)*p.at(1) + 6 * p.at(5)*p.at(5) -
				8 * p.at(6);
			c[3] = 4 * p.at(0)*p.at(0) + p.at(1)*p.at(1) - 3 * p.at(0)*p.at(1) +
				2 * p.at(2)*p.at(2) + 5 * p.at(5) - 11 * p.at(6);
			for (int i = 0; i < 4; ++i)
				temp.insert(c[i]);
		}
		else if (m_what == 2)
		{
			temp.insert(p.at(0) + p.at(1) - 2);
			temp.insert(p.at(0) + 5 * p.at(1) - 5);
		}
		return temp;
	}
	std::string input() const
	{
		std::stringstream ss;
		ss << "benchmark" << m_what << "_input.txt";
		std::string filename = ss.str();
		std::ofstream fout(filename, std::ios::out);
		if (m_what == 1)
		{
			for (int line = 1; line <= 7; ++line)
			{
				fout << "-10 10";
				if (line != 7)	fout << '\n';
			}
		}
		else if (m_what == 2)
		{
			fout << "0 10\n0 10";
		}
		fout.close();
		return filename;
	}
	std::ofstream outstream;
private:
	int m_what;
};

template <class T> std::ostream& operator << (std::ostream& out, std::vector<T> vec)
{	// Just for convenience
	for (std::vector<T>::const_iterator it = vec.begin(); it != vec.end(); ++it)
		out << *it << ' ';
	return out;
}

int main()
{
	Benchmark bm(2);

	Pso* ppso = Pso::continuous(bm.input().c_str(), &bm, 200);
	ppso->run(3000, Pso::Term::gen_count);

	ppso->print_log(bm.outstream);

	bm.outstream << "\nResulting fitness value: " << ppso->get_best().first
		<< "\nParameters: " << ppso->get_best().second << '\n';

	delete ppso;

	system("pause");
	return 0;
}