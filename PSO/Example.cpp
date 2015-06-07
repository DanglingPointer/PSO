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
using namespace Opt;
// Creates discrete input file
std::string createDInputFile()
{
	std::string filename = "discrete_input.txt";
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
// Creates continuous input file
std::string createCInputFile()
{
	std::string filename = "cont_input.txt";
	std::ofstream file(filename, std::ios::out);

	for (int j = 1; j <= 3; j++)
	{
		file << "0 40";
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
	std::string inpFileName = createCInputFile();

	Fitness fit(25.23);
	IPso* ppso = IPso::continuous(inpFileName.c_str(), &fit, 50);
	
	ppso->run(100, IPso::Term::fit_dev);
	ppso->print_log(std::cout);

	std::cout << "\nResulting volume: " << ppso->get_best().first
		<< "\nParameters: " << ppso->get_best().second;
	
	delete ppso;

	system("pause");
	return 0;
}









//class Particle
//{
//public:
//	// Doesn't set initial fitness (use init_fit())
//	Particle(int numdim, int* maxpos, Particle*& gbest) :m_ndim(numdim), m_pgbest(gbest), m_pmaxpos(maxpos)
//	{
//		m_ppos = new int[m_ndim];
//		m_pvel = new int[m_ndim];
//		m_lbest.second = new int[m_ndim];
//
//		for (int i = 0; i < m_ndim; ++i)
//		{
//			*(m_ppos + i) = rand() % (*(m_pmaxpos + i));
//			*(m_lbest.second + i) = *(m_ppos + i);
//		}
//		for (int i = 0; i < m_ndim; ++i)
//			*(m_pvel + i) = (int)(rand() % *(m_pmaxpos + i) * std::pow(-1, rand())); // max 1/1 of side length
//	}
//	~Particle()
//	{
//		delete[] m_ppos;
//		delete[] m_pvel;
//		delete[] m_lbest.second;
//	}
//	int best_pos(int dimension) const
//	{
//		return *(m_lbest.second + dimension);
//	}
//	double best_fit() const
//	{
//		return m_lbest.first;
//	}
//	Particle& init_fit(const FitConv& fit)
//	{
//		m_lbest.first = fit(m_lbest.second, m_ndim);
//		if (m_lbest.first < m_pgbest->best_fit())
//			m_pgbest = this;
//		return *this;
//	}
//	Particle& update_fit(const FitConv& fit)
//	{
//		double fitness = fit(m_ppos, m_ndim);
//		if (fitness < m_lbest.first)
//		{
//			m_lbest.first = fitness;
//			for (int i = 0; i < m_ndim; ++i)
//				*(m_lbest.second + i) = *(m_ppos + i);
//		}
//		if (m_lbest.first < m_pgbest->best_fit())
//			m_pgbest = this;
//		return *this;
//	}
//	Particle& update_pos()
//	{
//		for (int i = 0; i < m_ndim; ++i)
//			*(m_ppos + i) += *(m_pvel + i);
//
//		for (int i = 0; i < m_ndim; ++i)
//			if (*(m_ppos + i) < 0)
//			{
//				*(m_ppos + i) = 0;
//				*(m_pvel + i) *= -1;
//			}
//			else if (*(m_ppos + i) >= *(m_pmaxpos + i))
//			{
//				*(m_ppos + i) = *(m_pmaxpos + i) - 1;
//				*(m_pvel + i) *= -1;
//			}
//		return *this;
//	}
//	// w = inertia weight, X = constriction factor (prevents explosion)
//	Particle& update_vel(double w, double X, int c1, int c2)
//	{
//		double r1 = (double)rand() / RAND_MAX;
//		double r2 = (double)rand() / RAND_MAX;
//		for (int dim = 0; dim < m_ndim; ++dim)
//			*(m_pvel + dim) = (int)(X*(w*(*(m_pvel + dim)) + c1*r1*(*(m_lbest.second + dim) - *(m_ppos + dim))
//									   + c2*r2*(m_pgbest->best_pos(dim) - *(m_ppos + dim))));
//		return *this;
//	}
//private:
//	Particle*& m_pgbest;				// global best
//	std::pair<double, int*> m_lbest;	// local best fitness & position
//
//	int m_ndim;							// number of dimensions 
//	int* m_pmaxpos;						// upper bounds for position, pointer to a predefined dynamic array
//
//	int* m_ppos;						// current position
//	int* m_pvel;						// current velocity
//};