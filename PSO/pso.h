#pragma once
#include<fstream>
#include<vector>
#include<exception>
#include<cmath>
#include<sstream>
#include<utility>
#include<set>

// Realization must be written by user
__interface IFitness
{
	double objective(const std::vector<double>& param);
	double constraint(const std::vector<double>& param);
};

class FeasibleSet
{
public:
	FeasibleSet(char* filename)
	{
		std::ifstream fin(filename, std::ios::in);
		if (!fin.is_open())
			throw std::invalid_argument("FeasibleSet::FeasibleSet()");

		for (m_ndim = 0; !fin.eof(); ++m_ndim) // no new line at end of the file!
			fin.ignore(10000, '\n');
		fin.seekg(0, fin.beg);

		m_psets = new std::vector<double>[m_ndim];
		std::string line;
		std::stringstream buff;
		double value;
		for (int i = 0; i < m_ndim && !fin.eof(); ++i)
		{
			std::getline(fin, line);
			buff.str(line);
			while (buff >> value)
				(m_psets + i)->push_back(value);

			line.clear;
			buff.str("");
		}
	}
	int ndim() const
	{
		return m_ndim;
	}
	int length(int dimNum) const // == maxpos[dimNum]
	{
		return (m_psets + dimNum)->size();
	}
	double value(int dimNum, int valNum) const
	{
		return (m_psets + dimNum)->at(valNum);
	}
private:
	std::vector<double>* m_psets;
	int m_ndim;
};

class FitConv
{
public:
	FitConv(IFitness* funcs, FeasibleSet* pfs) :m_pfuncs(funcs), m_pfs(pfs)
	{}
	double operator() (int* pos, int size) const
	{
		std::vector<double> values;
		for (int i = 0; i < size; ++i)
			values.push_back(m_pfs->value(i, *(pos + i)));
		return m_pfuncs->constraint(values) + m_pfuncs->constraint(values); // + possible weights etc
	}
private:
	IFitness* m_pfuncs;
	FeasibleSet* m_pfs;
};

class Particle
{
public:
	// Doesn't set initial fitness (use init_fit())
	Particle(int numdim, int* maxpos, Particle*& gbest) :m_ndim(numdim), m_pgbest(gbest), m_pmaxpos(maxpos)
	{
		m_ppos = new int[m_ndim];
		m_pvel = new int[m_ndim];

		for (int i = 0; i < m_ndim; ++i)
		{
			*(m_ppos + i) = rand() % *(m_pmaxpos + i);
			*(m_lbest.second + i) = *(m_ppos + i);
		}
		for (int i = 0; i < m_ndim; ++i)
			*(m_pvel + i) = rand() % (*(m_pmaxpos + i) / 5) * std::pow(-1, rand()); // max 1/5 of side length
	}
	~Particle()
	{
		delete[] m_ppos;
		delete[] m_pvel;
		delete[] m_lbest.second;
	}
	int best_pos(int dimension) const
	{
		return *(m_lbest.second + dimension);
	}
	double best_fit() const
	{
		return m_lbest.first;
	}
	Particle& init_fit(const FitConv& fit)
	{
		m_lbest.first = fit(m_ppos, m_ndim);
		return *this;
	}
	Particle& update_fit(const FitConv& fit)
	{
		double FitConv = fit(m_ppos, m_ndim);
		if (FitConv < m_lbest.first)
		{
			m_lbest.first = FitConv;
			for (int i = 0; i < m_ndim; ++i)
				*(m_lbest.second + i) = *(m_ppos + i);
		}
		if (m_lbest.first < m_pgbest->best_fit())
			m_pgbest = this;
		return *this;
	}
	Particle& update_pos()
	{
		for (int i = 0; i < m_ndim; ++i)
			*(m_ppos + i) += *(m_pvel + i);

		for (int i = 0; i < m_ndim; ++i)
			if (*(m_ppos + i) < 0)
			{
				*(m_ppos + i) = 0;
				*(m_pvel + i) *= -1;
			}
			else if (*(m_ppos + i) >= *(m_pmaxpos + i))
			{
				*(m_ppos + i) = *(m_pmaxpos + i) - 1;
				*(m_pvel + i) *= -1;
			}
		return *this;
	}
	// w = inertia weight, X = constriction factor
	Particle& update_vel(double w, double X, int c1, int c2)
	{
		double r1 = rand() / RAND_MAX;
		double r2 = rand() / RAND_MAX;
		for (int dim = 0; dim < m_ndim; ++dim)
			*(m_pvel + dim) = X*(w*(*(m_pvel + dim)) + c1*r1*(*(m_lbest.second + dim) - *(m_ppos + dim)) + c2*r2*(m_pgbest->best_pos(dim) - *(m_ppos + dim)));
		return *this;
	}
private:
	Particle*& m_pgbest;				// global best
	std::pair<double, int*> m_lbest;	// local best fitness & position

	int m_ndim;							// number of dimensions 
	int* m_pmaxpos;						// upper bounds for position, pointer to a predefined dynamic array

	int* m_ppos;						// current position
	int* m_pvel;						// current velocity
};

class Pso
{
public:
	Pso(char* filename, IFitness* funcs, int numPart) :m_fs(filename), m_fit(funcs, &m_fs), m_pgbest(nullptr), m_count(0)
	{
		c1 = c2 = 2;
		m_pmaxpos = new int(m_fs.ndim());
		for (int dim = 0; dim < m_fs.ndim(); ++dim)
			*(m_pmaxpos + dim) = m_fs.length(dim);

		Particle* tempptr = new Particle(m_fs.ndim(), m_pmaxpos, m_pgbest);
		m_pgbest = tempptr;
		m_parts.insert(tempptr);
		for (int n = 1; n < numPart; ++n)
			m_parts.insert(new Particle(m_fs.ndim(), m_pmaxpos, m_pgbest));
		for (std::set<Particle*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
			(*it)->init_fit(m_fit);

	}
	~Pso()
	{
		for (std::set<Particle*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
			delete *it;
		delete[] m_pmaxpos;
	}
	Pso& run(int numGen)
	{
		// ...to be written...
		for (int gen = 0; gen < numGen; ++gen)
			++m_count;
		return *this;
	}
	std::pair<double, std::vector<double>> get_best()
	{
		std::pair<double, std::vector<double>> temp;
		temp.first = m_pgbest->best_fit();
		for (int dim = 0; dim < m_fs.ndim(); ++dim)
			temp.second.push_back(m_fs.value(dim, m_pgbest->best_pos(dim)));
		return temp;
	}
	int c1;
	int c2;
private:
	double intertWeight() const
	{
		// depending on m_count
		// ...to be written...
		return 0.7;
	}
	double constrFactor() const
	{
		// ...to be written...
		return 0.5;
	}
	FeasibleSet m_fs;
	FitConv m_fit;
	std::set<Particle*> m_parts;
	Particle* m_pgbest;
	int* m_pmaxpos;
	int m_count;
};