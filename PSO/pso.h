#pragma once

// Particle swarm optimization for discrete parameter values
// defined by user and imported through a text file. All values
// for one parameter must be on the same line. No new line at
// the end of the file.
// Objective and constraint functions are to be implemented
// by deriving from IFitness. Otherwise only the class Pso
// is supposed to be used directly by user. 

#include<fstream>
#include<vector>
#include<exception>
#include<cmath>
#include<sstream>
#include<string>
#include<utility>
#include<set>
#include<list>
#include<iomanip>
#include<algorithm>

// Realization must be written by user
__interface IFitness
{
	double objective(const std::vector<double>& param);
	std::set<double> constraint(const std::vector<double>& param);
};
class ParamSet
{
public:
	ParamSet(char* filename)
	{
		std::ifstream fin(filename, std::ios::in);
		if (!fin.is_open())
			throw std::invalid_argument("ParamSet::ParamSet()");

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
			line.clear();
		}
		fin.close();
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
	FitConv(IFitness* funcs, ParamSet* pfs, int* count) :m_pfuncs(funcs), m_pfs(pfs), m_pcount(count)
	{}
	double operator() (int* pos, int size) const
	{
		std::vector<double> values;
		for (int i = 0; i < size; ++i)
			values.push_back(m_pfs->value(i, *(pos + i)));
		
		double alpha = rand() % 1900 + 100;
		double H = 0;
		std::set<double> *p = &(m_pfuncs->constraint(values));
		for (std::set<double>::const_iterator it = p->begin(); it != p->end(); ++it)
		{
			double gamma = (*it < 1) ? 1.0 : 2.0;
			double theta;
			if (*it < 0.00001)		theta = alpha;
			else if (*it < 0.001)	theta = 10 * alpha;
			else if (*it < 1)		theta = 100 * alpha;
			else					theta = 1000 * alpha;
			H += (theta * std::pow(*it, gamma));
		}
		return m_pfuncs->objective(values) * std::pow(*m_pcount, 0.1) * H;
	}
private:
	IFitness* m_pfuncs;
	ParamSet* m_pfs;
	int* m_pcount;
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
			*(m_pvel + i) = (int)(rand() % (*(m_pmaxpos + i) / 5) * std::pow(-1, rand())); // max 1/5 of side length
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
		if (m_lbest.first < m_pgbest->best_fit())
			m_pgbest = this;
		return *this;
	}
	Particle& update_fit(const FitConv& fit)
	{
		double fitness = fit(m_ppos, m_ndim);
		if (fitness < m_lbest.first)
		{
			m_lbest.first = fitness;
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
	// w = inertia weight, X = constriction factor (prevents explosion)
	Particle& update_vel(double w, double X, int c1, int c2)
	{
		double r1 = rand() / RAND_MAX;
		double r2 = rand() / RAND_MAX;
		for (int dim = 0; dim < m_ndim; ++dim)
			*(m_pvel + dim) = (int)( X*(w*(*(m_pvel + dim)) + c1*r1*(*(m_lbest.second + dim) - *(m_ppos + dim)) 
									   + c2*r2*(m_pgbest->best_pos(dim) - *(m_ppos + dim))) );
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
	Pso(char* filename, IFitness* funcs, int numPart) :m_fs(filename), m_pgbest(nullptr), m_count(0), m_fit(funcs, &m_fs, &m_count)
	{
		c1 = c2 = 2;
		m_pmaxpos = new int(m_fs.ndim());
		for (int dim = 0; dim < m_fs.ndim(); ++dim)
			*(m_pmaxpos + dim) = m_fs.length(dim);

		m_pgbest = new Particle(m_fs.ndim(), m_pmaxpos, m_pgbest);
		m_parts.insert(m_pgbest);
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
	enum class Term	// termination criterions
	{
		gen_count,	// generation count 
		fit_dev		// no gbest deviation for 'numGen' iterations
	};
	Pso& run(int numGen, Term criterion = Term::gen_count)
	{
		m_log.push_back(m_pgbest->best_fit());
		double X = constriction();
		double w;
		if (criterion == Term::gen_count)
		{
			for (int gen = 0; gen < numGen; ++gen)
			{
				w = chaotic_w(numGen);
				for (std::set<Particle*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
					(*it)->update_pos();
				for (std::set<Particle*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
				{
					(*it)->update_vel(w, X, c1, c2);
					(*it)->update_fit(m_fit);
				}
				++m_count;
				m_log.push_back(m_pgbest->best_fit());
			}
		}
		else if (criterion == Term::fit_dev)
		{
			double prev_fit = m_pgbest->best_fit();
			double curnt_fit;
			for (int count = 0; count < numGen; ++count)
			{
				w = rand_w();
				for (std::set<Particle*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
					(*it)->update_pos();
				for (std::set<Particle*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
				{
					(*it)->update_vel(w, X, c1, c2);
					(*it)->update_fit(m_fit);
				}
				curnt_fit = m_pgbest->best_fit();
				if (prev_fit != curnt_fit)
					count = 0;
				prev_fit = curnt_fit;
				++m_count;
				m_log.push_back(m_pgbest->best_fit());
			}
		}
		return *this;
	}
	Pso& reset_counter()
	{
		m_count = 0;
		return *this;
	}
	std::pair<double, std::vector<double>> get_best() const
	{
		std::pair<double, std::vector<double>> temp;
		temp.first = m_pgbest->best_fit();
		for (int dim = 0; dim < m_fs.ndim(); ++dim)
			temp.second.push_back(m_fs.value(dim, m_pgbest->best_pos(dim)));
		return temp;
	}
	void print_log(std::ostream& out) const
	{
		out << "Generation" << " || " << "Best fitness" << '\n';
		out << "---------- || ------------\n";
		int gen = 1;
		for (std::list<double>::const_iterator it = m_log.begin(); it != m_log.end(); ++it, ++gen)
			out << std::setw(6) << gen << "    " << " || " << "    " << *it << '\n';
		out << "---------- || ------------\n\n";
	}
	int c1;
	int c2;
private:
	double rand_w() const
	{
		return 0.5 + 0.5*(rand() / RAND_MAX);
	}
	double chaotic_w(int maxGen) const
	{
		double z = rand() / RAND_MAX;
		z = 4 * z*(1 - z);
		return (0.9 - 0.4)*(maxGen - m_count) / maxGen + 0.4*z;
	}
	double constriction() const
	{
		double phi = std::max(c1 + c2, 4);
		return 2 / std::abs(2 - phi - std::sqrt(phi*phi - 4 * phi));
	}
	ParamSet m_fs;
	FitConv m_fit;
	std::set<Particle*> m_parts;
	Particle* m_pgbest;
	int* m_pmaxpos;
	int m_count;
	std::list<double> m_log;
};