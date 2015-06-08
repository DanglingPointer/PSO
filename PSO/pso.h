#pragma once

// Particle swarm optimization for discrete or continuous
// parameter values imported from a text file.
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

namespace Opt
{
	// ===============USER==INTERFACE============================================
	
	// Realization must be written by user (what to optimize)
	__interface IFitness
	{
		// Function to minimize
		double objective(const std::vector<double>& param);
		// Feasible when (-inf, 0], infeasible when (0, +inf)
		std::set<double> constraint(const std::vector<double>& param);
	};

	class Pso
	{
	public:
		// Create new discrete PSO
		static Pso* discrete(const char* filename, IFitness* funcs, int numPart);
		// Create new continuous PSO
		static Pso* continuous(const char* filename, IFitness* funcs, int numParts);

		virtual ~Pso()
		{}
		// Termination criterions
		enum class Term
		{
			gen_count,	// generation count
			fit_dev		// no deviation in gbest during 'numGen' iterations
		};
		// Carry out optimization
		virtual Pso* run(int numGen, Pso::Term criterion) = 0;
		// Nullify iteration counter
		virtual Pso* reset_counter() = 0;
		// Optimization results
		virtual std::pair<double, std::vector<double>> get_best() const = 0;
		// Best fitness history
		virtual void print_log(std::ostream& out) const = 0;
		// Local learning factor
		int c1;
		// Global learning factor
		int c2;
	};

	// ======================GENERIC=============================================

	template<class T> class FitConv
	{
	public:
		virtual ~FitConv<T>(){}
		virtual double operator () (T* pos, int size) const = 0;
	};
	template <class T> class Particle
	{
	protected:
		// Doesn't set initial fitness (use init_fit())
		Particle<T>(int numdim, Particle<T>*& gbest) : m_ndim(numdim), m_pgbest(gbest)
		{}
	public:
		virtual ~Particle<T>()
		{
			delete[] m_ppos;
			delete[] m_pvel;
			delete[] m_lbest.second;
		}
		T best_pos(int dimension) const
		{
			return *(m_lbest.second + dimension);
		}
		double best_fit() const
		{
			return m_lbest.first;
		}
		void init_fit(const FitConv<T>& fit)
		{
			m_lbest.first = fit(m_lbest.second, m_ndim);
			if (m_lbest.first < m_pgbest->best_fit())
				m_pgbest = this;
		}
		void update_fit(const FitConv<T>& fit)
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
		}
		virtual void update_pos() = 0;
		// w = inertia weight, X = constriction factor
		void update_vel(double w, double X, int c1, int c2)
		{
			double r1 = (double)rand() / RAND_MAX;
			double r2 = (double)rand() / RAND_MAX;
			for (int dim = 0; dim < m_ndim; ++dim)
				*(m_pvel + dim) = (T)(X*(w*(*(m_pvel + dim)) + c1*r1*(*(m_lbest.second + dim) - *(m_ppos + dim))
										 + c2*r2*(m_pgbest->best_pos(dim) - *(m_ppos + dim))));
		}
	protected:
		Particle*& m_pgbest;			// global best
		std::pair<double, T*> m_lbest;	// local best fitness & position
		int m_ndim;						// number of dimensions 
		T* m_ppos;						// current position
		T* m_pvel;						// current velocity
	};
	template<class Val_t, class Par_t> class Pso_base :public Pso // Val_t = int or double, Par_t = ParamSet or ParamBounds
	{
	protected:
		Pso_base(const char* filename, IFitness* funcs, int numPart) :m_count(0), m_fs(filename), m_pgbest(nullptr), m_pfit(nullptr)
		{
			c1 = c2 = 2;
		}
	public:
		virtual ~Pso_base()
		{
			for (std::set<Particle<Val_t>*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
				delete *it;
			delete m_pfit;
		}
		void print_log(std::ostream& out) const
		{
			out << "\nGeneration" << " || " << "Best fitness" << '\n';
			out << "---------- || ------------\n";
			int gen = 1;
			for (std::list<double>::const_iterator it = m_log.begin(); it != m_log.end(); ++it, ++gen)
				out << std::setw(6) << gen << "    " << " || " << "    " << *it << '\n';
			out << "---------- || ------------\n\n";
		}
	protected:
		void run_generic(int numGen, Pso::Term criterion)
		{
			m_log.push_back(m_pgbest->best_fit());
			double X = constriction();
			double w;
			if (criterion == Term::gen_count)
			{
				for (int gen = 0; gen < numGen; ++gen)
				{
					w = chaotic_w(numGen);
					for (std::set<Particle<Val_t>*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
						(*it)->update_pos();
					for (std::set<Particle<Val_t>*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
					{
						(*it)->update_vel(w, X, c1, c2);
						(*it)->update_fit(*m_pfit);
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
					w = const_w();
					for (std::set<Particle<Val_t>*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
						(*it)->update_pos();
					for (std::set<Particle<Val_t>*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
					{
						(*it)->update_vel(w, X, c1, c2);
						(*it)->update_fit(*m_pfit);
					}
					curnt_fit = m_pgbest->best_fit();
					if (prev_fit != curnt_fit)
						count = 0;
					prev_fit = curnt_fit;
					++m_count;
					m_log.push_back(m_pgbest->best_fit());
				}
			}
		}
		void reset_counter_generic()
		{
			m_count = 0;
		}
		double const_w() const
		{
			//return 0.5 + 0.5*((double)rand() / RAND_MAX); // rand_w()
			return 0.5;
		}
		double chaotic_w(int maxGen) const
		{
			double z = (double)rand() / RAND_MAX;
			z = 4 * z*(1 - z);
			return (0.9 - 0.4)*(maxGen - m_count) / maxGen + 0.4*z;
		}
		double constriction() const
		{
			double phi = std::max(c1 + c2, 4);
			return 2 / std::abs(2 - phi - std::sqrt(phi*phi - 4 * phi));
		}
		int m_count;
		Par_t m_fs;
		FitConv<Val_t>* m_pfit;
		std::set<Particle<Val_t>*> m_parts;
		Particle<Val_t>* m_pgbest;
		std::list<double> m_log;
	};
	
	// ======================DISCRETE=PSO========================================

	class ParamSet
	{
	public:
		ParamSet(const char* filename)
		{
			std::ifstream fin(filename, std::ios::in);
			if (!fin.is_open())
				throw std::invalid_argument("ParamSet::ParamSet()");

			for (m_ndim = 0; !fin.eof(); ++m_ndim) // no new line at end of the file!
				fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			if (m_ndim == 0)
			{
				fin.close();
				throw std::invalid_argument("ParamSet::ParamSet()");
			}
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
				buff.clear();
			}
			fin.close();
			for (int i = 0; i < m_ndim; ++i)
				if ((m_psets + i)->empty())
					throw std::invalid_argument("ParamSet::ParamSet()");
		}
		~ParamSet()
		{
			delete[] m_psets;
		}
		int ndim() const
		{
			return m_ndim;
		}
		int length(int dimNum) const
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
	class DFitConv :public FitConv<int>
	{
	public:
		DFitConv(IFitness* funcs, ParamSet* pfs, int* count) :m_pfuncs(funcs), m_pfs(pfs), m_pcount(count)
		{}
		double operator() (int* pos, int size) const
		{
			std::vector<double> values;
			for (int i = 0; i < size; ++i)
				values.push_back(m_pfs->value(i, *(pos + i)));

			double alpha = rand() % 1900 + 100;
			double H = 0;
			std::set<double> con = m_pfuncs->constraint(values);

			for (std::set<double>::const_iterator it = con.begin(); it != con.end(); ++it)
			{
				double gamma = (*it < 1) ? 1.0 : 2.0;
				double theta;
				if (*it < 0.00001)		theta = alpha;
				else if (*it < 0.001)	theta = 10 * alpha;
				else if (*it < 1)		theta = 100 * alpha;
				else					theta = 1000 * alpha;
				H += (theta * std::pow(std::max(0.0, *it), gamma));
			}
			return m_pfuncs->objective(values) + std::pow((*m_pcount + 1), 0.1) * H;
		}
	private:
		IFitness* m_pfuncs;
		ParamSet* m_pfs;
		int* m_pcount;
	};
	class DParticle :public Particle<int>
	{
	public:
		DParticle(int numdim, int* maxpos, Particle<int>*& gbest) :Particle<int>(numdim, gbest), m_pmaxpos(maxpos)
		{
			m_ppos = new int[m_ndim];
			m_pvel = new int[m_ndim];
			m_lbest.second = new int[m_ndim];
			
			for (int i = 0; i < m_ndim; ++i)
			{
				*(m_ppos + i) = rand() % (*(m_pmaxpos + i));
				*(m_lbest.second + i) = *(m_ppos + i);
			}
			for (int i = 0; i < m_ndim; ++i)
				*(m_pvel + i) = (int)(rand() % *(m_pmaxpos + i) * std::pow(-1, rand()));
		}
		void update_pos()
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
		}
	private:
		int* m_pmaxpos;		// upper bounds for position
	};

	typedef Pso_base<int, ParamSet> Discr_Pso_base;

	class DiscretePso :public Discr_Pso_base
	{
	public:
		DiscretePso(const char* filename, IFitness* funcs, int numPart) :Discr_Pso_base(filename, funcs, numPart)
		{
			m_pfit = new DFitConv(funcs, &m_fs, &m_count);

			m_pmaxpos = new int[m_fs.ndim()];
			for (int dim = 0; dim < m_fs.ndim(); ++dim)
				*(m_pmaxpos + dim) = m_fs.length(dim);

			m_pgbest = new DParticle(m_fs.ndim(), m_pmaxpos, m_pgbest);
			m_parts.insert(m_pgbest);
			for (int n = 1; n < numPart; ++n)
				m_parts.insert(new DParticle(m_fs.ndim(), m_pmaxpos, m_pgbest));
			for (std::set<Particle<int>*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
				(*it)->init_fit(*m_pfit);
		}
		~DiscretePso()
		{
			delete[] m_pmaxpos;
		}
		Pso* run(int numGen, Pso::Term criterion)
		{
			Discr_Pso_base::run_generic(numGen, criterion);
			return static_cast<Pso*>(this);
		}
		Pso* reset_counter()
		{
			Discr_Pso_base::reset_counter_generic();
			return static_cast<Pso*>(this);
		}
		std::pair<double, std::vector<double>> get_best() const
		{
			std::pair<double, std::vector<double>> temp;
			temp.first = m_pgbest->best_fit();
			for (int dim = 0; dim < m_fs.ndim(); ++dim)
				temp.second.push_back(m_fs.value(dim, m_pgbest->best_pos(dim)));
			return temp;
		}
	private:
		int* m_pmaxpos;
	};
	inline Pso* Pso::discrete(const char* filename, IFitness* funcs, int numPart)
	{
		return static_cast<Pso*>(new DiscretePso(filename, funcs, numPart));
	}
	
	// ======================CONTINUOUS=PSO======================================

	class ParamBounds
	{
	public:
		ParamBounds(const char* filename) :m_pbounds(nullptr)
		{
			std::ifstream fin(filename, std::ios::in);
			if (!fin.is_open())
				throw std::invalid_argument("ParamBounds::ParamBounds()");

			for (m_ndim = 0; !fin.eof(); ++m_ndim) // no new line at end of the file!
				fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			if (m_ndim == 0)
			{
				fin.close();
				throw std::invalid_argument("ParamSet::ParamSet()");
			}
			fin.seekg(0, fin.beg);

			m_pbounds = new std::pair<double, double>[m_ndim];

			for (int dim = 0; dim < m_ndim; ++dim)
				fin >> (m_pbounds + dim)->first >> (m_pbounds + dim)->second;
			fin.close();
		}
		~ParamBounds()
		{
			delete[] m_pbounds;
		}
		int ndim() const
		{
			return m_ndim;
		}
		double length(int dim) const
		{
			return (m_pbounds + dim)->second - (m_pbounds + dim)->first;
		}
		double min(int dim) const
		{
			return (m_pbounds + dim)->first;
		}
		double max(int dim) const
		{
			return (m_pbounds + dim)->second;
		}
	private:
		int m_ndim;
		std::pair<double, double>* m_pbounds;
	};
	class CFitConv :public FitConv<double>
	{
	public:
		CFitConv(IFitness* funcs, ParamBounds* pfs, int* count) :m_pfuncs(funcs), m_pfs(pfs), m_pcount(count)
		{}
		double operator() (double* pos, int size) const
		{
			std::vector<double> values;
			for (int i = 0; i < size; ++i)
				values.push_back(*(pos + i));

			double alpha = rand() % 1900 + 100;
			double H = 0;
			std::set<double> con = m_pfuncs->constraint(values);

			for (std::set<double>::const_iterator it = con.begin(); it != con.end(); ++it)
			{
				double gamma = (*it < 1) ? 1.0 : 2.0;
				double theta;
				if (*it < 0.00001)		theta = alpha;
				else if (*it < 0.001)	theta = 10 * alpha;
				else if (*it < 1)		theta = 100 * alpha;
				else					theta = 1000 * alpha;
				H += (theta * std::pow(std::max(0.0, *it), gamma));
			}
			return m_pfuncs->objective(values) + std::pow((*m_pcount + 1), 0.1) * H;
		}
	private:
		IFitness* m_pfuncs;
		ParamBounds* m_pfs;
		int* m_pcount;
	};
	class CParticle :public Particle<double>
	{
	public:
		CParticle(ParamBounds* pb, Particle<double>*& gbest) :Particle<double>(pb->ndim(), gbest), m_pb(pb)
		{
			m_ppos = new double[m_ndim];
			m_pvel = new double[m_ndim];
			m_lbest.second = new double[m_ndim];

			for (int i = 0; i < m_ndim; ++i)
			{
				*(m_ppos + i) = ((double)rand() / RAND_MAX) * m_pb->length(i) + m_pb->min(i);
				*(m_lbest.second + i) = *(m_ppos + i);
			}
			for (int i = 0; i < m_ndim; ++i)
				*(m_pvel + i) = ((double)rand() / RAND_MAX) * m_pb->length(i) * std::pow(-1, rand());
		}
		void update_pos()
		{
			for (int i = 0; i < m_ndim; ++i)
				*(m_ppos + i) += *(m_pvel + i);

			for (int i = 0; i < m_ndim; ++i)
				if (*(m_ppos + i) < m_pb->min(i))
				{
					*(m_ppos + i) = m_pb->min(i);
					*(m_pvel + i) *= -1;
				}
				else if (*(m_ppos + i) > m_pb->max(i))
				{
					*(m_ppos + i) = m_pb->max(i);
					*(m_pvel + i) *= -1;
				}
		}
	private:
		ParamBounds* m_pb;
	};

	typedef Pso_base<double, ParamBounds> Cont_Pso_base;
	
	class ContinuousPso :public Cont_Pso_base
	{
	public:
		ContinuousPso(const char* filename, IFitness* funcs, int numPart) :Cont_Pso_base(filename, funcs, numPart)
		{
			m_pfit = new CFitConv(funcs, &m_fs, &m_count);

			m_pgbest = new CParticle(&m_fs, m_pgbest);
			m_parts.insert(m_pgbest);
			for (int n = 1; n < numPart; ++n)
				m_parts.insert(new CParticle(&m_fs, m_pgbest));
			for (std::set<Particle<double>*>::iterator it = m_parts.begin(); it != m_parts.end(); ++it)
				(*it)->init_fit(*m_pfit);
		}
		Pso* run(int numGen, Pso::Term criterion)
		{
			Cont_Pso_base::run_generic(numGen, criterion);
			return static_cast<Pso*>(this);
		}
		Pso* reset_counter()
		{
			Cont_Pso_base::reset_counter_generic();
			return static_cast<Pso*>(this);
		}
		std::pair<double, std::vector<double>> get_best() const
		{
			std::pair<double, std::vector<double>> temp;
			temp.first = m_pgbest->best_fit();
			for (int dim = 0; dim < m_fs.ndim(); ++dim)
				temp.second.push_back(m_pgbest->best_pos(dim));
			return temp;
		}
	};
	inline Pso* Pso::continuous(const char * filename, IFitness * funcs, int numParts)
	{
		return static_cast<Pso*>(new ContinuousPso(filename, funcs, numParts));
	}
}