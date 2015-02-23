#ifndef PARTSWARMOPT_H
#define PARTSWARMOPT_H

// Class for carrying out optimisation using the particle
// swarm optimization algorithm (PSO)

#include<fstream>
#include<string>
#include<iostream>
#include<vector>
#include<cmath>
#include<ctime>
using namespace std;


// void debug_print(const vector<double> vec); // DEBUGGING


class Udm
{ // User-defined methods and constants
public:

// Prototype of the objective function
// Precondition: parameter values in the same order as in the input file
// Postcondition: objective function value
static double Objective(const vector<double>& parameterValues);

// Prototype of the constraint function
// Precondition: parameter values in the same order as in the input file
// Postcondition: constraints values, (0,1] = infeasible, [-1,0] = feasible
// If no constraints applied the function returns an empty vector
static vector<double> Constraints(const vector<double>& parameterValues);

static const int c1 = 1; // local learning factor
static const int c2 = 1; // global learning factor
static const int NOG = 500; // assumed number of generations for calculating inertia weight
};


class Feasible_set
{ // Not to be used by user
public:
	Feasible_set()
	{
		m_numberOfSets = 0;
	}

	// Adds new parameter range at the end of m_data
	void add_set(const vector<double>& newset)
	{
		unsigned int size = newset.size();
		for (unsigned int i = 0; i < size; i++)
		{
			m_data.push_back(newset[i]);
		}
		m_numbersInEachSet[m_numberOfSets] = size;
		m_numberOfSets++;
	}

	// Retrieve parameter value (value number from 1)
	double get_value(int setNumber, int valueNumber) const
	{
		int position = -1; // because m_data starting index is 0
		for (int i = 0; i < setNumber - 1; i++)
		{
			position += m_numbersInEachSet[i];
		}
		position += valueNumber;
		return m_data[position];
	}

	// Retrieve set length (one of the sets)
	int get_length(int setnumber) const
	{
		return m_numbersInEachSet[setnumber - 1];
	}

	// Get the number of dimensions
	int numberOfRanges() const
	{
		return m_numberOfSets;
	}

private:
	vector<double> m_data;
	int m_numbersInEachSet[500];
	int m_numberOfSets;
};


class Pso
{ // To be used by user
public:
	Pso()
	{ }

	// Full constructor, takes input file name/path and number of particles
	Pso(string filename, unsigned int numberOfPart) :m_inputSource(filename), m_numberOfParticles(numberOfPart)
	{
		m_seed = static_cast<int>(time(0));
		readIn(filename);
	}

	// Initialize particles
	void initialize()
	{
		initializeAll();
	}

	// Runs new generation for all particles
	Pso& operator ++ ()
	{
		m_currentGeneration++;

		for (unsigned int i = 0; i < m_numberOfParticles; i++)
		{
			iterate(m_theParticles[i]);
		}
		determineGBest();
		return *this;
	}

	// Print output using << operator
	friend ostream& operator <<(ostream& out, const Pso& opt);

	// Outputs current global best
	double bestFitness() const
	{
		return m_currentGBest;
	}

	// Outputs current best position
	vector<int> bestPosision() const
	{
		return m_gbestPosition;
	}

	// Accessor+mutator to position of a particle
	vector<int>& currentPosition(unsigned int particleNumber)
	{
		return m_theParticles[particleNumber-1].position;
	}

	// Outputs global best parameter values
	vector<double> bestValues() const
	{
		vector<double> values(m_numberOfDimensions);
		for (unsigned int i = 0; i < m_numberOfDimensions; i++)
		{ // Retrieving parameter values at the current position
			values[i] = m_myranges.get_value(i + 1, m_gbestPosition[i]);
		}
		return values;
	}

private:
	struct particle
	{
		vector<int>position;
		vector<int>velocity;
		vector<int>pbestPosition;
		double pbest;
	};


	string m_inputSource;
	Feasible_set m_myranges;
	int m_seed;

	unsigned int m_numberOfDimensions;
	unsigned int m_numberOfParticles;

	int m_currentGeneration;
	double m_currentGBest;
	vector<int> m_gbestPosition;
	vector<int> m_maxVelocityAllowed;
	vector<particle> m_theParticles;


	// Reads data ranges into Feasible_set 'm_myranges'
	void readIn(string filename)
	{
		ifstream sourceFile;
		sourceFile.open(filename, ios::in);
		if (!sourceFile.is_open()) {
			cerr << "Error: Input file cannot be open\n";
			exit(1);
		}
		m_myranges = Feasible_set();
		vector<double> datavec;
		int numberOfValues, value;
		while (sourceFile >> numberOfValues)
		{ // reading a line
			for (int i = 0; i < numberOfValues; i++)
			{
				sourceFile >> value;
				datavec.push_back(value);
			}
			m_myranges.add_set(datavec);
			datavec.clear();
		}
		m_numberOfDimensions = m_myranges.numberOfRanges();
	}

	// Calculates inertia weight at the moment
	double inertiaWeight()
	{
		return 0.7;
	}

	// Runs one generation for one particle
	void iterate(particle& mypart)
	{
		double r1;
		double r2;
		for (unsigned int i = 0; i < m_numberOfDimensions; i++)
		{ // Setting new velocity
			srand(m_seed);
			r1 = rand() / static_cast<double>(RAND_MAX);
			r2 = rand() / static_cast<double>(RAND_MAX);
			mypart.velocity[i] = static_cast<int>(inertiaWeight()*mypart.velocity[i] + Udm::c1*r1*(m_gbestPosition[i] - mypart.position[i])
				+ Udm::c2*r2*(mypart.pbestPosition[i] - mypart.position[i])+0.5); // round up
			if (mypart.velocity[i] > m_maxVelocityAllowed[i])
			{ // Avoids explosions
				mypart.velocity[i] = m_maxVelocityAllowed[i];
			}
			m_seed++;
		}
		for (unsigned int i = 0; i < m_numberOfDimensions; i++)
		{ // Setting new position
			mypart.position[i] += mypart.velocity[i];
			while ( mypart.position[i] > m_myranges.get_length(i + 1) )
			{ // Wall collisions
				mypart.position[i]--;
				mypart.velocity[i]--;
			}
			while (mypart.position[i] < 1)
			{ // Wall collisions
				mypart.position[i]++; 
				mypart.velocity[i]++;
			}
		}
		vector<double> values(m_numberOfDimensions);
		for (unsigned int i = 0; i < m_numberOfDimensions; i++)
		{ // Retrieving parameter values at the current position
			values[i] = m_myranges.get_value(i + 1, mypart.position[i]);
		}
		vector<double> constraints = Udm::Constraints(values);

		double objective = Udm::Objective(values);
		double fitness = objective;
		for (unsigned int i = 0; i < constraints.size(); i++)
		{ // Final penalty function value
			if (constraints[i] > 0)
			{
				fitness += exp(constraints[i]) * fabs(objective);
			}
		}
		if (fitness < mypart.pbest)
		{ // New local best
			mypart.pbest = fitness;
			mypart.pbestPosition = mypart.position;
		}
	}

	void determineGBest()
	{
		for (unsigned int i = 0; i < m_numberOfParticles; i++)
		{
			if (m_theParticles[i].pbest < m_currentGBest)
			{ // New global best
				m_currentGBest = m_theParticles[i].pbest;
				m_gbestPosition = m_theParticles[i].pbestPosition;
			}
		}
	}

	void initializeAll()
	{
		m_currentGeneration = 1;

		for (unsigned int i = 1; i <= m_numberOfDimensions; i++)
		{ // Setting velocity limit
			m_maxVelocityAllowed.push_back(static_cast<int>(0.2*m_myranges.get_length(i)+0.5));
		}

		vector<particle> temp_part(m_numberOfParticles);
		vector<int> temp_dim(m_numberOfDimensions);
		for (unsigned int i = 0; i < m_numberOfParticles; i++)
		{
			temp_part[i].position = temp_dim;
			temp_part[i].velocity = temp_dim;
		}
		m_theParticles = temp_part;

		double r1;
		int r2;
		for (unsigned int i = 0; i < m_numberOfParticles; i++)
		{ // Initial random velocity and position
			for (unsigned int j = 0; j < m_numberOfDimensions; j++)
			{
				srand(m_seed);
				r1 = rand() / static_cast<double>(RAND_MAX);
				r2 = rand();
				m_theParticles[i].velocity[j] = static_cast<int>(m_maxVelocityAllowed[j] * r1 *pow(-1, rand()));
				m_theParticles[i].position[j] = r2 % m_myranges.get_length(j + 1) + 1;
				m_seed++;
			}
		}
		// Calculating fitness
		vector<double> values(m_numberOfDimensions);
		for (unsigned int i = 0; i < m_numberOfParticles; i++)
		{
			for (unsigned int j = 0; j < m_numberOfDimensions; j++)
			{ // Retrieving parameter values at the current position
				values[j] = m_myranges.get_value(j + 1, m_theParticles[i].position[j]);
			}

			// debug_print(values); // DEBUGGING

			vector<double> constraints = Udm::Constraints(values);
			double objective = Udm::Objective(values);
			double fitness = objective;
			for (unsigned int i = 0; i < constraints.size(); i++)
			{ // Final penalty function value
				if (constraints[i] > 0)
				{
					fitness += exp(constraints[i]) * fabs(objective);
				}
			}
			constraints.clear();
			m_theParticles[i].pbest = fitness;
			m_theParticles[i].pbestPosition = m_theParticles[i].position;
		}
		// Setting global best
		m_currentGBest = m_theParticles[0].pbest;
		m_gbestPosition = m_theParticles[0].pbestPosition;
		determineGBest();
	}
};

ostream& operator <<(ostream& out, const Pso& opt)
{
	out << "Current generation: " << opt.m_currentGeneration << "\n"
		<< "Current global best: " << opt.m_currentGBest << "\n"
		<< "Global best parameter values: ";
	vector<double> values = opt.bestValues();
	for (unsigned int i = 0; i < opt.m_numberOfDimensions; i++)
	{
		out << values[i] << " ";
	}
	out << endl << endl;
	return out;
}


// DEBUGGING
//void debug_print(const vector<double> vec)
//{
//	cout << endl << endl;
//	if (vec.size() == 0)
//		cout << "Empty vector";
//	else
//	{
//		for (unsigned int i = 0; i < vec.size(); i++)
//		{
//			cout << vec[i] << " ";
//		}
//	}
//	cout << endl << endl;
//}

#endif