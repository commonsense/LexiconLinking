/* LPBoost 1.5 class Linear Programming Boosting implementation
 * Copyright (C) 2007  Sebastian Nowozin
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <limits>
#include <algorithm>

/* Open Solver Interface (OSI), http://projects.coin-or.org/Osi
 * include files
 */
#include <OsiClpSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>
#include <CoinPackedVector.hpp>

#include <LPBoost15OSI.h>

namespace LPBoost
{

void
LPBoost15OSI::BuildClassifier (std::vector<const SeqMine::Sequence*>& S,
	std::vector<double>& Ylabels,
	SeqMine::PrefixSpanLPBoost::SearchStrategy searchStrategy)
{
	/* Initialize Wcur uniformly, but class-balanced, such that the weights
	 * sum to one.
	 */
	std::vector<double> Wcur (Ylabels.size ());
	unsigned int posCount = std::count (Ylabels.begin (), Ylabels.end (), 1.0);
	unsigned int negCount = std::count (Ylabels.begin (), Ylabels.end (), -1.0);
	for (unsigned int n = 0 ; n < Wcur.size () ; ++n)
		Wcur[n] = 1.0 / (2.0 * ((Ylabels[n] >= 0.0) ? ((double) posCount)
			: ((double) negCount)));

	/* nu regularization parameter
	 */
	double D_pos = 1.0 / (nu * posCount);
	double D_neg = 1.0 / (nu * negCount);

	/* The two soft margins
	 */
	rho1 = rho2 = 0.0;

	/* Variables: [lambdamu, gamma, iota]
	 *
	 * Initially iota does not contain any elements and it grows with the
	 * number of constraints.
	 */
	unsigned int iota_count = 0;
	OsiSolverInterface* si = new OsiClpSolverInterface;

	CoinPackedMatrix* matrix = new CoinPackedMatrix (false, 0, 0);
	matrix->setDimensions (0, Ylabels.size()+1+iota_count);

	/* Add the constraints
	 *    \sum_{n=1}^{N} \lambda_n = 1,
	 *    \sum_{m=1}^{M} \mu_m = 1.
	 */
	CoinPackedVector lambdasum;
	CoinPackedVector musum;
	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
		if (Ylabels[n] >= 0.0)
			lambdasum.insert (n, 1.0);
		else
			musum.insert (n, 1.0);
	}

	lambdasum.insert (Ylabels.size (), 0.0);
	musum.insert (Ylabels.size (), 0.0);
	matrix->appendRow (lambdasum);
	matrix->appendRow (musum);

	/* As rho_1, rho_2 are unconstrained we have equality constraints here.
	 * (Unconstrained primal variables become equality constraints in the
	 * dual.)
	 */
	double lambdamusumLB[] = { 1.0, 1.0, };
	double lambdamusumUB[] = { 1.0, 1.0, };

	/* Primal variables (lambda, gamma) lower/upper bounds
	 */
	double varLB[Ylabels.size()+1];
	double varUB[Ylabels.size()+1];
	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
		varLB[n] = 0;

		if (Ylabels[n] >= 0.0)
			varUB[n] = D_pos;
		else
			varUB[n] = D_neg;
	}

	/* gamma unconstrained.
	 */
	varLB[Ylabels.size()] = -si->getInfinity ();
	varUB[Ylabels.size()] = 0.0;
	//varUB[Ylabels.size()] = si->getInfinity ();

	/* Number of stumps in the LP so far
	 */
	unsigned int stump_count = 0;
	stumps.clear ();

	/* Set objective function: f(lambdamu,gamma) = gamma
	 *
	 * the variable ordering we use is [lambdamu, gamma]
	 */
	double objective[Ylabels.size() + 1];
	objective[Ylabels.size()] = 1.0;	// max: gamma
	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n)
		objective[n] = 0.0;	// 0.0*lambdamu_n

	si->loadProblem (*matrix, varLB, varUB, objective,
		lambdamusumLB, lambdamusumUB);
	si->setObjSense (-1);	// maximize (-1), minimize (1)

	/* theta is the minimum required gain (see lpboost.m)
	 */
	double theta = std::numeric_limits<double>::min ();

	for (unsigned int iter = 0 ; true ; ++iter) {
		std::cout << "LPBoost 1.5-class, iter " << iter
			<< ", theta " << theta << std::endl;

		/* Find the pattern stump with the highest gain for the current
		 * weights Wcur.
		 */
		SeqMine::PrefixSpanLPBoost15 boost (options, Ylabels, Wcur, theta,
			pricingCount, exploration_max, pos_minimum_support, searchStrategy);
		boost.Mine (S);

		/* In case no patterns were found, we are done.
		 */
		if (boost.bestPatterns.size () == 0)
			break;

		/* Convergence check: the new best sequence has to be better than our
		 * current performance.
		 */
		std::cout << "best gain: " << boost.bestPatterns.begin ()->first << std::endl;
		if (boost.bestPatterns.begin ()->first <= theta + convergenceEpsilon) {
			std::cout << "Converged, there is no better pattern." << std::endl;

			break;
		}

		/* Add the new stumps to the stump set H
		 */
		stump_count += boost.bestPatterns.size ();

		/* Add the mined patterns (new constraints)
		 */
		for (std::multimap<double, std::pair<SeqMine::SequenceT, double> >::iterator
			iv = boost.bestPatterns.begin () ; iv != boost.bestPatterns.end () ; ++iv)
		{
			/* Add to global list of stumps
			 */
			stumps.push_back (iv->second);

			CoinPackedVector row;
			for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
				/* Predicate: does sample n contain the current pattern?
				 */
				bool sampleActive = (iv->second.first.Occurence.find (n)
					!= iv->second.first.Occurence.end ());
				double stumpResponse = (sampleActive ? 1.0 : 0.0);

				row.insert (n, Ylabels[n] * stumpResponse);
			}
			row.insert (Ylabels.size (), 1.0);	// gamma
			row.insert (Ylabels.size ()+1+(stumps.size () - 1), -1.0);	// iota_i

			/* Add another iota variable with all zero coefficients in the
			 * constraint matrix
			 */
			CoinPackedVector iotaVec;
			si->addCol (iotaVec, 0.0, si->getInfinity (), -alphaBoundC);

			/* Each iteration we add one row of this kind:
			 *
			 * \sum_{n=1}^{\ell} y_n h(x_n ; t) lambdamu_n - \iota_i + \gamma <= 0
			 */
			si->addRow (row, -si->getInfinity(), 0.0);
		}

		/* Solve the LP
		 */
		si->initialSolve ();
		if (si->isProvenOptimal () == false) {
			std::cerr << "Solve failed." << std::endl;

			std::cerr << "Problem: " << si->getNumCols () << " variables, "
				<< si->getNumRows () << " rows." << std::endl;
			std::cerr << "         problem sense: " << si->getObjSense () << std::endl;

			std::cerr << "STATUS:  numerical difficulties: "
				<< (si->isAbandoned () ? "YES" : "no") << std::endl;
			std::cerr << "STATUS:       primal infeasible: "
				<< (si->isProvenPrimalInfeasible () ? "YES" : "no") << std::endl;
			std::cerr << "STATUS:         dual infeasible: "
				<< (si->isProvenDualInfeasible () ? "YES" : "no") << std::endl;
			std::cerr << "STATUS: iteration limit reached: "
				<< (si->isIterationLimitReached () ? "YES" : "no") << std::endl;

			assert (0);
		}

		const double* primal = si->getColSolution ();
		const double* dual = si->getRowPrice ();

		/* The solution contains the primal variables (lambdamu, gamma), as
		 * well as the dual variables (the original primal variables, but we
		 * are only interested in alpha and rho, where we get alpha directly
		 * from the Lagrangian multipliers of the constraints and rho from the
		 * objective function)
		 */

		/* Alpha, the stump multipliers
		 */
		alpha.resize (stump_count);
		std::cout << "alpha:";
		for (unsigned int n = 0 ; n < stump_count ; ++n) {
			alpha[n] = dual[2+n];
			std::cout << " " << alpha[n];
		}
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "iota:";
		for (unsigned int n = 0 ; n < stump_count ; ++n) {
			std::cout << " " << primal[Ylabels.size()+1+n];
		}
		std::cout << std::endl;
		std::cout << std::endl;

		/* Track changes in weights from iteration to iteration
		 */
		std::vector<double> Wlast (Wcur);

		double wdiff = 0.0;
		double wsum = 0.0;
		std::cout << "weights:";
		for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
			Wcur[n] = primal[n];
			wdiff += fabs (Wcur[n] - Wlast[n]);
			wsum += Wcur[n];

			std::cout << " " << Wcur[n] << "(" << Ylabels[n] << ")";
		}
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "weight sum: " << wsum << std::endl;
		std::cout << "weight change: " << wdiff << std::endl;

		/* Two soft margin values, rho1 and rho2
		 */
		rho1 = -dual[0];	// flipped constraint, hence negation
		rho2 = dual[1];

		/* We use max: \gamma,
		 * but our constraints are of the form
		 *
		 *   \sum_{n=1}^{\ell} y_n h(x_n ; t, \omega) \lambda_n + \gamma <= 0,
		 *
		 * such that -\gamma is the right hand side in the nu-LPBoost
		 * formulations of [Demirez2002] and [Raetsch2001] and hence
		 * \theta=-\gamma, the gain we have to beat using our mining.
		 */
		theta = -primal[Ylabels.size ()];	// gamma is <= 0
		std::cout << "theta = " << theta
			<< ", obj = " << si->getObjValue () << std::endl;
		std::cout << std::endl;

		PrintCurrentTrainingStatus (S, Ylabels);
		std::cout << "(rho1, rho2) = (" << rho1 << ", " << rho2 << ")" << std::endl;

		//assert (0);

		if (wdiff < 1e-6) {
			std::cout << "Stopped due to identical weights." << std::endl;

			break;
		}
	}

	std::cout << "Final classification" << std::endl;
	PrintCurrentTrainingStatus (S, Ylabels);

	delete matrix;
	delete si;
}


double
LPBoost15OSI::ValidationLossGeneral (const std::vector<const SeqMine::Sequence*>& S,
	const std::vector<double>& Y, double& posErr, double& negErr) const
{
	double loss = 0.0;
	std::vector<double> sOut (Y.size ());

	/* Sum the classifier output.
	 */
	for (unsigned int sn = 0 ; sn < alpha.size () ; ++sn) {
		for (unsigned int n = 0 ; n < Y.size () ; ++n) {
			bool does_contain = (options.maxgap < 0) ?
				S[n]->contains (stumps[sn].first) :
				S[n]->containsGap (stumps[sn].first, options.maxgap);

			sOut[n] += alpha[sn] * (does_contain ? 1.0 : 0.0);
		}
	}

	/* Compute the empirical loss
	 */
	posErr = 0.0;
	negErr = 0.0;
	for (unsigned int n = 0 ; n < Y.size () ; ++n) {
		double cl = ClassDecision (sOut[n]);
		if (cl * Y[n] <= 0.0) {
			loss += 1.0;

			if (Y[n] > 0.0)
				posErr += 1.0;
			else
				negErr += 1.0;
		}
	}
	loss /= (double) Y.size ();

	posErr /= (double) std::count (Y.begin (), Y.end (), 1.0);
	negErr /= (double) std::count (Y.begin (), Y.end (), -1.0);

	return (loss);
}

void
LPBoost15OSI::SetValidationSet (std::vector<const SeqMine::Sequence*>& Sval,
	std::vector<double>& Yval)
{
	doValidation = true;

	this->Sval = Sval;
	this->Yval = Yval;
}

double
LPBoost15OSI::ClassDecision (double classifierOutput) const
{
	double p = classifierOutput - (rho1 + rho2)/2.0;

	if (p >= 0.0)
		return (1);
	else
		return (-1);
}

void
LPBoost15OSI::PrintCurrentTrainingStatus (
	const std::vector<const SeqMine::Sequence*>& Strain,
	std::vector<double>& Ylabels) const
{
	double trainingLoss = 0.0;

	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
#if 0
		std::cout << "sample " << n << ": ";
#endif

		double sum = 0.0;
		for (unsigned int d = 0 ; d < alpha.size () ; ++d) {
#if 0
			std::cout << " " << alpha[d] << "*";
#endif

			bool active = (stumps[d].first.Occurence.find (n) !=
				stumps[d].first.Occurence.end ());
			double response = (active ? 1.0 : 0.0);
#if 0
			std::cout << response;
#endif

			sum += alpha[d] * response;
		}
#if 0
		std::cout << " = " << sum << " (true " << Ylabels[n] << ")";
#endif

		if (ClassDecision (sum) * Ylabels[n] <= 1e-4) {
#if 0
			std::cout << " INCORRECT";
#endif
			trainingLoss += 1.0;
		}
#if 0
		std::cout << std::endl;
#endif
	}
	trainingLoss /= Ylabels.size ();
	std::cout << "TOTAL TRAINING LOSS: " << trainingLoss
		<< "   for (" << std::count (Ylabels.begin (), Ylabels.end (), 1.0)
		<< " pos, " << std::count (Ylabels.begin (), Ylabels.end (), -1.0)
		<< " neg) samples" << std::endl;
	double posErr, negErr;
	std::cout << "   2. TRAINING LOSS: "
		<< ValidationLossGeneral (Strain, Ylabels, posErr, negErr) << std::endl;
	std::cout << "         (pos, neg): " << posErr << ", " << negErr << std::endl;

	if (doValidation) {
		std::cout << "TOTAL VALIDATN LOSS: "
			<< ValidationLossGeneral (Sval, Yval, posErr, negErr)
			<< "  for (" << std::count (Yval.begin (), Yval.end (), 1.0)
			<< " pos, " << std::count (Yval.begin (), Yval.end (), -1.0)
			<< " neg) samples" << std::endl;
		std::cout << "         (pos, neg): " << posErr << ", " << negErr << std::endl;
	}
}


void
LPBoost15OSI::DumpGY (std::ostream& os, std::vector<const SeqMine::Sequence*>& S,
	std::vector<double>& Y, int type) const
{
	for (unsigned int n = 0 ; n < Y.size () ; ++n) {
		os << Y[n];

		for (unsigned int d = 0 ; d < alpha.size () ; ++d) {
			bool active;

			// FIXME: maxgap predicate, maybe coalesce with LPBoostOSI.cpp
			if (type == 1) {
				active = (stumps[d].first.Occurence.find (n) !=
					stumps[d].first.Occurence.end ());
			} else {
				active = S[n]->contains (stumps[d].first);
			}
			double response = (active ? 1.0 : 0.0);

			os << " " << response;
		}
		os << std::endl;
	}
}

}


