/* LPBoost Linear Programming Boosting implementation
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
//#include <OsiGlpkSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>
#include <CoinPackedVector.hpp>

#include <LPBoostOSI.h>

namespace LPBoost
{

void
LPBoostOSI::BuildClassifier (std::vector<const SeqMine::Sequence*>& S,
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
	double D = 1.0 / (nu * Ylabels.size ());

	/* Variables: [lambda, gamma]
	 */
	OsiSolverInterface* si = new OsiClpSolverInterface;
	//OsiSolverInterface* si = new OsiGlpkSolverInterface;

	CoinPackedMatrix* matrix = new CoinPackedMatrix (false, 0, 0);
	matrix->setDimensions (0, Ylabels.size()+1);

	/* Add the \sum_{n=1}^{\ell} \lambda_n >= 1 constraint
	 */
	CoinPackedVector lambdasum;
	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n)
		lambdasum.insert (n, 1.0);

	lambdasum.insert (Ylabels.size (), 0.0);
	matrix->appendRow (lambdasum);

	/* This is delicate: as we leave \rho unconstrained in the LPBoost primal,
	 * in order to dualize it we represent it by \rho := \rho_2 - \rho_1, such
	 * that we arrive in the canonical LP form.  Then, in the dual, the lambda
	 * sum is an equality constraint, such that we have a distribution over
	 * samples:
	 *
	 *   \sum_{n=1}^{\ell} \lambda_n = 1.
	 *
	 * See the footnote on page 8 in [Demiriz2002].
	 */
	double lambdasumLB = 1.0;
	double lambdasumUB = 1.0;
	//double lambdasumUB = si->getInfinity ();

	/* Primal variables (lambda, gamma) lower/upper bounds
	 */
	double varLB[Ylabels.size()+1];
	double varUB[Ylabels.size()+1];
	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
		varLB[n] = 0;
		varUB[n] = D;
	}

	/* gamma unconstrained.
	 */
	varLB[Ylabels.size()] = -si->getInfinity ();
	varUB[Ylabels.size()] = si->getInfinity ();

	/* Number of stumps in the LP so far
	 */
	unsigned int stump_count = 0;
	stumps.clear ();

	/* Set objective function: f(lambda,gamma) = gamma
	 *
	 * the variable ordering we use is [lambda, gamma]
	 */
	double objective[Ylabels.size() + 1];
	objective[Ylabels.size()] = 1.0;	// max: gamma
	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n)
		objective[n] = 0.0;	// 0.0*lambda_n

	si->loadProblem (*matrix, varLB, varUB, objective, &lambdasumLB, &lambdasumUB);
	si->setObjSense (-1);	// maximize

	/* theta is the minimum required gain (see lpboost.m)
	 */
	double theta = std::numeric_limits<double>::min ();

	for (unsigned int iter = 0 ; true ; ++iter) {
		std::cout << "LPBoost 2-class, iter " << iter
			<< ", theta " << theta << std::endl;

		/* Find the pattern stump with the highest gain for the current
		 * weights Wcur.
		 */
		SeqMine::PrefixSpanLPBoost boost (options, Ylabels, Wcur, theta,
			pricingCount, exploration_max, searchStrategy);
		boost.directional = directional;
		boost.Mine (S);

		/* In case no patterns were found, we are done.
		 */
		if (boost.bestPatterns.size () == 0)
			break;

#ifdef DEBUG
		/* Mining returned with some successful sequences that gain us
		 * something.  Lets take a look at those and where they appear.
		 */
		std::cout << "Obtained sequences" << std::endl;
		for (std::multimap<double, std::pair<SeqMine::SequenceT, double> >::iterator
			iv = boost.bestPatterns.begin () ; iv != boost.bestPatterns.end () ; ++iv)
		{
			std::cout << "gain " << iv->first << " (ybase "
				<< iv->second.second << "), seq: "
				<< iv->second.first << std::endl;

			std::cout << "   occurence: " << std::endl;
			for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
				bool sampleActive = (iv->second.first.Occurence.find (n)
					!= iv->second.first.Occurence.end ());

				if (sampleActive == false)
					continue;

				std::cout << n << "(" << Ylabels[n] << ") ";
			}
			std::cout << std::endl;
		}
#endif

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
				double stumpResponse =
					(sampleActive ? 1.0 : -1.0) * iv->second.second;

				row.insert (n, Ylabels[n] * stumpResponse);
			}
			row.insert (Ylabels.size(), 1.0);	// gamma

			/* Each iteration we add one row of this kind:
			 *
			 * \sum_{n=1}^{\ell} y_n h(x_n ; t,\omega) \lambda_n + \gamma <= 0
			 */
			si->addRow (row, -si->getInfinity(), 0.0);
		}

		/* DEBUG: output as MPS
		 */
		char lpfilename[256];
		sprintf (lpfilename, "lp-debug-%d", iter);
		#if 1
		si->writeMps (lpfilename, "mps", si->getObjSense ());
		#endif

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

		//int colN = si->getNumCols ();
		const double* primal = si->getColSolution ();
		//int rowN = si->getNumRows ();
		const double* dual = si->getRowPrice ();

		/* The solution contains the primal variables (lambda, gamma), as
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
			alpha[n] = dual[1+n];
			std::cout << " " << alpha[n];
		}
		std::cout << std::endl;
		std::cout << std::endl;

		double wsum = 0.0;
		std::cout << "weights:";
		for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
			Wcur[n] = primal[n];
			wsum += Wcur[n];

			std::cout << " " << Wcur[n] << "(" << Ylabels[n] << ")";
		}
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "weight sum: " << wsum << std::endl;

		/* We use max: \gamma,
		 * but our constraints are of the form
		 *
		 *   \sum_{n=1}^{\ell} y_n h(x_n ; t, \omega) \lambda_n + \gamma <= 0,
		 *
		 * such that -\gamma is the right hand side in the nu-LPBoost
		 * formulations of [Demirez2002] and [Raetsch2001] and hence
		 * \theta=-\gamma, the gain we have to beat using our mining.
		 */
		theta = -primal[Ylabels.size ()];
		std::cout << "theta = " << theta << std::endl;
		std::cout << std::endl;

		PrintCurrentTrainingStatus (S, Ylabels);
		bool valid = VerifyCurrentStatus (Ylabels.size (), Wcur, D, alpha);
		assert (valid == true);

		//assert (0);
	}

	delete matrix;
	delete si;
}


bool
LPBoostOSI::VerifyCurrentStatus (unsigned int samplecount,
	std::vector<double>& Wcur, double D, std::vector<double>& alpha) const
{
	bool res = true;

	double Wsum = 0.0;
	for (std::vector<double>::iterator iv = Wcur.begin () ;
		iv != Wcur.end () ; ++iv)
	{
		Wsum += *iv;
		if (*iv < -1e-4 || *iv > (D+1e-4)) {
			std::cerr << "VerifyCurrentStatus: weight " << *iv
				<< " not within [0 ; " << D << "]" << std::endl;
			res = false;
		}
	}

	if (fabs (Wsum - 1.0) > 1e-4) {
		std::cerr << "VerifyCurrentStatus: weight sum " << Wsum
			<< " != 1.0" << std::endl;
		res = false;
	}

	double alphasum = 0.0;
	for (std::vector<double>::iterator iv = alpha.begin () ;
		iv != alpha.end () ; ++iv)
	{
		alphasum += *iv;
		if (*iv < -1e-4 || *iv > (1.0+1e-4)) {
			std::cerr << "VerifyCurrentStatus: alpha " << *iv
				<< " not within [0 ; 1]" << std::endl;
			res = false;
		}
	}

	if (fabs (alphasum - 1.0) > 1e-4) {
		std::cerr << "VerifyCurrentStatus: alphasum " << alphasum
			<< " != 1.0" << std::endl;
		res = false;
	}

	return (res);
}


#if 0
double
LPBoostOSI::ValidationLoss (void) const
{
	assert (doValidation);

	double loss = 0.0;
	std::vector<double> sOut (Yval.size ());

	/* Sum the classifier output.
	 */
	for (unsigned int sn = 0 ; sn < alpha.size () ; ++sn) {
		SeqMine::SequenceTSearch st (stumps[sn].first);

		for (unsigned int n = 0 ; n < Yval.size () ; ++n)
			sOut[n] += alpha[sn] * stumps[sn].second *
				(st.containsSequence (*Sval[n]) ? 1.0 : -1.0);
	}

	/* Compute the empirical loss
	 */
	for (unsigned int n = 0 ; n < Yval.size () ; ++n) {
		double margin = sOut[n] * Yval[n];
		if (margin <= 1e-4)
			loss += 1.0;
	}
	loss /= (double) Yval.size ();

	return (loss);
}
#endif


double
LPBoostOSI::ValidationLossGeneral (const std::vector<const SeqMine::Sequence*>& S,
	const std::vector<double>& Y) const
{
	double loss = 0.0;
	std::vector<double> sOut (Y.size ());

	/* Sum the classifier output.
	 */
	for (unsigned int sn = 0 ; sn < alpha.size () ; ++sn) {
		#if 0
		SeqMine::SequenceTSearch st (stumps[sn].first);

		for (unsigned int n = 0 ; n < Y.size () ; ++n)
			sOut[n] += alpha[sn] * stumps[sn].second *
				(st.containsSequence (*S[n]) ? 1.0 : -1.0);
		#endif
		for (unsigned int n = 0 ; n < Y.size () ; ++n) {
			bool does_contain = (options.maxgap < 0) ?
				S[n]->contains (stumps[sn].first) :
				S[n]->containsGap (stumps[sn].first, options.maxgap);

			sOut[n] += alpha[sn] * stumps[sn].second *
				(does_contain ? 1.0 : (directional ? 0.0 : -1.0));
		}
	}

	/* Compute the empirical loss
	 */
	for (unsigned int n = 0 ; n < Y.size () ; ++n) {
		double margin = sOut[n] * Y[n];
		if (margin <= 1e-4)
			loss += 1.0;
	}
	loss /= (double) Y.size ();

	return (loss);
}

void
LPBoostOSI::SetValidationSet (std::vector<const SeqMine::Sequence*>& Sval,
	std::vector<double>& Yval)
{
	doValidation = true;

	this->Sval = Sval;
	this->Yval = Yval;
}


void
LPBoostOSI::PrintCurrentTrainingStatus (
	const std::vector<const SeqMine::Sequence*>& Strain,
	std::vector<double>& Ylabels) const
{
	double trainingLoss = 0.0;
	double realizedMargin = std::numeric_limits<double>::max ();

	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
		std::cout << "sample " << n << ": ";

		double sum = 0.0;
		for (unsigned int d = 0 ; d < alpha.size () ; ++d) {
			std::cout << " " << alpha[d] << "*";

			bool active = (stumps[d].first.Occurence.find (n) !=
				stumps[d].first.Occurence.end ());

			double response;
			if (directional) {
				response = stumps[d].second * (active ? 1.0 : 0.0);
			} else {
				response = stumps[d].second * (active ? 1.0 : -1.0);
			}
			std::cout << response;

			sum += alpha[d] * response;
		}
		std::cout << " = " << sum << " (true " << Ylabels[n] << ")";

		double singleMargin = Ylabels[n] * sum;

		/* Realized margin:
		 * \rho(\alpha) = \min_{n=1,...,\ell} y_n \sum_{\alpha_j} \alpha_j h_j(x_n)
		 */
		if (singleMargin < realizedMargin)
			realizedMargin = singleMargin;

		if (singleMargin <= 1e-4) {
			std::cout << " INCORRECT";
			trainingLoss += 1.0;
		}
		std::cout << std::endl;
	}
	trainingLoss /= Ylabels.size ();
	std::cout << "TOTAL TRAINING LOSS: " << trainingLoss << std::endl;
	std::cout << "   2. TRAINING LOSS: "
		<< ValidationLossGeneral (Strain, Ylabels) << std::endl;
	std::cout << "    realized margin: " << realizedMargin << std::endl;

	if (doValidation) {
		std::cout << "TOTAL VALIDATN LOSS: "
			<< ValidationLossGeneral (Sval, Yval) << std::endl;
	}
}


void
LPBoostOSI::DumpGY (std::vector<const SeqMine::Sequence*>& S,
	std::vector<double>& Y, int type) const
{
	for (unsigned int n = 0 ; n < Y.size () ; ++n) {
		std::cout << Y[n];

		for (unsigned int d = 0 ; d < alpha.size () ; ++d) {
			bool active;

			// TODO: maxgap predicate
			if (type == 1) {
				active = (stumps[d].first.Occurence.find (n) !=
					stumps[d].first.Occurence.end ());
			} else {
				active = S[n]->contains (stumps[d].first);
			}
			double response = stumps[d].second * (active ? 1.0 : -1.0);

			std::cout << " " << response;
		}
		std::cout << std::endl;
	}
}

}


