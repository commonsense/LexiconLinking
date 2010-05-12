
#ifndef	LPBOOST15OSI_H
#define	LPBOOST15OSI_H

/* Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 28th March 2007
 */

#include <vector>
#include <iostream>
#include <assert.h>

#include <Sequence.h>
#include <PrefixSpanLPBoost15.h>
#include <PrefixSpanLPBoost.h>
#include <PrefixSpan.h>

namespace LPBoost
{

/* 1.5 class LPBoost
 *
 * Solve over variables \alpha, \rho_1, \rho_2, \xi_1, \xi_2,
 *
 *    min    \rho_2 - \rho_1 + \frac{1}{\nu N} \sum_{n=1}^{N} \xi_{1,n}
 *                           + \frac{1}{\nu M} \sum_{m=1}^{M} \xi_{2,m}
 *    sb.t.  \sum_{t \in T} \alpha_t h(x_{1,n} ; t) - \rho_1 + \xi_{1,n} >= 0,
 *                                                  n = 1, ..., N
 *           \sum_{t \in T} \alpha_t h(x_{2,m} ; t) - \rho_2 - \xi_{2,m} <= 0,
 *                                                  m = 1, ..., M
 *           \sum_{t \in T} \alpha_t = 1,   (or <= 1)
 *           (optional) \alpha_t <= C,   t \in T.
 *
 *    vars, bounds
 *           \alpha_t >= 0,    t \in T,
 *           \xi_{1,n} >= 0,   n = 1, ..., N,
 *           \xi_{2,m} >= 0,   m = 1, ..., M,
 *           \rho_1, \rho_2 \in \R.
 *
 * The dual LP is solved over the variables \lambda, \mu and \gamma,
 *
 *    max    \gamma - C(\sum_{i=1}^K \iota_i)
 *    sb.t.  \sum_{n=1}^{N} \lambda_n h(x_{1,n} ; t)
 *              - \sum_{m=1}^{M} \mu_m h(x_{2,m} ; t) + \gamma <= 0,   t \in T,
 *           \sum_{n=1}^{N} \lambda_n = 1,
 *           \sum_{m=1}^{M} \mu_m = 1,
 *           0 <= \lambda_n <= \frac{1}{\nu N},   n = 1, ..., N,
 *           0 <= \mu_m <= \frac{1}{\nu M},   m = 1, ..., M.
 *           (optional) \gamma <= 0
 *
 * which is the actual problem we solve.
 */
class
LPBoost15OSI
{
private:
	SeqMine::PrefixSpanOptions options;

	double convergenceEpsilon;

	double nu;
	unsigned int pricingCount;
	unsigned int pos_minimum_support;
	unsigned int exploration_max;

	/* alpha_i <= C bound, extra handling.
	 * If alphaBoundC = 1.0, the original 1.5-class case is recovered.
	 */
	double alphaBoundC;

	/* Validation set.
	 */
	bool doValidation;
	std::vector<const SeqMine::Sequence*> Sval;
	std::vector<double> Yval;

	/* Consistency check on Wcur and alpha.
	 */
	bool VerifyCurrentStatus (unsigned int samplecount,
		std::vector<double>& Wcur, double D, std::vector<double>& alpha) const;

	void PrintCurrentTrainingStatus (
		const std::vector<const SeqMine::Sequence*>& Strain,
		std::vector<double>& Ylabels) const;

public:
	std::vector<double> alpha;
	double rho1, rho2;
	std::vector<std::pair<SeqMine::SequenceT, double> > stumps;

	LPBoost15OSI (SeqMine::PrefixSpanOptions options,
		double nu, unsigned int pricingCount,
		double convergenceEpsilon = 1e-4,
		unsigned int pos_minimum_support = 0,
		unsigned int exploration_max = 0,
		double alphaBoundC = 1.0)
		: options (options),
			convergenceEpsilon (convergenceEpsilon), nu (nu),
			pricingCount (pricingCount),
			pos_minimum_support (pos_minimum_support),
			exploration_max (exploration_max),
			alphaBoundC (alphaBoundC)
	{
		assert (alphaBoundC > 1e-3 && alphaBoundC <= 1.0);

		doValidation = false;
	}

	/* Train a 1.5 class LPBoost classifier.
	 *
	 * Input
	 *    S: N+M training sequences, but not necessarily in order.
	 *    Ylabels: N+M vector of labels in { -1, 1 }.
	 *    searchStrategy: see PrefixSpanLPBoost.h
	 *
	 * After training, the public variables 'alpha' and 'stumps' are the final
	 * classifier state.
	 */
	void BuildClassifier (std::vector<const SeqMine::Sequence*>& S,
		std::vector<double>& Ylabels,
		SeqMine::PrefixSpanLPBoost::SearchStrategy searchStrategy =
			SeqMine::PrefixSpanLPBoost::DepthFirstSearch);

	/* Set a validation set to evaluate the error on after each boosting
	 * iteration.  This method has to be called before BuildClassifier.
	 */
	void SetValidationSet (std::vector<const SeqMine::Sequence*>& Sval,
		std::vector<double>& Yval);

	/* Compute the validation score for an arbitrary set of sequences.
	 */
	double ValidationLossGeneral (const std::vector<const SeqMine::Sequence*>& S,
		const std::vector<double>& Y, double& posErr, double& negErr) const;

	/* Auxiliary debug function, printing the indicator array for a set of
	 * sequences.  Only for debug purposes.
	 */
	void DumpGY (std::ostream& os, std::vector<const SeqMine::Sequence*>& S,
		std::vector<double>& Y, int type = 2) const;

	/* sign (\sum_{t \in T} \alpha_t h(x ; t) - (rho1 - rho2)/2)
	 */
	double ClassDecision (double classifierOutput) const;
};

};

#endif


