
#ifndef	LPBOOSTOSI_H
#define	LPBOOSTOSI_H

#include <vector>
#include <iostream>

#include <Sequence.h>
#include <PrefixSpanLPBoost.h>
#include <PrefixSpan.h>

namespace LPBoost
{

class
LPBoostOSI
{
private:
	SeqMine::PrefixSpanOptions options;

	double convergenceEpsilon;

	double nu;
	unsigned int pricingCount;
	unsigned int exploration_max;

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
	std::vector<std::pair<SeqMine::SequenceT, double> > stumps;

	/* Use one-sided decision stumps?
	 */
	bool directional;

	LPBoostOSI (SeqMine::PrefixSpanOptions options,
		double nu, unsigned int pricingCount,
		double convergenceEpsilon = 1e-4, unsigned int exploration_max = 0)
		: options (options),
			convergenceEpsilon (convergenceEpsilon),
			nu (nu), pricingCount (pricingCount),
			exploration_max (exploration_max),
			doValidation (false), directional (false)
	{
	}

	void BuildClassifier (std::vector<const SeqMine::Sequence*>& S,
		std::vector<double>& Ylabels,
		SeqMine::PrefixSpanLPBoost::SearchStrategy searchStrategy =
			SeqMine::PrefixSpanLPBoost::DepthFirstSearch);

	/* Set a validation set to evaluate the error on after each boosting
	 * iteration.
	 */
	void SetValidationSet (std::vector<const SeqMine::Sequence*>& Sval,
		std::vector<double>& Yval);

	/* Compute the validation score.
	 */
	//double ValidationLoss (void) const;

	double ValidationLossGeneral (
		const std::vector<const SeqMine::Sequence*>& S,
		const std::vector<double>& Y) const;

	void DumpGY (std::vector<const SeqMine::Sequence*>& S,
		std::vector<double>& Y, int type) const;
};

};

#endif


