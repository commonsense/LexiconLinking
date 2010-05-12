/* mexprefixspan.cpp - MEX interface for (weighted) PrefixSpan
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 26rd February 2007
 *
 * Usage
 * =====
 * 1. For frequent subsequence mining
 *
 *    [subseqs, occurence] = mexprefixspan (S, minsup, size_req);
 *
 * 2. For weighted subsequence mining (sequence boosting)
 *
 *    [subseqs, ybase, GY] = mexprefixspan (S, minsup, size_req, ...
 *        Y, weights, tau, pricing_count, htype);
 *
 ****
 * Input
 *
 * S: (1,n) cellarray of n sequences, each sequence being a cellarray of
 *    uint32 (1,p_i) arrays.
 * minsup: The minimum required subsequence support.
 * size_req: (1,2) real vector containing lower/upper bounds on subsequence
 *    size (total item length).  Use [0 0] for no limit, [n 0] for a lower
 *    bound, [0 n] for an upper bound.
 *
 * Y: (n,1) real array of labels {-1, 1}.
 * weights: (n,1) array of reals in [0 ; 1].
 * tau: minimum required gain (if unsure, set to zero).
 *
 * pricing_count (optional): number of stumps to return (>= 1). (default: 1)
 * htype (optional): Boosting type, 2 for standard 2-class LPBoost, 1 for
 *    1.5-class LPBoost. (default: 2)
 *
 ****
 * Output
 *
 * subseqs: (1,p) cellarray of p sequences, format as S.
 * occurence: (1,p) cellarray of uint32 real arrays with the indices of
 *    sequences in S in which the subsequence appears.
 * ybase: (1,p) stump directions (1 for positive occurence, -1 for negative
 *    occurence).
 * GY: (n,p) stump response matrix for the training set.
 */

#include <assert.h>
#include <limits.h>
#include <mex.h>

#include <vector>
#include <iostream>

#include <Sequence.h>
#include <PrefixSpan.h>
#include <PrefixSpanLPBoost.h>
#include <PrefixSpanLPBoost15.h>


static mxArray*
writeSequence (const SeqMine::Sequence& seq)
{
	mxArray* elems = mxCreateCellMatrix (1, seq.lengthElements ());
	for (unsigned int eidx = 0 ; eidx < seq.lengthElements () ; ++eidx) {
		mxArray* elem = mxCreateNumericMatrix (1, seq.lengthSet (eidx),
			mxUINT32_CLASS, 0);
		uint32_T* elemP = (uint32_T*) mxGetPr (elem);

		for (unsigned int iidx = 0 ; iidx < seq.lengthSet (eidx) ; ++iidx)
			elemP[iidx] = seq.get (eidx, iidx);

		mxSetCell (elems, eidx, elem);
	}

	return (elems);
}


/* Simply modification to normal PrefixSpan class: record all frequent
 * subsequences.
 */
class
PrefixSpanRecord : public SeqMine::PrefixSpan
{
	void Report (const SeqMine::Sequence& seq)
	{
		std::cout << "PrefixSpanRecord::Report, with: " << seq << std::endl;

		if (seq.length () < options.minlen)
			return;

		collected.push_back (writeSequence (seq));
		collectedOccurence.push_back (seq.Occurence);
	}

public:
	std::vector<mxArray*> collected;
	std::vector<std::set<unsigned int> > collectedOccurence;

	PrefixSpanRecord (SeqMine::PrefixSpanOptions& options)
		: PrefixSpan (options)
	{
	}
};


static void
cleanupSequences (std::vector<const SeqMine::Sequence*>& cleanupS)
{
	for (unsigned int n = 0 ; n < cleanupS.size () ; ++n)
		delete cleanupS[n];

	cleanupS.clear ();
}


/* Parse a single sequence from Matlab's data structures into ours.
 */
static SeqMine::SequenceT*
parseSequence (const mxArray* singleSeq, unsigned int seqId)
{
	if (mxIsCell (singleSeq) == 0) {
		mexErrMsgTxt ("Input sequence is not a valid sequence.\n");

		return (NULL);
	}

	/* Create sequence object and fill it in, doh!
	 */
	SeqMine::SequenceT* res = new SeqMine::SequenceT (seqId);

	const mxArray** singleSeqP = (const mxArray**) mxGetPr (singleSeq);
	for (unsigned int eidx = 0 ; eidx < mxGetN (singleSeq) ; ++eidx) {
		unsigned int elementLength = mxGetN (singleSeqP[eidx]);
		std::vector<unsigned int> pvec (elementLength);

		if (mxIsUint32 (singleSeqP[eidx]) == 0) {
			mexErrMsgTxt ("Sequence items are not of type uint32, doh!\n");
			delete res;

			return (NULL);
		}

		uint32_T* elemP = (uint32_T*) mxGetPr (singleSeqP[eidx]);
		for (unsigned int iidx = 0 ; iidx < elementLength ; ++iidx)
			pvec[iidx] = elemP[iidx];

		res->appendElement (pvec);
	}

	return (res);
}


static void
handleOutput (int nlhs, mxArray* plhs[],
	const std::vector<double>& Y,
	const SeqMine::BestPatternList& bestPatterns, double negResponse)
{
	unsigned int num_patterns = bestPatterns.size ();

	/* Convert results back into Matlab format:
	 *
	 * plhs[0], subseqs: (1,p) cellarray
	 * plhs[1], ybase: stump directions
	 * plhs[2], GY: stump response
	 */
	mxArray* seqP = plhs[0] = mxCreateCellMatrix (1, num_patterns);

	plhs[1] = mxCreateNumericMatrix (1, num_patterns, mxDOUBLE_CLASS, 0);
	double* ybaseP = (double*) mxGetPr (plhs[1]);

	plhs[2] = mxCreateNumericMatrix (Y.size (), num_patterns, mxDOUBLE_CLASS, 0);
	double* GYP = (double*) mxGetPr (plhs[2]);

	/* Mining returned with some successful sequences that gain us
	 * something.  Lets take a look at those and where they appear.
	 */
	unsigned int sidx = 0;
	for (std::multimap<double, std::pair<SeqMine::SequenceT, double> >::const_iterator
		iv = bestPatterns.begin () ; iv != bestPatterns.end () ; ++iv)
	{
		/* Output subsequence itself
		 */
		mxSetCell (seqP, sidx, writeSequence (iv->second.first));

		/* Output ybase and GY column
		 */
		ybaseP[sidx] = iv->second.second;
		for (unsigned int n = 0 ; n < Y.size () ; ++n) {
			bool sampleActive = (iv->second.first.Occurence.find (n)
				!= iv->second.first.Occurence.end ());

			double response = sampleActive ? 1.0 : negResponse;
			response *= iv->second.second;

			GYP[sidx*Y.size()+n] = response;
		}

		sidx += 1;	// Next, please
	}
}


void
mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nlhs > 3) {
		mexErrMsgTxt ("Too many left hand side parameters.\n");
		return;
	}

	if (nrhs != 3 && (nrhs < 7 || nrhs > 8)) {
		mexErrMsgTxt ("Wrong number of right hand side parameters (must be 3, 7 or 8).\n");

		return;
	}

	if (mxIsCell (prhs[0]) == 0) {
		mexErrMsgTxt ("S must be a cellarray of graph structures.\n");

		return;
	}

	unsigned int minsup = (unsigned int) mxGetScalar (prhs[1]);
	if (mxGetM (prhs[2]) != 1 || mxGetN (prhs[2]) != 2) {
		mexErrMsgTxt ("size_req must be a (1,2) real vector.\n");

		return;
	}
	double* size_req = (double *) mxGetPr (prhs[2]);
	unsigned int maxpat_min = (unsigned int) size_req[0];
	unsigned int maxpat_max = (unsigned int) size_req[1];

	/* Parse in all sequences
	 */
	std::vector<const SeqMine::Sequence*> S;
	const mxArray** S_P = (const mxArray**) mxGetPr (prhs[0]);

	for (unsigned int n = 0 ; n < mxGetN (prhs[0]) ; ++n) {
		SeqMine::SequenceT* seq_cur = parseSequence (S_P[n], n);
		if (seq_cur == NULL) {
			cleanupSequences (S);

			return;
		}
		S.push_back (seq_cur);
	}

	SeqMine::PrefixSpanOptions options (minsup, maxpat_min);

	if (nrhs == 3) {
		/* Frequent subsequence mining
		 */
		PrefixSpanRecord pspan = PrefixSpanRecord (options);
		pspan.Mine (S);

		unsigned int num_patterns = pspan.collected.size ();
		mxArray* seqP = plhs[0] = mxCreateCellMatrix (1, num_patterns);
		mxArray* occP = plhs[1] = mxCreateCellMatrix (1, num_patterns);

		unsigned int sidx = 0;
		for (unsigned int n = 0 ; n < num_patterns ; ++n) {
			mxSetCell (seqP, sidx, pspan.collected[n]);

			mxArray* occCur = mxCreateNumericMatrix (1,
				pspan.collectedOccurence[n].size (), mxUINT32_CLASS, 0);
			uint32_T* occCurP = (uint32_T*) mxGetPr (occCur);

			unsigned int pidx = 0;
			for (std::set<unsigned int>::iterator si =
				pspan.collectedOccurence[n].begin () ;
				si != pspan.collectedOccurence[n].end () ; ++si)
			{
				occCurP[pidx] = *si + 1;	// Matlab index correction
				pidx += 1;
			}
			mxSetCell (occP, sidx, occCur);

			sidx += 1;
		}
	} else {
		/* Weighted subsequence mining
		 */
		assert (mxGetM (prhs[3]) == mxGetM (prhs[4]));
		assert (mxGetN (prhs[3]) == 1 && mxGetN (prhs[4]) == 1);

		/* Labels
		 */
		std::vector<double> Y (mxGetM (prhs[3]));
		const double* Y_P = (const double*) mxGetPr (prhs[3]);
		for (unsigned int n = 0 ; n < mxGetM (prhs[3]) ; ++n)
			Y[n] = Y_P[n];

		/* Weights
		 */
		std::vector<double> W (Y.size ());
		const double* W_P = (const double*) mxGetPr (prhs[4]);
		for (unsigned int n = 0 ; n < mxGetM (prhs[4]) ; ++n)
			W[n] = W_P[n];

		/* Minimum required gain
		 */
		double tau = mxGetScalar (prhs[5]);
		unsigned int pricingCount = 1;

		if (nrhs >= 7)
			pricingCount = (unsigned int) mxGetScalar (prhs[6]);
		assert (pricingCount > 0);

		unsigned int htype = 2;
		if (nrhs >= 8)
			htype = (unsigned int) mxGetScalar (prhs[7]);
		assert (htype == 1 || htype == 2);

		unsigned int maximumMiningIterations = 2000*1000;
		bool memoryConserving = true;

		/* Either call 1.5 or 2-class LPBoost.
		 */
		if (htype == 1) {
			SeqMine::PrefixSpanLPBoost15 boost (options, Y, W, tau, pricingCount,
				maximumMiningIterations);
			boost.Mine (S);

			handleOutput (nlhs, plhs, Y, boost.bestPatterns, 0.0);
		} else if (htype == 2) {
			SeqMine::PrefixSpanLPBoost boost (options, Y, W, tau, pricingCount,
				maximumMiningIterations);
			boost.Mine (S);

			handleOutput (nlhs, plhs, Y, boost.bestPatterns, -1.0);
		}
	}

	/* Clean up our allocated sequences
	 */
	cleanupSequences (S);
}

