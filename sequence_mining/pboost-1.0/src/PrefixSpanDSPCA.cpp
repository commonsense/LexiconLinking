
#include <algorithm>
#include <limits>

#include <PrefixSpanDSPCA.h>

namespace SeqMine
{

bool
PrefixSpanDSPCA::TestExtendSequence (const Sequence& seq) const
{
	// TODO
	return (seq.Occurence.size () >= options.minsup);
}

double
PrefixSpanDSPCA::ExplorationPriority (const Sequence& seq) const
{
	// TODO
	return (0);
}

double
PrefixSpanDSPCA::ExplorationPriorityBound (const Sequence& seq) const
{
	// TODO
	return (1);
}


void
PrefixSpanDSPCA::Report (const Sequence& seq)
{
	unsigned int freq = seq.Occurence.size ();

	if (freq < options.minsup)
		return;

	/* Check minimum length requirement.
	 */
	if (seq.length () < options.minlen)
		return;

	bool curVerbose = false;
	reportCount += 1;
	if (reportCount % 10000 == 1) {
		std::cout << "[" << reportCount << "] " << seq
			<< ", freq " << freq << ", minsup " << options.minsup << std::endl;

		curVerbose = true;
	}

	if (forbiddenSequences.count (seq) > 0) {
		//std::cout << "Skipping over forbidden sequence: " << seq << std::endl;

		return;
	}

	/* Allocate a new covariance matrix in case its not there yet
	 */
	if (Cv.size () == components) {
		Cv.push_back (gmm::dense_matrix<double> (0, 0));
		Cpatterns.push_back (component_patterns_type ());
	}

	/* The current covariance submatrix C.  In the following we construct a
	 * matrix C2 with
	 *   C2 = [C b ; b' d]
	 *
	 * and measure its maximum Eigenvalue \lambda_{max}(C2).  The vector b is
	 * the cross-variance with the C-patterns, d is its self-variance.
	 */
	gmm::dense_matrix<double>& C = Cv[components];
	component_patterns_type& Cpat = Cpatterns[components];

	gmm::dense_matrix<double> C2 (C);
	gmm::resize (C2, gmm::mat_nrows (C) + 1, gmm::mat_ncols (C) + 1);

	/* Compute cross-variance elements
	 */
	for (unsigned int n = 0 ; n < gmm::mat_nrows (C) ; ++n) {
		double bi = 0.0;

		std::set<unsigned int> shared;
		std::set_intersection (Cpat[n].Occurence.begin (),
			Cpat[n].Occurence.end (), seq.Occurence.begin (),
			seq.Occurence.end (), std::inserter (shared, shared.begin ()));

		// O_i \cap O_j
		bi += ((double) shared.size ());
		// -2/N |O_i| |O_j| + 1/N |O_i| |O_j|
		bi -= (1.0 / ((double) sequenceCount)) *
			Cpat[n].Occurence.size () * seq.Occurence.size ();
		bi /= (double) (sequenceCount - 1);

		C2(gmm::mat_nrows (C), n) = C2(n, gmm::mat_nrows (C)) = bi;
	}

	/* Compute self-variance
	 */
	double bn = 0.0;
	bn += ((double) seq.Occurence.size ());
	bn -= (1.0 / (double) sequenceCount) *
		seq.Occurence.size () * seq.Occurence.size ();
	bn /= (double) (sequenceCount - 1);

	C2(gmm::mat_nrows (C), gmm::mat_nrows (C)) = bn;

	/* Compute the Eigenvalues
	 */
	std::vector<double> eigval (gmm::mat_nrows (C2));
	gmm::symmetric_qr_algorithm (C2, eigval);

	if (curVerbose) {
		std::cout << "eval:";
		for (unsigned int n = 0 ; n < eigval.size () ; ++n)
			std::cout << " " << eigval[n];

		std::cout << std::endl;
	}

	double lambdaMax = *std::max_element (eigval.begin (), eigval.end ());
#ifdef	DEBUG
	std::cout << "lambdaMax = " << lambdaMax << ", bestSeqLambdaMax = "
		<< bestSeqLambdaMax << std::endl;
#endif

	if (lambdaMax > bestSeqLambdaMax) {
		foundPattern = true;

		bestSeqLambdaMax = lambdaMax;
		bestSeq = seq;
		bestC = C2;

		std::cout << "new best " << lambdaMax << ": " << seq << std::endl;
	}
}


void
PrefixSpanDSPCA::PCA (std::vector<const Sequence*>& S)
{
	forbiddenSequences.clear ();
	Cv.clear ();
	Cpatterns.clear ();
	Ceval.clear ();

	sequenceCount = S.size ();

	/* Get PCA basis iteratively
	 */
	for (components = 0 ; components < basisCount ; ++components) {
		/* Grow covariance matrix iteratively
		 */
		for (unsigned int el = 0 ; el < basisNonZero ; ++el) {
			bestSeqLambdaMax = -std::numeric_limits<double>::max ();
			gmm::resize (bestC, 0, 0);

			std::cout << "MINE" << std::endl;
			reportCount = 0;

			foundPattern = false;
			Mine3 (S);

			/* Check if we actually found an Eigenvalue-increasing pattern.
			 * If this is not the case, there might be no orthogonal patterns
			 * left.
			 */
			if (foundPattern == false) {
				std::cout << "Failed to grow basis vector, component finished." << std::endl;

				break;
			}

			Cpatterns[components].push_back (bestSeq);
			forbiddenSequences.insert (bestSeq);

			Cv[components] = bestC;
		}

		std::cerr << "gmm::mat_nrows (Cv[components]) = "
			<< gmm::mat_nrows (Cv[components]) << std::endl;
		std::cerr << "basisNonZero = " << basisNonZero << std::endl;
		std::cerr << "Cpatterns[components].size () = "
			<< Cpatterns[components].size () << std::endl;
		assert (foundPattern == false || gmm::mat_nrows (Cv[components]) == basisNonZero);

		if (gmm::mat_nrows (Cv[components]) == 0) {
			std::cout << "Failed to find any non-zero coefficient, PCA finished." << std::endl;

			break;
		}

		/* We a new basis set, compute the basis vector
		 */
		std::vector<double> eigval (gmm::mat_nrows (Cv[components]));
		gmm::dense_matrix<double> eigvect (gmm::mat_nrows (Cv[components]),
			gmm::mat_nrows (Cv[components]));
		gmm::symmetric_qr_algorithm (Cv[components], eigval, eigvect);

		int max_idx = -1;
		double eigvalMax = -std::numeric_limits<double>::max ();
		for (unsigned int i = 0 ; i < eigval.size () ; ++i) {
			if (eigval[i] <= eigvalMax)
				continue;

			eigvalMax = eigval[i];
			max_idx = i;
		}

		std::cout << "BASIS VECTOR with lambda " << eigvalMax << std::endl;
		std::vector<double> Ceval_cur (gmm::mat_nrows (Cv[components]));
		for (unsigned int el = 0 ; el < gmm::mat_nrows (Cv[components]) ; ++el) {
 			Ceval_cur[el] = eigvect(el,max_idx);

			std::cout << eigvect(el,max_idx) << ": "
				<< Cpatterns[components][el] << std::endl;
		}
		Ceval.push_back (Ceval_cur);

		std::cout << std::endl;
	}
}

}


