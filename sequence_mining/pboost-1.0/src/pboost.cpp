/* pboost Discriminative Subsequence Mining
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

/* pboost.cpp - Simple textual frontend to LPBoosting based on weighted
 * subsequence mining.
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 12th June 2007
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <set>

#include <LPBoostOSI.h>
#include <LPBoost15OSI.h>
#include <DDAG.h>
#include <Classifier.h>

#include <boost/program_options.hpp>


namespace po = boost::program_options;

/* Parse a single sequence file 'seqfilename'.
 */
static const SeqMine::Sequence*
parseSequenceFile (char* seqfilename, unsigned int seqId)
{
	std::cout << "new sequence" << std::endl;
	SeqMine::SequenceT* seq = new SeqMine::SequenceT (seqId);

	std::ifstream in (seqfilename);

	/* Read in the sequence: each line a bag, each bag a set of ordered
	 * indices.
	 */
	while (in.eof () == false) {
		//std::cout << "    new elem" << std::endl;
		char line[32768];

		line[0] = '\0';
		in.getline (line, sizeof (line) - 1);

		/* Skip over empty elements: they are not important to our sequence
		 * mining.
		 */
		if (strlen (line) == 0) {
			//std::cout << "        skipped." << std::endl;
			continue;
		}
		std::istringstream is (line);

		std::vector<unsigned int> elem;
		while (is.eof () == false) {
			unsigned int item;

			is >> item;
			elem.push_back (item);
			//std::cout << "        read item " << item << std::endl;
		}

		seq->appendElement (elem);
	}
	in.close ();

	return (seq);
}

/* Parse the input set description file 'filename', read all contained
 * sequence files and build the training set 'S', with associated labels
 * stored in 'Ylabels'.
 */
static bool
parseSequenceMetaFile (const char* filename,
	std::vector<const SeqMine::Sequence*>& S,
	std::vector<double>& Ylabels)
{
	std::ifstream in (filename);

	/* Read in pairs of filenames and labels
	 */
	S.clear ();
	Ylabels.clear ();

	for (unsigned int seqId = 0 ; in.eof () == false ; ++seqId) {
		char seqFilename[256];
		double label;

		in.width (sizeof (seqFilename) - 1);
		seqFilename[0] = '\0';
		in >> seqFilename;
		if (strlen (seqFilename) == 0)
			break;
		in >> label;
		std::cout << "  reading \"" << seqFilename << "\", label " << label << std::endl;

		Ylabels.push_back (label);
		S.push_back (parseSequenceFile (seqFilename, seqId));
	}
	in.close ();

	std::cout << "Successfully read in " << Ylabels.size ()
		<< " sequence files." << std::endl;

	return (true);
}


void
trainingSplit (std::vector<const SeqMine::Sequence*>& S,
	std::vector<double>& Ylabels, double trainPart,
	std::vector<const SeqMine::Sequence*>& Strain,
	std::vector<double>& Ytrain,
	std::vector<const SeqMine::Sequence*>& Sval,
	std::vector<double>& Yval)
{
	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
		double rd = (double) random () / ((double) RAND_MAX);

		if (rd < trainPart) {
			Strain.push_back (new SeqMine::SequenceT (*S[n],
				(unsigned int) Ytrain.size ()));
			Ytrain.push_back (Ylabels[n]);
		} else {
			Sval.push_back (new SeqMine::SequenceT (*S[n],
				(unsigned int) Yval.size ()));
			Yval.push_back (Ylabels[n]);
		}
	}
}


static void
ProduceMulticlassSplit (unsigned int pos_idx, unsigned int neg_idx,
	const std::vector<const SeqMine::Sequence*>& Smulti,
	const std::vector<double>& Ymulti,
	std::vector<const SeqMine::Sequence*>& Sout,
	std::vector<double>& Yout)
{
	Sout.clear ();
	Yout.clear ();

	unsigned int nsidx = 0;
	for (unsigned int n = 0 ; n < Ymulti.size () ; ++n) {
		std::cout << "Sample " << n << " label " << Ymulti[n] << std::endl;

		if (Ymulti[n] == ((double) pos_idx) || Ymulti[n] == ((double) neg_idx)) {
			std::cout << "   "
				<< ((Ymulti[n] == ((double) pos_idx)) ? "positive" : "negative")
				<< std::endl;

			Sout.push_back (new SeqMine::SequenceT (*Smulti[n], nsidx));
			Yout.push_back ((Ymulti[n] == ((double) pos_idx)) ? 1.0 : -1.0);
			nsidx += 1;
		}
	}
}


int
main (int argc, char* argv[])
{
	/* Parse command line options
	 */
	po::options_description generic ("Generic options");
	generic.add_options ()
		("help", "Produce help message")    
		;

	double nu;
	double convEpsilon;
	double alphaBoundC;
	unsigned int support_min;
	unsigned int pos_support_min;
	unsigned int length_min;
	unsigned int length_max;
	int maxgap;
	unsigned int search_strategy;
	unsigned int classifier_type;
	unsigned int pricing_count;
	unsigned int exploration_max;
	bool directional;
	std::string outputFilename;

	po::options_description config ("Parameters");
	config.add_options ()
		("output", po::value<std::string>(&outputFilename)->default_value ("output.txt"),
			"Classifier output file")
		("classifier-type",
			po::value<unsigned int>(&classifier_type)->default_value (2),
			"Classifier type, 1: 1.5-class LPBoost, 2: 2-class LPBoost")
		("nu",
			po::value<double>(&nu)->default_value (0.1),
			"LPBoost nu parameter, 0.0 < nu < 1.0")
		("conv-eps", po::value<double>(&convEpsilon)->default_value (1e-3),
			"LPBoost convergence epsilon")
		("pricing-count",
			po::value<unsigned int>(&pricing_count)->default_value (1),
			"Number of patterns to add in each iteration")
		("search-strategy",
			po::value<unsigned int>(&search_strategy)->default_value (3),
			"Search strategy, 1: A-Star, 2: DFS, 3: id A-Star")
		("exploration-max",
			po::value<unsigned int>(&exploration_max)->default_value (1000*1000),
			"Maximum number of examined patterns in each iteration or "
				"zero for exhaustive search")
			// TODO: replace environment variables
			;

	po::options_description configReg ("Regularization options");
	configReg.add_options ()
		("support-min,S",
			po::value<unsigned int>(&support_min)->default_value (2),
			"Minimum support of the pattern in the entire training set")
		("length-min,l",
			po::value<unsigned int>(&length_min)->default_value (1),
			"Minimum pattern sequence length")
		("length-max,L",
			po::value<unsigned int>(&length_max)->default_value (0),
			"Maximum pattern sequence length or zero for no maximum")
		("maxgap,G", po::value<int>(&maxgap)->default_value (-1),
			"Maximum allowed gap or -1 for infinite gaps")
			;

	po::options_description config15 ("1.5-class LPBoost options");
	config15.add_options ()
		("alpha-bound", po::value<double>(&alphaBoundC)->default_value (1.0),
			"alpha_i <= C.  So far only implemented for 1.5-class boosting")
		("pos-support-min",
			po::value<unsigned int>(&pos_support_min)->default_value (0),
			"Minimum support of the pattern in the positive-labelled "
				"training set (only for 1.5-class LPBoost)")
		;

	po::options_description config2 ("2-class LPBoost options");
	config2.add_options ()
		("directional",
			po::value<bool>(&directional)->default_value (false),
			"Use one-sided directional (omega/0)-stumps")
		;

	po::options_description hidden ("Hidden options");
	hidden.add_options ()
		("input", po::value<std::vector<std::string> >(), "Sequence input file")
		;

	po::options_description visible_options;
	visible_options.add (generic).add (config).add (configReg).add (config2).add (config15);

	po::options_description all_options;
	all_options.add (config).add (generic).add (configReg).add (config2).add (config15).add (hidden);

	po::positional_options_description p;
	p.add ("input", -1);

	/* Parse
	 */
	po::variables_map vm;
	po::store (po::command_line_parser (argc, argv).options (all_options).positional (p).run (), vm);
	po::notify (vm);

	if (vm.count ("help") || vm.count ("input") == 0 ||
		vm["input"].as<std::vector<std::string> >().size () != 1)
	{
		std::cerr << "pboost $Id: pboost.cpp 1097 2007-09-24 09:17:53Z nowozin $" << std::endl;
		std::cerr << "================================================================================" << std::endl;
		std::cerr << "Copyright (C) 2007 -- Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>" << std::endl;
		std::cerr << std::endl;
		std::cerr << "Usage: pboost [options] input.txt" << std::endl;
		std::cerr << std::endl;
		std::cerr << "Required parameters:" << std::endl;
		std::cerr << "  input.txt     Training sequences" << std::endl;
		std::cerr << visible_options << std::endl;

		std::cerr << "Input files:" << std::endl;
		std::cerr << "  input.txt contains lines \"sequence1.txt 1\", where 1/-1 is the class label." << std::endl;
		std::cerr << "      If the class labels are 0,...,K then a multi-class problem is solved by" << std::endl;
		std::cerr << "      DDAG decomposition into 1-vs-1 classifiers." << std::endl;
		std::cerr << std::endl;

		std::cerr << "Environment variables:" << std::endl;
		std::cerr << "  PBOOST_TRAIN_GY, PBOOST_TEST_GY, PBOOST_ALPHA environment variables "
			<< std::endl;
		std::cerr << "      can be set to filenames such that output will be written to them."
			<< std::endl;
		std::cerr << "  PBOOST_CLASSIFIER will produce a sequence meta file of the final classifier." 
			<< std::endl;
		std::cerr << std::endl;

		exit (EXIT_FAILURE);
	}

	assert (nu > 0.0 && nu < 1.0);
	assert (support_min >= 1);
	assert (alphaBoundC >= 1.0 || classifier_type == 1);

	if (classifier_type != 1 && pos_support_min != 0) {
		std::cerr << "Minimum support for positive class is available only with 1.5-class LPBoost." << std::endl;

		exit (EXIT_FAILURE);
	}

	/* Parse in the input data file.
	 */
	std::vector<const SeqMine::Sequence*> S;
	std::vector<double> Ylabels;

	std::vector<std::string> input = vm["input"].as<std::vector<std::string> >();
	bool parseResult = parseSequenceMetaFile (input[0].c_str (), S, Ylabels);
	assert (parseResult);

	/* Check if we are solving a multi-class problem
	 */
	bool multiclass = false;
	bool seenNegative = false;
	std::set<int> seenLabels;

	for (unsigned int n = 0 ; n < Ylabels.size () ; ++n) {
		if (Ylabels[n] < 0)
			seenNegative = true;

		if (Ylabels[n] != 1.0 && Ylabels[n] != -1.0)
			multiclass = true;

		seenLabels.insert ((int) Ylabels[n]);
	}

	/* Build the decision dag
	 */
	unsigned int classCount = seenLabels.size ();	// Number of classes

	/* Multiclass input/output label translation tables for DDAG.
	 *
	 * Map from input class label to internal class label, use
	 * mclassInputLabel; for the reverse direction, use find.
	 */
	std::vector<int> mclassOutputLabel;	// output label, eg. [1 -1] or [0 1 2]

	if (multiclass) {
		if (seenNegative) {
			std::cerr << "   OOPS, there are negative class labels, please see the usage." << std::endl;
			exit (EXIT_FAILURE);
		}

		std::cerr << "Solving multi-class problem, "
			<< classCount << " classes" << std::endl;

		for (unsigned int n = 0 ; n < classCount ; ++n)
			mclassOutputLabel.push_back (n);	// identity
	} else {
		classCount = 2;
		std::cerr << "Solving two-class problem" << std::endl;

		mclassOutputLabel.push_back (1);	// positive class gets mapped to index 0
		mclassOutputLabel.push_back (-1);	// negative class gets mapped to index 1
	}
	Classifier::DDAG ddag = Classifier::DDAG (classCount);
	std::cerr << "Maximum gap constraint parameter: " << maxgap << std::endl;

	std::vector<const SeqMine::Sequence*> Strain = S;
	std::vector<double> Ytrain = Ylabels;


	/* Build the LPBoost classifier and train it.
	 */
	SeqMine::PrefixSpanLPBoost::SearchStrategy searchStrategy =
		SeqMine::PrefixSpanLPBoost::DepthFirstSearch;

	searchStrategy = (SeqMine::PrefixSpanLPBoost::SearchStrategy) search_strategy;
	std::cerr << "Selected search strategy: " << searchStrategy << std::endl;

	assert (classifier_type == 1 || classifier_type == 2);
	assert (pricing_count > 0 && pricing_count < 1000);

	/* The corresponding 1-vs-1 classifiers
	 */
	unsigned int trun = 0;
	std::vector<Classifier::TwoClassClassifier*> resultClassifiers;

	for (std::vector<std::pair<int, int> >::iterator classIV = ddag.table.begin () ;
		classIV != ddag.table.end () ; ++classIV)
	{
		/* Skip class result targets
		 */
		if (classIV->first == classIV->second) {
			resultClassifiers.push_back (NULL);

			continue;
		}

		trun += 1;

		/* Produce the split for this 1-vs-1 class decomposition.
		 */
		std::vector<const SeqMine::Sequence*> Ssplit;
		std::vector<double> Ysplit;

		ProduceMulticlassSplit (classIV->first, classIV->second,
			Strain, Ytrain, Ssplit, Ysplit);

		std::cerr << "1-vs-1 training "
			<< trun << "/" << ((ddag.K*(ddag.K-1)) / 2) << ", class "
			<< classIV->first << " (pos, "
			<< std::count (Ysplit.begin (), Ysplit.end (), 1.0) << ") vs. "
			<< classIV->second << " (neg, "
			<< std::count (Ysplit.begin (), Ysplit.end (), -1.0) << ")"
			<< std::endl;

		/* Results are put here
		 */
		std::vector<SeqMine::SequenceT*> currentStumpSequences;
		std::vector<double> currentStumpAlpha;
		double currentClassThresh = 0.0;

		SeqMine::PrefixSpanOptions poptions (support_min, length_min, length_max, maxgap);
		if (classifier_type == 2) {
			/* 2-class LPBoost
			 */
			LPBoost::LPBoostOSI boost (poptions, nu, pricing_count,
				convEpsilon, exploration_max);

			/* Use one-sided directional stumps?
			 */
			if (directional)
				boost.directional = true;

			//boost.SetValidationSet (Sval, Yval);
			boost.BuildClassifier (Ssplit, Ysplit, searchStrategy);

			/* Output stumps
			 */
			std::cout << "Final decision stumps and their weights:" << std::endl;
			for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
				std::cout << "ybase " << boost.stumps[n].second
					<< " weight " << boost.alpha[n] << ": "
					<< boost.stumps[n].first << std::endl;
			}

			std::cout << "GY indicator matrix, training set (prefixspan knowledge)" << std::endl;
			boost.DumpGY (Ssplit, Ysplit, 1);
			std::cout << std::endl;

			std::cout << "GY indicator matrix, training set (contains)" << std::endl;
			boost.DumpGY (Ssplit, Ysplit, 0);
			std::cout << std::endl;

			for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
				if (n > 0)
					std::cout << " ";

				std::cout << boost.alpha[n];
			}
			std::cout << std::endl;

			/* Store results
			 */
			currentClassThresh = 0.0;
			for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
				if (boost.alpha[n] < 1e-8)
					continue;

				currentStumpAlpha.push_back (boost.stumps[n].second * boost.alpha[n]);
				currentStumpSequences.push_back (new SeqMine::SequenceT (boost.stumps[n].first));
			}

#if 0
			if (getenv ("PBOOST_CLASSIFIER") != NULL) {
				std::ofstream pout (getenv ("PBOOST_CLASSIFIER"));

				pout << "2" << std::endl;	// 2-class classifier
				pout << (directional ? "1" : "2") << std::endl;	// directionality
				for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
					char subSeqFilename[256];
					snprintf (subSeqFilename, sizeof (subSeqFilename) - 1, "%s.seq%d",
						getenv ("PBOOST_CLASSIFIER"), n);

					std::cout << "writing file " << subSeqFilename << std::endl;

					pout << subSeqFilename << " "
						<< (boost.stumps[n].second * boost.alpha[n])
						<< std::endl;

					boost.stumps[n].first.writeSequenceFile (subSeqFilename);
				}
				pout << std::endl;
				pout.close ();
			}
#endif
		} else if (classifier_type == 1) {
			/* 1.5-class LPBoost
			 */
			LPBoost::LPBoost15OSI boost (poptions, nu, pricing_count, convEpsilon,
				pos_support_min, exploration_max, alphaBoundC);
			boost.BuildClassifier (Ssplit, Ysplit, searchStrategy);

			/* Output stumps
			 */
			std::cout << "Final decision stumps and their weights:" << std::endl;
			for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
				std::cout << " weight " << boost.alpha[n] << ": "
					<< boost.stumps[n].first << std::endl;
			}

			if (getenv ("PBOOST_TRAIN_GY") != NULL) {
				std::ofstream pout (getenv ("PBOOST_TRAIN_GY"));
				boost.DumpGY (pout, Ssplit, Ysplit, 0);
				pout.close ();
			} else {
				std::cout << "GY indicator matrix, training set (contains)" << std::endl;
				boost.DumpGY (std::cout, Ssplit, Ysplit, 0);
				std::cout << std::endl;
			}

			if (getenv ("PBOOST_ALPHA") != NULL) {
				std::ofstream pout (getenv ("PBOOST_ALPHA"));
				for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
					if (n > 0)
						pout << " ";

					pout << boost.alpha[n];
				}
				pout << std::endl;
				pout.close ();
			} else {
				for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
					if (n > 0)
						std::cout << " ";

					std::cout << boost.alpha[n];
				}
				std::cout << std::endl;
			}

			/* Store results
			 */
			currentClassThresh = (boost.rho1 + boost.rho2) / 2.0;
			for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
				if (boost.alpha[n] < 1e-8)
					continue;

				currentStumpAlpha.push_back (boost.stumps[n].second * boost.alpha[n]);
				currentStumpSequences.push_back (new SeqMine::SequenceT (boost.stumps[n].first));
			}

#if 0
			if (getenv ("PBOOST_CLASSIFIER") != NULL) {
				std::ofstream pout (getenv ("PBOOST_CLASSIFIER"));

				/* Output class decision threshold
				 */
				pout << "1" << std::endl;	// 1.5 class classifier
				pout << ((boost.rho1 + boost.rho2) / 2.0) << std::endl;

				for (unsigned int n = 0 ; n < boost.alpha.size () ; ++n) {
					char subSeqFilename[256];
					snprintf (subSeqFilename, sizeof (subSeqFilename) - 1, "%s.seq%d",
						getenv ("PBOOST_CLASSIFIER"), n);

					pout << subSeqFilename << " "
						<< (boost.stumps[n].second * boost.alpha[n])
						<< std::endl;

					boost.stumps[n].first.writeSequenceFile (subSeqFilename);
				}
				pout << std::endl;
				pout.close ();
			}
#endif
		}

		/* Store the resulting classifier
		 */
		resultClassifiers.push_back (new Classifier::TwoClassClassifier
			(currentStumpSequences, currentStumpAlpha, currentClassThresh,
				classifier_type, directional, maxgap));
	}

	// To trigger flush/close by destructors
	{
		std::ofstream ofs (outputFilename.c_str ());
		boost::archive::text_oarchive oa (ofs);

		Classifier::MulticlassClassifier mclassifier (ddag, resultClassifiers);
		oa << mclassifier;
#if 0
		if (resultClassifiers.size () == 1) {
			oa << *resultClassifiers[0];
		} else {
			/* Multiclass classifier
			 */
			Classifier::MulticlassClassifier mclassifier (ddag, resultClassifiers);
			oa << mclassifier;
		}
#endif
	}
}

