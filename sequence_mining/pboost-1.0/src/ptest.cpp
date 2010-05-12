/* ptest.cpp - Test sequences using a learnt classifier
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 19th March 2007
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <stdlib.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/program_options.hpp>

#include <Sequence.h>
#include <DDAG.h>
#include <Classifier.h>

namespace po = boost::program_options;

/* Parse a single sequence file 'seqfilename'.
 */
static SeqMine::Sequence*
parseSequenceFile (char* seqfilename)
{
	SeqMine::SequenceT* seq = new SeqMine::SequenceT (0);

	std::cerr << "opening " << seqfilename << std::endl;
	std::ifstream in (seqfilename);
	assert (in.fail () == false);

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
 * sequence files and build the training set 'S', with associated weights
 * alpha stored in 'alpha'.
 */
static bool
parseClassifierMetaFile (char* filename, std::vector<SeqMine::Sequence*>& S,
	std::vector<double>& alphaV, unsigned int& cltype, double& classthresh,
	bool& directional)
{
	std::ifstream in (filename);
	assert (in.fail () == false);

	/* A real classifier, opposed to a test sequence list
	 */
	if (cltype == 1) {
		in >> cltype;
		assert (cltype == 1 || cltype == 2);

		directional = false;
		classthresh = 0.0;
		if (cltype == 1) {
			in >> classthresh;
		} else {
			unsigned int dirnum;
			in >> dirnum;
			assert (dirnum >= 1 && dirnum <= 2);

			if (dirnum == 1)
				directional = true;
		}
	}

	/* Read in pairs of filenames and labels
	 */
	S.clear ();
	alphaV.clear ();

	for (unsigned int seqId = 0 ; in.eof () == false ; ++seqId) {
		char seqFilename[256];
		double alpha;

		in.width (sizeof (seqFilename) - 1);
		seqFilename[0] = '\0';
		in >> seqFilename;
		if (strlen (seqFilename) == 0)
			break;
		in >> alpha;
		std::cout << "  reading \"" << seqFilename
			<< "\", weight " << alpha << std::endl;

		alphaV.push_back (alpha);
		S.push_back (parseSequenceFile (seqFilename));
	}
	in.close ();

	std::cout << "Successfully read in " << alphaV.size ()
		<< " sequence files." << std::endl;

	return (true);
}


int
main (int argc, char* argv[])
{
	std::string confusionFilename;
	std::string dumpClassifierFilename;
	int sequenceMatchesSample;	// ID of sample to verbosely print matches

	po::options_description config ("Parameters");
	config.add_options ()
		("confusion-matrix",
			po::value<std::string>(&confusionFilename),
			"Output confusion matrix")
		("output-seq-matches",
			po::value<int>(&sequenceMatchesSample)->default_value (-1),
			"Output sequence matches for a single sample")
		("output-match-matrix,O",
			po::value<std::vector<unsigned int> >(),
			"Select classifier and stump to dump")
		("dump-classifier",
			po::value<std::string>(&dumpClassifierFilename),
			"Dump a textual description of the classifier")
		;

	po::options_description hidden ("Hidden options");
	hidden.add_options ()
		("input", po::value<std::vector<std::string> >(), "input files")
		;

	po::options_description visible_options;
	visible_options.add (config);

	po::options_description all_options;
	all_options.add (config).add (hidden);

	po::positional_options_description p;
	p.add ("input", -1);

	po::variables_map vm;
	po::store (po::command_line_parser (argc, argv).options
		(all_options).positional (p).run (), vm);
	po::notify (vm);

	if (vm.count ("help") || vm.count ("input") == 0 ||
		vm["input"].as<std::vector<std::string> >().size () != 2)
	{
		std::cerr << "ptest $Id: ptest.cpp 1087 2007-09-24 08:21:00Z nowozin $" << std::endl;
		std::cerr << "================================================================================" << std::endl;
		std::cerr << "Copyright (C) 2007 -- Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>" << std::endl;
		std::cerr << std::endl;
		std::cerr << "Usage: ptest [options] serialized-classifier.txt testmetafile.txt" << std::endl;
		std::cerr << std::endl;
		std::cerr << "Required parameters:" << std::endl;
		std::cerr << "  serialized-classifier.txt    A classifier file produced with pboost."
			<< std::endl;
		std::cerr << "  testmetafile.txt             A meta file with test sequences and true labels"
			<< std::endl;
		std::cerr << visible_options << std::endl;
		std::cerr << std::endl;

		exit (EXIT_FAILURE);
	}

	std::vector<std::string> input = vm["input"].as<std::vector<std::string> >();

	/* Read in classifier
	 */
	Classifier::MulticlassClassifier mclassifier;
    {
        std::ifstream ifs (input[0].c_str (), std::ios::binary);
        boost::archive::text_iarchive ia (ifs);

        ia >> mclassifier;
    }

	/* Write out textual description of the classifier
	 */
	if (vm.count ("dump-classifier")) {
		std::ofstream clout (dumpClassifierFilename.c_str ());

		clout << "# DDAG 1-vs-1 decomposition, " << mclassifier.ddag.K << " classes." << std::endl;
		clout << "# maximum gap constraint: " << mclassifier.classifiers[0]->maxgap << std::endl;

		for (unsigned int tidx = 0 ; tidx < mclassifier.ddag.table.size () ; ++tidx) {
			if (mclassifier.ddag.table[tidx].first == mclassifier.ddag.table[tidx].second)
				continue;

			clout << "### class " << mclassifier.ddag.table[tidx].first
				<< " vs class " << mclassifier.ddag.table[tidx].second << std::endl;

			clout << "# " << tidx << ". "
				<< mclassifier.classifiers[tidx]->cltype << "-class classifier, "
				<< mclassifier.classifiers[tidx]->alpha.size () << " stumps" << std::endl;

			/* Dump individual stumps
			 */
			for (unsigned int sidx = 0 ; sidx < mclassifier.classifiers[tidx]->alpha.size () ; ++sidx) {
				clout << sidx << " " << mclassifier.classifiers[tidx]->alpha[sidx]
					<< " " << *(mclassifier.classifiers[tidx]->S[sidx]) << std::endl;
			}
		}
		clout.close ();
	}

	/* Read in test sequences
	 */
	std::vector<SeqMine::Sequence*> Stest;
	std::vector<double> SY;
	double dummy;
	bool bdummy;
	unsigned int tcltype = 0;
	parseClassifierMetaFile ((char *) input[1].c_str (), Stest, SY, tcltype, dummy, bdummy);

	/* Optionally build a output match matrix
	 */
	std::vector<std::vector<int> > OMM;	// indexed: [sample][tbin]
	std::ofstream OMMfile;
	std::vector<unsigned int> clid_sid;
	if (vm.count ("output-match-matrix")) {
		OMMfile.open ("omm.txt");

		clid_sid = vm["output-match-matrix"].as<std::vector<unsigned int> >();

		/* Ensure safety.
		 */
		assert (clid_sid[0] >= 0 && clid_sid[0] < mclassifier.ddag.table.size ());
		assert (mclassifier.ddag.table[clid_sid[0]].first !=
			mclassifier.ddag.table[clid_sid[0]].second);
		assert (clid_sid[1] >= 0 &&
			clid_sid[1] < mclassifier.classifiers[clid_sid[0]]->alpha.size ());

		std::cout << "clid " << clid_sid[0] << ", sid " << clid_sid[1] << std::endl;

		/* Allocate an all-zero matrix.
		 */
		unsigned int max_elength = 0;
		for (unsigned int n = 0 ; n < Stest.size () ; ++n) {
			unsigned int cur_elength = Stest[n]->lengthElements ();

			if (cur_elength > max_elength)
				max_elength = cur_elength;
		}

		for (unsigned int n = 0 ; n < Stest.size () ; ++n)
			OMM.push_back (std::vector<int> (max_elength));
	}

	/* Verbose sequence matches
	 */
	std::ofstream SSmatches;
	if (vm.count ("output-seq-matches"))
		SSmatches.open ("sequence-matches.txt");

	if (mclassifier.classifiers[0]->maxgap >= 0)
		std::cout << "# maximum gap constraint: " << mclassifier.classifiers[0]->maxgap
			<< std::endl;

	/* Test
	 */
	std::vector<std::vector<unsigned int> > confusionMatrix (mclassifier.ddag.K);
	for (unsigned int cl = 0 ; cl < mclassifier.ddag.K ; ++cl)
		confusionMatrix[cl] = std::vector<unsigned int> (mclassifier.ddag.K);

	double loss = 0.0;
	for (int n = 0 ; ((unsigned int) n) < SY.size () ; ++n) {
		int curY = 0;

		/* Walk down the DDAG
		 */
		unsigned int tidx = 0;
		do {
			if (mclassifier.ddag.table[tidx].first == mclassifier.ddag.table[tidx].second) {
				curY = mclassifier.ddag.table[tidx].first;

				break;
			}

			const Classifier::TwoClassClassifier* current = mclassifier.classifiers[tidx];
			double sum = 0.0;

			for (int j = 0 ; ((unsigned int) j) < current->alpha.size () ; ++j) {
				double bval = 0.0;

				/* Precise, sub-sequence matching
				 */
				std::vector<unsigned int> pmatches;
				bool matches;
				if (current->maxgap < 0)
					matches = Stest[n]->contains (*(current->S[j]), &pmatches);
				else
					matches = Stest[n]->containsGap (*(current->S[j]), current->maxgap);

				/* Should we collect the individual subsequence matches?
				 */
				if (vm.count ("output-seq-matches") && sequenceMatchesSample == n) {
					/* Sequence matches
					 */
					SSmatches << mclassifier.ddag.table[tidx].first << " "
						<< mclassifier.ddag.table[tidx].second << " "
						<< current->alpha[j] << " " << (matches ? "" : "not-")
						<< "match " << *(current->S[j]);

					/* Optionally also output matching sequence element indices
					 */
					if (matches) {
						SSmatches << " at";
						for (std::vector<unsigned int>::const_iterator pi = pmatches.begin () ;
							pi != pmatches.end () ; ++pi)
						{
							SSmatches << " " << *pi;
						}
					}

					SSmatches << std::endl;
				}

				/* Optionally store in the output match matrix
				 */
				if (vm.count ("output-match-matrix") &&
					tidx == clid_sid[0] && ((unsigned int) j) == clid_sid[1])
				{
					/* Directional indicator function
					 */
					for (std::vector<unsigned int>::const_iterator civ = pmatches.begin () ;
						civ != pmatches.end () ; ++civ)
					{
						OMM[n][*civ] = (matches ? 1 : 0) * (current->alpha[j] > 0 ? 1 : -1);
					}
				}

				bool pb;
				if (current->maxgap < 0)
					pb = Stest[n]->contains (*(current->S[j]));
				else
					pb = Stest[n]->containsGap (*(current->S[j]), current->maxgap);

				if (current->cltype == 2) {
					if (current->directional)
						bval = pb ? 1.0 : 0.0;
					else
						bval = pb ? 1.0 : -1.0;
				} else if (current->cltype == 1) {
					bval = pb ? 1.0 : 0.0;
				}

				bval *= current->alpha[j];
				sum += bval;
			}

			int classresult = (sum >= current->classthresh) ? 1 : -1;

			/* Output 1-vs-1 classifier result
			 */
			std::cout << "   [" << mclassifier.ddag.table[tidx].first
				<< " vs. " << mclassifier.ddag.table[tidx].second << "] ("
				<< sum << ">=" << current->classthresh << ") => "
				<< classresult << std::endl;

			/* Advance in DDAG
			 */
			tidx += mclassifier.ddag.relativeOffset (
				mclassifier.ddag.table[tidx].first,
				mclassifier.ddag.table[tidx].second, classresult);
		} while (1);

		/* curY contains the class index.
		 */
		bool correct = (curY == ((int) SY[n])) ? true : false;
		std::cout << "Sample " << (n+1) << "/" << SY.size () << ": "
			<< (correct ? "CORRECT" : "INCORRECT") << std::endl;

		loss += correct ? 0.0 : 1.0;
		confusionMatrix[(int) SY[n]][curY] += 1;
	}

	/* Close the verbose subsequence match file.
	 */
	if (sequenceMatchesSample >= 0) {
		SSmatches.close ();

		std::cout << "Subsequence matches for sample " << sequenceMatchesSample
			<< " written to \"sequence-matches.txt\"." << std::endl;
	}

	/* Write out the output match matrix
	 */
	if (vm.count ("output-match-matrix")) {
		for (std::vector<std::vector<int> >::const_iterator iv = OMM.begin () ;
			iv != OMM.end () ; ++iv)
		{
			if (iv != OMM.begin ())
				OMMfile << std::endl;

			for (std::vector<int>::const_iterator riv = iv->begin () ;
				riv != iv->end () ; ++riv)
			{
				if (riv != iv->begin ())
					OMMfile << " ";

				OMMfile << *riv;
			}
		}
		OMMfile.close ();
	}

	/* Print the overall test sample loss
	 */
	std::cout << "LOSS (" << SY.size () << " samples): "
		<< (loss/((double) SY.size ())) << std::endl;

	/* Write the confusion matrix
	 */
	if (vm.count ("confusion-matrix")) {
		std::cout << "Writing confusion matrix to \"" << confusionFilename.c_str ()
			<< "\"" << std::endl;

		std::ofstream CMout (confusionFilename.c_str ());
		for (unsigned int r = 0 ; r < mclassifier.ddag.K ; ++r) {
			for (unsigned int c = 0 ; c < mclassifier.ddag.K ; ++c) {
				if (c > 0)
					CMout << " ";

				CMout << confusionMatrix[r][c];
			}
			CMout << std::endl;
		}
		CMout.close ();
	}
}


