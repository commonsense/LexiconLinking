/* pspan Subsequence Mining
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

/* pspan.cpp - PrefixSpan frequent subsequence mining
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 4th July 2007
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iterator>

#include <boost/program_options.hpp>

#include <PrefixSpan.h>
#include <PrefixSpanK.h>
#include <PrefixSpanDSPCA.h>

namespace po = boost::program_options;

using namespace SeqMine;


static bool
parseSequenceMetaFile (const char* filename, std::vector<const SeqMine::Sequence*>& S)
{
	std::ifstream in (filename);
	
	/* Read in filenames
	 */
	S.clear ();
	


	for (unsigned int seqId = 0 ; in.eof () == false ; ++seqId) {
		char seqFilename[256];

		in.width (sizeof (seqFilename) - 1);
		seqFilename[0] = '\0';
		in >> seqFilename;
		if (strlen (seqFilename) == 0)
			break;

		std::cout << "  reading \"" << seqFilename << "\"" << std::endl;
		std::cout << "done 1" << std::endl;
		SequenceT* seq = SequenceT::readFromFile (seqFilename);
		std::cout << "done 1" << std::endl;
		seq->Occurence.insert (seqId);
		std::cout << "done 1" << std::endl;
		S.push_back (seq);
	}
	in.close ();

	std::cout << "Successfully read in " << S.size ()
		<< " sequence files." << std::endl;

	return (true);
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

	/* Output file
	 */
	std::string outputFilename;
	std::string projectFilename;
	unsigned int topK;
	unsigned int dspcaBasisNonZero;
	unsigned int dspcaBasisCount;
	unsigned int minsup = 0;
	double minsupPercent;
	unsigned int length_min, length_max;
	int maxgap;

	po::options_description config ("Parameters");
	config.add_options ()
		("output", po::value<std::string>(&outputFilename)->default_value ("output.txt"),
			"Frequent subsequence output file")
		("project,p", po::value<std::string>(&projectFilename),
			"Project samples on found sequences.")
		("verbose,v", "Verbose mode, dumping found subsequences.")
		("verify,V", "Verify all projections explictly (slow).")
		("topk,K", po::value<unsigned int>(&topK)->default_value (20),
			"Mine the top K frequent subsequences.")
		("dspca", po::value<unsigned int>(&dspcaBasisNonZero),
			"Perform direct sparse PCA with the given number of "
			"nonzero coefficients per basis vector")
		("dspcacount", po::value<unsigned int>(&dspcaBasisCount)->default_value (5),
			"The number of PCA basis vectors to extract")
		("minsup,S", po::value<unsigned int>(&minsup), "Minimum support (count)")
		("minsuppct,s", po::value<double>(&minsupPercent),
			"Minimum support tau, 0<tau<=1 (percentile)")
		("length-min",
			po::value<unsigned int>(&length_min)->default_value (1),
			"Minimum pattern sequence length")
		("length-max",
			po::value<unsigned int>(&length_max)->default_value (0),
			"Maximum pattern sequence length or zero for no maximum")
		("maxgap,G", po::value<int>(&maxgap)->default_value (-1),
			"Maximum allowed gap or -1 for infinite gaps")
		;

	po::options_description hidden ("Hidden options");
	hidden.add_options ()
		("input", po::value<std::vector<std::string> >(), "Sequence input file")
		;

	po::options_description visible_options;
	visible_options.add (generic).add (config);

	po::options_description all_options;
	all_options.add (config).add (generic).add (hidden);

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
		std::cerr << "pspan, $Id: pspan.cpp 1097 2007-09-24 09:17:53Z nowozin $" << std::endl;
		std::cerr << "================================================================================" << std::endl;
		std::cerr << "Copyright (C) 2007 -- Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>" << std::endl;
		std::cerr << std::endl;
		std::cerr << "Usage: pspan [options] input.txt" << std::endl;
		std::cerr << std::endl;
		std::cerr << "Required parameters:" << std::endl;
		std::cerr << "  input.txt           List of input sequences." << std::endl;
		std::cerr << visible_options << std::endl;
		std::cerr << std::endl;

		exit (EXIT_FAILURE);
	}

	/* Read sequence files
	 */
	std::vector<std::string> input = vm["input"].as<std::vector<std::string> >();
	std::vector<const SeqMine::Sequence*> S;

	bool seqresult = parseSequenceMetaFile (input[0].c_str (), S);
	if (seqresult == false) {
		std::cerr << "Error while reading in sequence files." << std::endl;

		exit (EXIT_FAILURE);
	}

	/* Get mode and calculate details
	 */
	if (vm.count ("minsup") && vm.count ("minsuppct")) {
		std::cerr << "The \"minsup\" and \"minsuppct\" options "
			"are mutually exclusive." << std::endl;

		exit (EXIT_FAILURE);
	}

	bool minsup_mode = false;
	if (vm.count ("minsup") || vm.count ("minsuppct"))
		minsup_mode = true;

	if (vm.count ("minsuppct")) {
		minsup = (unsigned int) (minsupPercent * ((double) S.size ()) + 0.5);
		std::cerr << "Minimum support is " << minsupPercent << ", set to "
			<< minsup << std::endl;
	}

	std::vector<std::pair<unsigned int, SequenceT> > subseq;
	PrefixSpanOptions poptions (minsup, length_min, length_max, maxgap);

	if (vm.count ("dspca") == 0) {
		if (minsup_mode) {
			PrefixSpan pspan (poptions);
			if (vm.count ("verify"))
				pspan.SetVerification (true);

			pspan.Mine (S);

			std::copy (pspan.subsequences.begin (), pspan.subsequences.end (),
				std::back_insert_iterator<std::vector<std::pair<unsigned int, SequenceT> > >
					(subseq));
		} else if (vm.count ("topk")) {
			PrefixSpanK pspan (topK, poptions);
			if (vm.count ("verify"))
				pspan.SetVerification (true);

			pspan.Mine (S);

			std::cout << "The minimum support threshold for K=" << topK
				<< " is " << pspan.MinimumSupportThreshold () << std::endl;

			std::copy (pspan.mostFrequent.begin (), pspan.mostFrequent.end (),
				std::back_insert_iterator<std::vector<std::pair<unsigned int, SequenceT> > >
					(subseq));
		}
	} else {
		PrefixSpanDSPCA pspan (dspcaBasisNonZero, dspcaBasisCount, poptions);

		pspan.PCA (S);

		/* Do PCA projection (different from normal projection):
		 *
		 * Xred = (X-repmat(cmean,size(X,1),1))*comp';
		 */
		if (vm.count ("project")) {
			std::ofstream projout (projectFilename.c_str ());

			for (unsigned int n = 0 ; n < S.size () ; ++n) {
				for (unsigned int c = 0 ; c < pspan.components ; ++c) {
					double f = 0.0;

					for (unsigned int el = 0 ; el < pspan.Cpatterns[c].size () ; ++el) {
						bool contains = (pspan.Cpatterns[c][el].Occurence.find (n) !=
							pspan.Cpatterns[c][el].Occurence.end ());

						f += pspan.Ceval[c][el] * ((contains ? 1.0 : 0.0) -
							((double) pspan.Cpatterns[c][el].Occurence.size ()) /
								((double) S.size ()));
					}

					if (c > 0)
						projout << " ";

					projout << f;
				}
				projout << std::endl;
			}
			projout.close ();
		}

		exit (EXIT_SUCCESS);
	}

	std::cout << "Mined " << subseq.size () << " sequences." << std::endl;

	/* Write output file
	 */
	{
		std::ofstream ssout (outputFilename.c_str ());
		for (std::vector<std::pair<unsigned int, SequenceT> >::iterator
			iv = subseq.begin () ; iv != subseq.end () ; ++iv)
		{
			ssout << iv->second << std::endl;
		}
		ssout.close ();
	}

	/* Dump found sequences
	 */
	if (vm.count ("verbose")) {
		for (std::vector<std::pair<unsigned int, SequenceT> >::iterator
			iv = subseq.begin () ; iv != subseq.end () ; ++iv)
		{
			std::cout << iv->first << " times: " << iv->second << std::endl;
		}
	}

	/* Do projection on subsequences and output file.
	 */
	if (vm.count ("project")) {
		std::ofstream projout (projectFilename.c_str ());

		for (unsigned int n = 0 ; n < S.size () ; ++n) {
			for (std::vector<std::pair<unsigned int, SequenceT> >::iterator
				iv = subseq.begin () ; iv != subseq.end () ; ++iv)
			{
				if (iv != subseq.begin ())
					projout << " ";

				if (iv->second.Occurence.find (n) != iv->second.Occurence.end ())
					projout << "1";
				else
					projout << "0";
			}
			projout << std::endl;
		}
		projout.close ();
	}

	/* Cleanup
	 */
	for (std::vector<const Sequence*>::iterator si = S.begin () ;
		si != S.end () ; ++si)
	{
		delete (*si);
	}
}


