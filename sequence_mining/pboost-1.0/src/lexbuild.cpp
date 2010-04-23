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
	unsigned int minsup = 0;
	unsigned int length_min, length_max;
	int maxgap;

	po::options_description config ("Parameters");
	config.add_options ()
		("output", po::value<std::string>(&outputFilename)->default_value ("output.txt"),
			"Frequent subsequence output file")
		("project,p", po::value<std::string>(&projectFilename),
			"Project samples on found sequences.")
		("lexicon,L", "Lexicon Construction Mode")
		("verbose,v", "Verbose mode, dumping found subsequences.")
        ("verify,V", "Verify all projections explictly (slow).");
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
		std::cerr << "pspan, $Id: pspan.cpp + lexicon builder  $" << std::endl;
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

    /* Lexicon Contruction Code 
        Dustin Smith 2010-04-24 
    */
    std::cout << "Constructing a lexicon...\n";
    int it = 1;
    
    unsigned int sidx = 0;
    unsigned int max_len = 0;
	for (std::vector<const Sequence*>::iterator iv = S.begin () ; iv != S.end () ; ++iv, ++sidx)
	    for (unsigned int n = 0 ; n < (*iv)->lengthElements () ; ++n) 
            if ((*iv)->lengthSet(n) > max_len) 
                max_len =(*iv)->lengthSet(n);
                
    /* parameters */
    int lowestLength = 3;
    int topK = 10;
    length_min = max_len-3;
    length_max = 0;
    maxgap = 0;
	
	            
    while(1) {
        
        std::cout << "Lexicon construction Iteration " << it << " with length " << length_min <<  endl;
        
        /* min support percentage mode */
    
        std::vector<std::pair<unsigned int, SequenceT> > subseq;
    	PrefixSpanOptions poptions (minsup, length_min, length_max, maxgap);
    	PrefixSpanK pspan (topK, poptions);
    	if (vm.count ("verify"))
    		pspan.SetVerification (true);

    	pspan.Mine (S);

        //minsup = (unsigned int) (minsupPercent * ((double) S.size ()) + 0.5);
    	//std::cerr << "Minimum support is " << minsupPercent << ", set to "	<< minsup << std::endl;

       	std::copy (pspan.mostFrequent.begin (), pspan.mostFrequent.end (),
        		std::back_insert_iterator<std::vector<std::pair<unsigned int, SequenceT> > >
        			(subseq));

        std::cout << "Mined " << subseq.size () << " sequences." << std::endl;
    	std::cout << "The minimum support threshold for K=" << topK	<< " is " << pspan.MinimumSupportThreshold () << std::endl;

        if (subseq.size() == 0) {
            
            std::cout << "No sequences, lowering length to " << length_min-1;
            length_min = length_min -1;

            
        } else {
        
 
					
        	std::ofstream ssout (outputFilename.c_str ());
        	for (std::vector<std::pair<unsigned int, SequenceT> >::iterator
        		iv = subseq.begin () ; iv != subseq.end () ; ++iv)
        	{
        		ssout << iv->second << std::endl;
        	}
        	ssout.close ();

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
        }
        
        it ++;
        if(lowestLength >= length_min) {
            break; // exit loop
        } else if (it > 50) {
            break;
        }
        
    } // end while
  
   /* Cleanup
     */

  
    for (std::vector<const Sequence*>::iterator si = S.begin () ;
    	si != S.end () ; ++si)
    {
    	delete (*si);
    }
  
  
}


