
#include <iostream>
#include <PrefixSpan.h>
#include <PrefixSpanK.h>
#include <Lexicon.cpp>

using namespace SeqMine;

int  main (int argc, char* argv[])
{
	std::vector<unsigned int> el = std::vector<unsigned int> ();

    Lexicon lexicon = Lexicon("../../story_keys.txt");
    // def read in file 
    

	// a(abc)(ac)d(cf)
	SequenceT s10 = SequenceT ();
	el.clear (); el.push_back (1); s10.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (2); el.push_back (3); s10.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (3); s10.appendElement (el);
	el.clear (); el.push_back (2); s10.appendElement (el);
	el.clear (); el.push_back (3); el.push_back (6); s10.appendElement (el);

	// (ad)c(bc)(ae)
	SequenceT s20 = SequenceT ();
	el.clear (); el.push_back (1); el.push_back (4); s20.appendElement (el);
	el.clear (); el.push_back (3); s20.appendElement (el);
	el.clear (); el.push_back (2); el.push_back (3); s20.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (3); s20.appendElement (el);

	// (ef)(ab)(df)cb
	SequenceT s30 = SequenceT ();
	el.clear (); el.push_back (5); el.push_back (6); s30.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (2); s30.appendElement (el);
	el.clear (); el.push_back (4); el.push_back (6); s30.appendElement (el);
	el.clear (); el.push_back (2); s30.appendElement (el);
	el.clear (); el.push_back (3); s30.appendElement (el);

	// eg(af)cbc
	SequenceT s40 = SequenceT ();
	el.clear (); el.push_back (1); s40.appendElement (el);
	el.clear (); el.push_back (2); s40.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (6); s40.appendElement (el);
	el.clear (); el.push_back (3); s40.appendElement (el);
	el.clear (); el.push_back (2); s40.appendElement (el);
	el.clear (); el.push_back (3); s40.appendElement (el);

	std::vector<const Sequence*> S = std::vector<const Sequence*> ();
	S.push_back (&s10);
	S.push_back (&s20);
	S.push_back (&s30);
	S.push_back (&s40);
	//PrefixSpanOptions poptions (minsup, length_min, length_max, maxgap);
    PrefixSpanOptions poptions (2, 2, 10, 2);
	//PrefixSpan pspan = PrefixSpan (poptions);
	//pspan.Mine (S);
	
	
   std::vector<std::pair<unsigned int, SequenceT> > subseq;
    	//PrefixSpanOptions poptions (minsup, length_min, length_max, maxgap);
   PrefixSpanK pspan (2, poptions);
   pspan.SetVerification (true);
   pspan.Mine (S);

        //minsup = (unsigned int) (minsupPercent * ((double) S.size ()) + 0.5);
    	//std::cerr << "Minimum support is " << minsupPercent << ", set to "	<< minsup << std::endl;

       	std::copy (pspan.mostFrequent.begin (), pspan.mostFrequent.end (),
        		std::back_insert_iterator<std::vector<std::pair<unsigned int, SequenceT> > >
        			(subseq));

        std::cout << "Mined " << subseq.size () << " sequences." << std::endl;
    	std::cout << "The minimum support threshold is " << pspan.MinimumSupportThreshold () << std::endl;
    
    
    

		for (std::vector<std::pair<unsigned int, SequenceT> >::iterator
			iv = subseq.begin () ; iv != subseq.end () ; ++iv)
		{
			std::cout << iv->first << " times: " << iv->second << std::endl;
		}


	return 0;
}


