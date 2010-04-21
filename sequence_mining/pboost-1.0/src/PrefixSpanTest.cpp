
#include <iostream>

#include <PrefixSpan.h>

using namespace SeqMine;

int
main (int argc, char* argv[])
{
	std::vector<unsigned int> el = std::vector<unsigned int> ();

	// a(abc)(ac)d(cf)
	SequenceT s10 = SequenceT ();
	el.clear (); el.push_back (1); s10.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (2); el.push_back (3); s10.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (3); s10.appendElement (el);
	el.clear (); el.push_back (4); s10.appendElement (el);
	el.clear (); el.push_back (3); el.push_back (6); s10.appendElement (el);

	// (ad)c(bc)(ae)
	SequenceT s20 = SequenceT ();
	el.clear (); el.push_back (1); el.push_back (4); s20.appendElement (el);
	el.clear (); el.push_back (3); s20.appendElement (el);
	el.clear (); el.push_back (2); el.push_back (3); s20.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (5); s20.appendElement (el);

	// (ef)(ab)(df)cb
	SequenceT s30 = SequenceT ();
	el.clear (); el.push_back (5); el.push_back (6); s30.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (2); s30.appendElement (el);
	el.clear (); el.push_back (4); el.push_back (6); s30.appendElement (el);
	el.clear (); el.push_back (3); s30.appendElement (el);
	el.clear (); el.push_back (2); s30.appendElement (el);

	// eg(af)cbc
	SequenceT s40 = SequenceT ();
	el.clear (); el.push_back (5); s40.appendElement (el);
	el.clear (); el.push_back (7); s40.appendElement (el);
	el.clear (); el.push_back (1); el.push_back (6); s40.appendElement (el);
	el.clear (); el.push_back (3); s40.appendElement (el);
	el.clear (); el.push_back (2); s40.appendElement (el);
	el.clear (); el.push_back (3); s40.appendElement (el);

	std::vector<const Sequence*> S = std::vector<const Sequence*> ();
	S.push_back (&s10);
	S.push_back (&s20);
	S.push_back (&s30);
	S.push_back (&s40);

	PrefixSpan pspan = PrefixSpan (2, 0);
	pspan.Mine (S);
}


