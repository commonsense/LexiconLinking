#!/usr/bin/env python
"""  This is a phython-based interface to the web-based version of the start parser.

Currently, it is initialized to use the genesis server with the action of computing
the 'lf'.
"""

import urllib2, urllib, re

from collections import defaultdict

def remove_nested_phrases(sent):
    """ removes parenthetical and comma seperated remarks.  Parenthetical remarks
    can be much larger (to be removed) than comma separated statements.

    In contemporary English, the comma is used to stand for "any pause".  Consequently,
    this function restricts the filter to short (20 character or fewer) comma-delinated
    phrases, such as proper names, that should not be the most meaningful component of
    the sentence.
    """
    return re.sub(r'\([a-zA-Z \/]{1,90}\)\s+','',re.sub(r',[a-zA-Z \/]{1,20},','',sent)).replace(" , "," ")

def start_parse_clean(sent):
    return start_parse(remove_nested_phrases(sent).replace("  "," "))

def start_parse(sent):
    """ Takes a sentence, as a string, and parses it using the start parser. 
    Ive found it helps parsing to first remove nested phrases and parenthetical 
    remarks """

    url = "http://start.csail.mit.edu/askstart.cgi"
    data = {'server': 'genesis',\
            'query': sent,\
            'pa': 'parse',\
            'action': 'compute-lf'}
    print "Querying", sent
    encoded = urllib.urlencode(data)
    request = urllib2.Request(url, encoded)
    response = urllib2.urlopen(request)
    page = response.read()
    if page.count("<PRE>") == 0 or page.count("</PRE>") == 0: 
        print "Error: non-parseable sentence: ", sent
        return {}
    pre_start = page.index("<PRE>")
    pre_end = page.index("</PRE>")
    results = defaultdict(list) 
    middle = page[pre_start+4:pre_end].strip()
    for entry in middle.split("\n"):
        # if the line is non-blank and begins with an open bracket
        if len(entry) > 0 and entry[0] == '[':
            # if the tuple is size 3
            if len(entry.split()) == 3:
                left, rel, right = entry[1:-1].split()
                results[left].append((rel,right))
            else:
                print "Error: strange tuple:", entry
    print results
    # remove the defaultdict stuff
    return dict(results)

print remove_nested_phrases("I, who want to live, went to (home) school")

print start_parse_clean("we lived in a village ( charles hill a dsfasdf asdf asd fasdfasd) out west almost into namibia.")

#print start_parse("I love to eat lots of foods")
