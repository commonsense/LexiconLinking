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

    In contemporary English, the comma is often used to stand in for any pause. To avoid 
    over-filtering, this function restricts the filter to short (20 character or fewer) 
    comma-delinated phrases, such as proper names, that should not be the most meaningful
    component of the sentence.
    """
    return re.sub(r'\([a-zA-Z \/]{1,90}\)\s+','',re.sub(r',[a-zA-Z \/]{1,20},','',sent)).replace(" , "," ")

def start_parse_clean(sent):
    return start_parse(remove_nested_phrases(sent).replace("  "," "))

def save_start_error(kind,sent):
    of = open('start_errors.log','a')
    of.write("%s\t%s\n" % (kind,sent))
    print "Error: %s with sentence %s" % (kind,sent)
    of.close()

def remove_plus(word):
    return word.count("+") == 0 and word or word[0:word.index("+")]

def start_parse(sent):
    """ Takes a sentence, as a string, and parses it using the start parser. 
    Ive found it helps parsing to first remove nested phrases and parenthetical 
    remarks """

    url = "http://start.csail.mit.edu/askstart.cgi"
    data = {'server': 'genesis',\
            'query': sent,\
            'pa': 'parse',\
            'action': 'compute-lf'}
    encoded = urllib.urlencode(data)
    try:
        request = urllib2.Request(url, encoded)
        response = urllib2.urlopen(request)
        page = response.read()
    except:
        save_start_error("URL Request", sent)
        return {}
    if page.count("<PRE>") == 0 or page.count("</PRE>") == 0: 
        save_start_error("Unparsable", sent)
        return {}
    pre_start = page.index("<PRE>")
    pre_end = page.index("</PRE>")
    results = defaultdict(list) 
    middle = page[pre_start+4:pre_end].strip()
    for entry in middle.split("\n"):
        # if the line is non-blank and begins with an open bracket
        if len(entry) > 0 and entry[0] == '[':
            # if the tuple is size 3
            split_entry = entry[1:-1].split()
            if len(split_entry) == 3:
                if '+' in split_entry[0]:
                    left, rel, right = map(lambda x: remove_plus(x), split_entry)
                    results[left].append((rel,right))
            else:
                save_start_error("Malformed Tuple", sent)
        else:
            save_start_error("No tuples", sent)
    # remove the defaultdict stuff
    return dict(results)

if __name__ == "__main__":
    # some tests parses 
    print start_parse_clean("i just keep thinking about what her behavior could be caused by") 
    print remove_nested_phrases("had too much lying around in my wardrobes , so took it all out , and did the process of elimination , i actually found a few forgotten numbers , and did some diy also , then rest of the clothes just popped them in a bag and off to the charity shop")
    print start_parse_clean("we lived in a village ( charles hill a dsfasdf asdf asd fasdfasd) out west almost into namibia.")

