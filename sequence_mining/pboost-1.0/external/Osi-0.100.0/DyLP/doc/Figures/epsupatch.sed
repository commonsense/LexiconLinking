# Patch file for .eps files produced by IslandDraw 2.3 for SunOS, running on
# Solaris in compatibility mode. We get garbage at the beginning and at the
# end, which has to be stripped. Note that this file must be run with the
# -n flag for sed --- it works by not printing the junk lines.
#							-- lh, 97.3.2 --

# Conveniently enough, this same set of actions does a great job of deleting
# the binary header and preview image of an EPSF file with TIFF preview.
#							-- lh, 98.12.10 --

# For IslandDraw 2.3, the first line(s) are
# "<<binary junk>>%!PS-Adobe-2.0 EPSF-2.0"
# We want to strip the binary junk, leaving the proper first line. For
# IslandDraw 4.10, there won't be any binary junk, but that's ok.

1,/%!PS-Adobe/s/^[^%]*%!PS-Adobe/%!PS-Adobe/p

# God only knows what this line is doing in the file, but it's not needed.

/userdict \/#copies [0-9][0-9]* put/d

# Except for the previous line, we want lines 2 through the final restore.
# After the final line ("restore"), there's yet more binary junk. It might be
# safe to search for "restore" on a line by itself, but it strikes me as more
# robust to look for "%%Trailer" and then print the rest of the lines down
# to the final restore. >> Order is important here! <<  The trick is to use
# %%Trailer twice, deleting one of them. You can't delete the line before
# printing it, eh?

/%%Trailer/,/restore$/p

2,/%%Trailer/{
/%%Trailer/d
p
}

