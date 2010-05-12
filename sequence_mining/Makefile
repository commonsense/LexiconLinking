

make:
	python preprocess.py 0
	./pspan inseq.txt  -v -K 50 --length-min 5 --project project.txt
	python interpret.py
